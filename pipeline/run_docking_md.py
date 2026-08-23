#!/usr/bin/env python3
"""
Unified AlphaFold Docking & xTB Molecular Dynamics Pipeline
===========================================================
This script compiles the parallel AlphaFold 3 execution and the 
parallel xTB geometry optimization & metrics calculation into a single script.

Workflow:
  1. Parse batch inputs and perform parallel folding search (database searches).
  2. Perform parallel AlphaFold 3 inference across multiple GPUs (1, 2, 3).
  3. Set up ORCA/xTB directories for all docked complex CIF files.
  4. Perform parallel xTB geometry optimizations (loose optimization under CPCM water solvent).
  5. Calculate interface residue RMSDs and compile all docking confidences & energies.

Usage:
  python run_docking_md_pipeline.py --run-all
"""

import os
import sys
import glob
import subprocess
import time
import re
import json
import csv
import shutil
import argparse
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

# Setup environment path to point to AF3 miniconda bin
AF3_ENV_BIN = "/verona/nobackup/rojashrj/miniconda3/envs/af3_312/bin"
os.environ["PATH"] = AF3_ENV_BIN + ":" + os.environ.get("PATH", "")

# Base paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(BASE_DIR)

# ------------------------------------------------------------------------
# Helper Functions for Alignment and Geometry
# ------------------------------------------------------------------------
def parse_pdb(pdb_path):
    coords = []
    atom_types = []
    if not os.path.exists(pdb_path):
        return None, None
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                serial = int(line[6:11].strip())
                atom_name = line[12:16].strip()
                resname = line[17:20].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coords.append([x, y, z])
                atom_types.append({
                    "serial": serial,
                    "atom_name": atom_name,
                    "resname": resname
                })
    return np.array(coords), atom_types

def read_xyz_coords(xyz_path):
    coords = []
    if not os.path.exists(xyz_path):
        return None
    with open(xyz_path, 'r') as f:
        lines = f.readlines()
        if not lines:
            return None
        for line in lines[2:]:
            parts = line.split()
            if len(parts) >= 4:
                x = float(parts[1])
                y = float(parts[2])
                z = float(parts[3])
                coords.append([x, y, z])
    return np.array(coords)

def kabsch_alignment(P, Q):
    c_P = np.mean(P, axis=0)
    c_Q = np.mean(Q, axis=0)
    P_c = P - c_P
    Q_c = Q - c_Q
    H = np.dot(P_c.T, Q_c)
    U, S, Vt = np.linalg.svd(H)
    R = np.dot(Vt.T, U.T)
    if np.linalg.det(R) < 0:
        Vt[2, :] *= -1
        R = np.dot(Vt.T, U.T)
    return c_P, c_Q, R

def apply_transform(coords, c_P, c_Q, R):
    return np.dot(coords - c_P, R.T) + c_Q

def calculate_rmsd_val(coords_a, coords_b):
    return np.sqrt(np.mean(np.sum((coords_a - coords_b) ** 2, axis=1)))

# ------------------------------------------------------------------------
# Step 1 & 2: AlphaFold 3 folding & docking
# ------------------------------------------------------------------------
def clean_job_for_af3_dialect(job):
    job_copy = json.loads(json.dumps(job))
    job_copy['dialect'] = "alphafold3"
    job_copy['version'] = 1
    for seq in job_copy['sequences']:
        if 'protein' in seq:
            if 'count' in seq['protein']:
                del seq['protein']['count']
        elif 'ligand' in seq:
            if 'count' in seq['ligand']:
                del seq['ligand']['count']
    return job_copy

def run_alphafold_pipeline(batch_json, output_dir, gpus):
    af3_dir = os.path.join(PARENT_DIR, "alphafold3")
    orig_cwd = os.getcwd()
    
    # Change directory to alphafold3/
    os.chdir(af3_dir)
    print(f"Working in alphafold3 directory: {os.getcwd()}")
    
    input_json_path = os.path.basename(batch_json)
    precompute_in_dir = "precompute_input_batch"
    precompute_out_dir = "output_precompute"
    
    # Clean previous temporary directories
    gpu_dirs = [f"populated_input_batch_gpu{g.strip()}" for g in gpus.split(",")]
    for d in [precompute_in_dir, precompute_out_dir] + gpu_dirs:
        if os.path.exists(d):
            shutil.rmtree(d)
            
    os.makedirs(precompute_in_dir, exist_ok=True)
    
    # 1. Load batch input
    with open(input_json_path, 'r') as f:
        jobs = json.load(f)
    print(f"Loaded {len(jobs)} jobs from {batch_json}")
    
    # 2. Identify unique protein sequence jobs for precomputation (MSAs and templates search)
    unique_seqs = {}
    for job in jobs:
        for seq in job['sequences']:
            if 'protein' in seq:
                prot_seq = seq['protein']['sequence']
                if prot_seq not in unique_seqs:
                    job_cleaned = clean_job_for_af3_dialect(job)
                    job_cleaned['modelSeeds'] = [1]
                    unique_seqs[prot_seq] = job_cleaned
                break
                
    representative_jobs = list(unique_seqs.values())
    print(f"Identified {len(representative_jobs)} unique protein sequences.")
    
    for idx, job in enumerate(representative_jobs):
        job_file = os.path.join(precompute_in_dir, f"{job['name']}.json")
        with open(job_file, 'w') as f:
            json.dump(job, f, indent=2)
            
    # 3. Run database search (data pipeline only) on unique sequences
    af3_env_python = os.path.join(AF3_ENV_BIN, "python")
    precompute_cmd = [
        af3_env_python,
        "run_alphafold.py",
        f"--input_dir={precompute_in_dir}",
        "--run_inference=False",
        "--model_dir=weights",
        "--db_dir=database",
        f"--output_dir={precompute_out_dir}",
        "--force_output_dir=True"
    ]
    
    print(f"Running AF3 database precomputation: {' '.join(precompute_cmd)}")
    start_time = time.time()
    subprocess.run(precompute_cmd, check=True)
    print(f"Precomputation completed in {time.time() - start_time:.2f} seconds.")
    
    # 4. Extract MSAs and templates
    msa_cache = {}
    for job in representative_jobs:
        job_name = job['name']
        data_json_path = os.path.join(precompute_out_dir, job_name, f"{job_name}_data.json")
        if not os.path.exists(data_json_path):
            raise FileNotFoundError(f"Expected data.json file not found at {data_json_path}")
        with open(data_json_path, 'r') as f:
            data = json.load(f)
        for seq in data['sequences']:
            if 'protein' in seq:
                prot_seq = seq['protein']['sequence']
                msa_cache[prot_seq] = {
                    'unpairedMsa': seq['protein'].get('unpairedMsa'),
                    'pairedMsa': seq['protein'].get('pairedMsa'),
                    'templates': seq['protein'].get('templates', [])
                }
                break
                
    # 5. Populate all jobs with the cached MSAs and templates
    populated_jobs = []
    for job in jobs:
        job_cleaned = clean_job_for_af3_dialect(job)
        for seq in job_cleaned['sequences']:
            if 'protein' in seq:
                prot_seq = seq['protein']['sequence']
                if prot_seq in msa_cache:
                    seq['protein']['unpairedMsa'] = msa_cache[prot_seq]['unpairedMsa']
                    seq['protein']['pairedMsa'] = msa_cache[prot_seq]['pairedMsa']
                    seq['protein']['templates'] = msa_cache[prot_seq]['templates']
                break
        populated_jobs.append(job_cleaned)
        
    # 6. Split populated jobs into GPU chunks and run parallel folding inference
    gpu_list = [g.strip() for g in gpus.split(",")]
    chunk_size = (len(populated_jobs) + len(gpu_list) - 1) // len(gpu_list)
    subprocesses = []
    
    for idx, gpu in enumerate(gpu_list):
        chunk = populated_jobs[idx * chunk_size : (idx + 1) * chunk_size]
        if not chunk:
            continue
        gpu_dir = f"populated_input_batch_gpu{gpu}"
        os.makedirs(gpu_dir, exist_ok=True)
        for job in chunk:
            job_file = os.path.join(gpu_dir, f"{job['name']}.json")
            with open(job_file, 'w') as f:
                json.dump(job, f, indent=2)
                
        run_cmd = [
            af3_env_python,
            "run_alphafold.py",
            f"--input_dir={gpu_dir}",
            "--run_data_pipeline=False",
            "--model_dir=weights",
            "--db_dir=database",
            f"--output_dir={output_dir}",
            "--force_output_dir=True"
        ]
        
        log_path = f"run_af3_gpu{gpu}.log"
        log_file = open(log_path, 'w')
        
        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = gpu
        env["PATH"] = os.environ["PATH"]
        
        print(f"Launching GPU {gpu} worker (logs at {log_path})...")
        p = subprocess.Popen(run_cmd, stdout=log_file, stderr=subprocess.STDOUT, env=env)
        subprocesses.append((p, log_file, gpu))
        
    print("Folding inference processes launched! Waiting for completion...")
    all_completed = False
    while not all_completed:
        time.sleep(30)
        all_completed = True
        for p, log, gpu in subprocesses:
            ret = p.poll()
            if ret is None:
                all_completed = False
            else:
                print(f"GPU {gpu} worker finished with exit code {ret}.")
                
    for p, log, gpu in subprocesses:
        log.close()
    print("Parallel folding inference processes have finished successfully!")
    
    # Restore original directory
    os.chdir(orig_cwd)

# ------------------------------------------------------------------------
# Step 3, 4 & 5: Setup & run xTB Parallel optimizations and compile metrics
# ------------------------------------------------------------------------
def setup_and_run_single(cif_path, target_dir_base, charge, cutoff):
    cif_abspath = os.path.abspath(cif_path)
    job_name = os.path.basename(os.path.dirname(cif_abspath))
    target_dir = os.path.join(PARENT_DIR, target_dir_base, job_name)
    
    xtb_out_path = os.path.join(target_dir, "04_xtb", "cluster_xtb.out")
    if os.path.exists(xtb_out_path):
        try:
            with open(xtb_out_path, 'r') as f:
                content = f.read()
            if "normal termination of xtb" in content:
                print(f"[{job_name}] Already completed successfully. Skipping.")
                return job_name, True, "Already completed"
        except Exception:
            pass

    print(f"[{job_name}] Setting up ORCA/xTB directories...")
    env = os.environ.copy()
    env["PATH"] = AF3_ENV_BIN + ":" + env.get("PATH", "")
    
    setup_cmd = [
        os.path.join(AF3_ENV_BIN, "python"),
        os.path.join(PARENT_DIR, "setup_orca_analysis.py"),
        cif_abspath,
        "-o", target_dir,
        "--charge", str(charge),
        "--cutoff", str(cutoff)
    ]
    
    try:
        subprocess.run(setup_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env)
    except subprocess.CalledProcessError as e:
        print(f"[{job_name}] ERROR: setup_orca_analysis.py failed: {e}")
        return job_name, False, "Setup failed"
        
    pipeline_sh = os.path.join(target_dir, "run_pipeline.sh")
    if not os.path.exists(pipeline_sh):
        print(f"[{job_name}] ERROR: run_pipeline.sh missing")
        return job_name, False, "run_pipeline.sh missing"
        
    env["OMP_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["VECLIB_MAXIMUM_THREADS"] = "1"
    env["NUMEXPR_NUM_THREADS"] = "1"
    
    log_path = os.path.join(target_dir, "pipeline.log")
    print(f"[{job_name}] Running xTB pipeline (logs at {log_path})...")
    
    start_time = time.time()
    try:
        with open(log_path, "w") as log_file:
            subprocess.run(["./run_pipeline.sh"], cwd=target_dir, stdout=log_file, stderr=subprocess.STDOUT, check=True, env=env)
        elapsed = time.time() - start_time
        print(f"[{job_name}] SUCCESS: xTB pipeline completed in {elapsed:.1f}s")
        return job_name, True, f"Completed in {elapsed:.1f}s"
    except subprocess.CalledProcessError as e:
        print(f"[{job_name}] ERROR: xTB pipeline failed. Check logs at {log_path}")
        return job_name, False, "xTB run failed"

def compile_metrics(cif_files, batch_json, target_dir_base, output_summary_csv):
    batch_json_path = os.path.join(PARENT_DIR, "alphafold3", os.path.basename(batch_json))
    job_map = {}
    if os.path.exists(batch_json_path):
        try:
            with open(batch_json_path, 'r') as f:
                jobs = json.load(f)
            for j in jobs:
                raw_name = j['name']
                sanitized_name = re.sub(r'[^a-zA-Z0-9_]', '_', raw_name)
                sanitized_name = re.sub(r'_+', '_', sanitized_name).strip('_')
                
                dataset = None
                for d in ["glaucoma_bulk", "glaucoma", "cardiomyopathy", "parkinsons", "hcc"]:
                    if d in sanitized_name:
                        dataset = d
                        break
                
                if dataset:
                    parts = sanitized_name.split('_')
                    idx = parts[-1]
                    suffix = f"_{dataset}_{idx}"
                    prefix = "docking_"
                    method = sanitized_name
                    if method.startswith(prefix):
                        method = method[len(prefix):]
                    if method.endswith(suffix):
                        method = method[:-len(suffix)]
                else:
                    method = "Unknown"
                    dataset = "Unknown"
                
                smiles = ""
                for seq in j['sequences']:
                    if 'ligand' in seq:
                        smiles = seq['ligand']['smiles']
                        break
                job_map[sanitized_name] = {"Method": method, "Dataset": dataset, "SMILES": smiles}
        except Exception as e:
            print(f"Warning: Could not parse batch JSON for metadata mapping: {e}")

    output_rows = []
    
    for cif_file in cif_files:
        job_name = os.path.basename(os.path.dirname(cif_file))
        target_dir = os.path.join(PARENT_DIR, target_dir_base, job_name)
        
        # Parse AF3 summary confidences
        af3_summary_json = os.path.join(PARENT_DIR, "alphafold3", "output_result5", job_name, f"{job_name}_summary_confidences.json")
        iptm = None
        ptm = None
        ranking_score = None
        
        if os.path.exists(af3_summary_json):
            try:
                with open(af3_summary_json, 'r') as f:
                    conf = json.load(f)
                iptm = conf.get("iptm")
                ptm = conf.get("ptm")
                ranking_score = conf.get("ranking_score")
            except Exception:
                pass
                
        # Parse xTB output
        xtb_out_path = os.path.join(target_dir, "04_xtb", "cluster_xtb.out")
        xtb_energy = None
        xtb_gap = None
        xtb_dipole = None
        
        if os.path.exists(xtb_out_path):
            try:
                with open(xtb_out_path, 'r') as f:
                    content = f.read()
                
                energy_match = re.search(r'TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh', content)
                if energy_match:
                    xtb_energy = float(energy_match.group(1))
                    
                gap_match = re.search(r'HOMO-LUMO GAP\s+(\d+\.\d+)\s+eV', content)
                if gap_match:
                    xtb_gap = float(gap_match.group(1))
                    
                dipole_pos = content.rfind("molecular dipole:")
                if dipole_pos != -1:
                    dipole_block = content[dipole_pos:dipole_pos+500]
                    full_line_match = re.search(r"full:\s+[\d.-]+\s+[\d.-]+\s+[\d.-]+\s+([\d.-]+)", dipole_block)
                    if full_line_match:
                        xtb_dipole = float(full_line_match.group(1))
            except Exception:
                pass
                
        # Calculate RMSDs
        pocket_bb_rmsd = None
        ligand_pocket_aligned_rmsd = None
        ligand_self_aligned_rmsd = None
        
        pdb_path = os.path.join(target_dir, "03_cluster", "cluster.pdb")
        xyz_path = os.path.join(target_dir, "04_xtb", "xtbopt.xyz")
        
        if os.path.exists(pdb_path) and os.path.exists(xyz_path):
            try:
                start_coords, atom_types = parse_pdb(pdb_path)
                opt_coords = read_xyz_coords(xyz_path)
                
                if start_coords is not None and opt_coords is not None and len(start_coords) == len(opt_coords):
                    ligand_indices = [i for i, at in enumerate(atom_types) if at["resname"] == "LIG"]
                    pocket_bb_indices = [
                        i for i, at in enumerate(atom_types)
                        if at["resname"] != "LIG" and at["atom_name"] in ["CA", "C", "N", "O"]
                    ]
                    
                    if pocket_bb_indices:
                        P_start = start_coords[pocket_bb_indices]
                        P_opt = opt_coords[pocket_bb_indices]
                        c_P, c_Q, R = kabsch_alignment(P_opt, P_start)
                        P_opt_aligned = apply_transform(P_opt, c_P, c_Q, R)
                        pocket_bb_rmsd = calculate_rmsd_val(P_opt_aligned, P_start)
                        
                        if ligand_indices:
                            L_start = start_coords[ligand_indices]
                            L_opt = opt_coords[ligand_indices]
                            L_opt_aligned = apply_transform(L_opt, c_P, c_Q, R)
                            ligand_pocket_aligned_rmsd = calculate_rmsd_val(L_opt_aligned, L_start)
                            
                    if ligand_indices:
                        L_start = start_coords[ligand_indices]
                        L_opt = opt_coords[ligand_indices]
                        c_L_opt, c_L_start, R_L = kabsch_alignment(L_opt, L_start)
                        L_opt_self_aligned = apply_transform(L_opt, c_L_opt, c_L_start, R_L)
                        ligand_self_aligned_rmsd = calculate_rmsd_val(L_opt_self_aligned, L_start)
            except Exception as e:
                print(f"Warning: Failed to calculate RMSD for {job_name}: {e}")
                
        meta = job_map.get(job_name, {"Method": "Unknown", "Dataset": "Unknown", "SMILES": "Unknown"})
        
        output_rows.append({
            "Job_Name": job_name,
            "SMILES": meta["SMILES"],
            "Method": meta["Method"],
            "test_dataset": meta["Dataset"],
            "iPTM": iptm,
            "PTM": ptm,
            "Ranking_Score": ranking_score,
            "xTB_Energy_Eh": xtb_energy,
            "xTB_Gap_eV": xtb_gap,
            "xTB_Dipole_Debye": xtb_dipole,
            "Pocket_Backbone_RMSD": pocket_bb_rmsd,
            "Ligand_RMSD_Pocket_Aligned": ligand_pocket_aligned_rmsd,
            "Ligand_Self_Aligned_RMSD": ligand_self_aligned_rmsd
        })
        
    with open(output_summary_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            "Job_Name", "SMILES", "Method", "test_dataset", 
            "iPTM", "PTM", "Ranking_Score", "xTB_Energy_Eh", "xTB_Gap_eV", "xTB_Dipole_Debye",
            "Pocket_Backbone_RMSD", "Ligand_RMSD_Pocket_Aligned", "Ligand_Self_Aligned_RMSD"
        ])
        writer.writeheader()
        writer.writerows(output_rows)
        
    print(f"Metrics written successfully to {os.path.abspath(output_summary_csv)}")

def main():
    parser = argparse.ArgumentParser(description="Unified AF3 Docking & xTB Optimization Pipeline")
    parser.add_argument("--batch_json", default="af3_batch_input_result5.json", help="Path to batch JSON input for AF3")
    parser.add_argument("--af3_out", default="output_result5", help="AlphaFold 3 output directory")
    parser.add_argument("--xtb_out", default="run_docking_result5", help="xTB runs base output directory")
    parser.add_argument("--gpus", default="1,2,3", help="CUDA GPUs list for parallel AF3 runs (comma-separated)")
    parser.add_argument("--workers", type=int, default=64, help="Number of parallel workers for xTB optimizations")
    parser.add_argument("--cutoff", type=float, default=4.0, help="Interface distance cutoff for cluster extraction")
    parser.add_argument("--charge", type=int, default=0, help="Total charge for cluster")
    parser.add_argument("--summary", default="docking_summary.csv", help="Output summary CSV name")
    parser.add_argument("--af3_env_bin", default="/verona/nobackup/rojashrj/miniconda3/envs/af3_312/bin", help="Path to AlphaFold 3 environment bin")
    parser.add_argument("--run-all", action="store_true", help="Execute the entire pipeline")
    
    args = parser.parse_args()
    
    global AF3_ENV_BIN
    AF3_ENV_BIN = args.af3_env_bin
    os.environ["PATH"] = AF3_ENV_BIN + ":" + os.environ.get("PATH", "")
    
    if not args.run_all:
        parser.print_help()
        sys.exit(1)
        
    # Step 1 & 2: AlphaFold 3 parallel folding/docking
    print("==================================================")
    print("1. RUNNING PARALLEL ALPHAFOLD 3 INFERENCE")
    print("==================================================")
    run_alphafold_pipeline(args.batch_json, args.af3_out, args.gpus)
    
    # Check outputs
    cif_pattern = os.path.join(PARENT_DIR, "alphafold3", args.af3_out, "*/*.cif")
    cif_files = [f for f in glob.glob(cif_pattern) if "seed-" not in f]
    print(f"Found {len(cif_files)} top-ranked mmCIF structures to process.")
    
    if not cif_files:
        print("Error: No mmCIF files found. Folding pipeline may have failed.")
        sys.exit(1)
        
    # Step 3 & 4: Setup & run xTB parallel optimization
    print("\n==================================================")
    print("2. EXECUTING PARALLEL XTB OPTIMIZATIONS")
    print("==================================================")
    results_summary = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(setup_and_run_single, cif_file, args.xtb_out, args.charge, args.cutoff): cif_file for cif_file in cif_files}
        for future in as_completed(futures):
            job_name, success, message = future.result()
            results_summary.append((job_name, success, message))
            
    success_count = sum(1 for _, s, _ in results_summary if s)
    print(f"\nxTB Optimizations completed! {success_count}/{len(cif_files)} succeeded.")
    
    # Step 5: Compile metrics
    print("\n==================================================")
    print("3. COMPILING ALIGNMENT AND ENERGY METRICS")
    print("==================================================")
    summary_path = os.path.join(PARENT_DIR, args.summary)
    compile_metrics(cif_files, args.batch_json, args.xtb_out, summary_path)
    print("Docking & MD Minimization Pipeline Completed Successfully!")

if __name__ == "__main__":
    main()
