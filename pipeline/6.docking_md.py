#!/usr/bin/env python3
"""
Step 6: AlphaFold 3 Docking & Semi-Empirical xTB Quantum Mechanical (QM) Analysis
Standalone pipeline for:
  1. AlphaFold 3 parallel batch inference (accelerated MSA caching + multi-GPU splitting).
  2. Active site pocket interface extraction (take_interface.py).
  3. Semi-empirical GFN2-xTB calculation: Total SCF Energy, HOMO-LUMO Gap, Dipole Moment, WBO.
  4. Full metric extraction & parsing:
     - AlphaFold 3 confidence (Ranking_Score, iPTM, PTM, mean_protein_plddt, mean_ligand_plddt)
     - Pocket-aligned ligand RMSD (Ligand_RMSD_Pocket_Aligned, Pocket_Backbone_RMSD, Ligand_Self_Aligned_RMSD)
     - Semi-empirical / QM descriptors (HOMO_LUMO_Gap_eV, Total_Dipole_Moment_Debye, Total_SCF_Energy_Eh, Min_WBO)
     - Metadata mapping (Method, Model_Index, test_dataset, SMILES)

All outputs default to TheraFlow/outputs/ and TheraFlow/pipeline/.
"""

import os
import sys
import re
import json
import time
import glob
import shutil
import argparse
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm


# =============================================================================
# 1. HELPER FUNCTIONS & DIALECT HANDLING
# =============================================================================

def clean_job_for_af3_dialect(job):
    """Ensures AF3 job dialect conformity (deletes forbidden 'count' keys, adds dialect/version)."""
    job_copy = json.loads(json.dumps(job))
    job_copy['dialect'] = "alphafold3"
    job_copy['version'] = 1
    
    for seq in job_copy.get('sequences', []):
        if 'protein' in seq:
            if 'count' in seq['protein']:
                del seq['protein']['count']
        elif 'ligand' in seq:
            if 'count' in seq['ligand']:
                del seq['ligand']['count']
    return job_copy


def run_af3_batch_pipeline(
    input_json_path="TheraFlow/outputs/af3_batch_input.json",
    af3_dir="/users/stgrad/rojashrj/Projects/omics2mol/alphafold3",
    output_dir="TheraFlow/outputs/af3_predictions",
    af3_python="/verona/nobackup/rojashrj/miniconda3/envs/af3_312/bin/python",
    gpus=["0", "1", "2", "3"],
    precompute_dir="/users/stgrad/rojashrj/Projects/omics2mol/alphafold3/output_precompute"
):
    """Executes AlphaFold 3 batch neural folding across multiple GPUs with cached MSAs."""
    input_json = Path(input_json_path).resolve()
    if not input_json.exists():
        raise FileNotFoundError(f"Input JSON not found: {input_json}")
        
    out_dir = Path(output_dir).resolve()
    os.makedirs(out_dir, exist_ok=True)
    
    with open(input_json, "r") as f:
        jobs = json.load(f)
    print(f"Loaded {len(jobs)} AF3 docking jobs from {input_json}")
    
    # 1. Load cached MSAs and templates from precompute directory
    msa_cache = {}
    precompute_path = Path(precompute_dir).resolve()
    if precompute_path.exists():
        for sub in os.listdir(precompute_path):
            sub_dir = precompute_path / sub
            if sub_dir.is_dir():
                for f_name in os.listdir(sub_dir):
                    if f_name.endswith("_data.json"):
                        data_json_path = sub_dir / f_name
                        try:
                            with open(data_json_path, "r") as df:
                                data = json.load(df)
                            for seq in data.get("sequences", []):
                                if "protein" in seq:
                                    prot_seq = seq["protein"]["sequence"]
                                    msa_cache[prot_seq] = {
                                        "unpairedMsa": seq["protein"].get("unpairedMsa"),
                                        "pairedMsa": seq["protein"].get("pairedMsa"),
                                        "templates": seq["protein"].get("templates", [])
                                    }
                                    break
                        except Exception as e:
                            print(f"Warning loading {data_json_path}: {e}")
    print(f"Cached MSAs/templates loaded for {len(msa_cache)} target protein sequences.")
    
    # 2. Populate jobs with cached MSAs
    populated_jobs = []
    for job in jobs:
        job_cleaned = clean_job_for_af3_dialect(job)
        for seq in job_cleaned.get("sequences", []):
            if "protein" in seq:
                prot_seq = seq["protein"]["sequence"]
                if prot_seq in msa_cache:
                    seq["protein"]["unpairedMsa"] = msa_cache[prot_seq]["unpairedMsa"]
                    seq["protein"]["pairedMsa"] = msa_cache[prot_seq]["pairedMsa"]
                    seq["protein"]["templates"] = msa_cache[prot_seq]["templates"]
                break
        populated_jobs.append(job_cleaned)
        
    # 3. Split across GPUs
    n_gpus = len(gpus)
    chunk_size = (len(populated_jobs) + n_gpus - 1) // n_gpus
    processes = []
    
    af3_path = Path(af3_dir).resolve()
    model_dir = af3_path / "weights"
    db_dir = af3_path / "database"
    
    for idx, gpu in enumerate(gpus):
        chunk = populated_jobs[idx * chunk_size : (idx + 1) * chunk_size]
        if not chunk:
            continue
            
        gpu_dir = out_dir / f"populated_input_batch_gpu{gpu}"
        if gpu_dir.exists():
            shutil.rmtree(gpu_dir)
        os.makedirs(gpu_dir, exist_ok=True)
        
        for j in chunk:
            j_file = gpu_dir / f"{j['name']}.json"
            with open(j_file, "w") as jf:
                json.dump(j, jf, indent=2)
                
        cmd = [
            af3_python, str(af3_path / "run_alphafold.py"),
            f"--input_dir={gpu_dir}",
            "--run_data_pipeline=False",
            f"--model_dir={model_dir}",
            f"--db_dir={db_dir}",
            f"--output_dir={out_dir}",
            "--force_output_dir=True"
        ]
        
        log_path = out_dir / f"run_af3_gpu{gpu}.log"
        log_file = open(log_path, "w")
        
        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = str(gpu)
        af3_env_bin = str(Path(af3_python).parent)
        env["PATH"] = af3_env_bin + ":" + env.get("PATH", "")
        
        print(f"Launching GPU {gpu} worker ({len(chunk)} jobs, logs at {log_path})...")
        p = subprocess.Popen(cmd, stdout=log_file, stderr=subprocess.STDOUT, env=env)
        processes.append((gpu, p, log_file))
        
    print(f"All {len(processes)} GPU inference workers launched. Monitoring execution...")
    for gpu, p, log_file in processes:
        p.wait()
        log_file.close()
        print(f"GPU {gpu} worker finished with exit code {p.returncode}")
        
    print("All AlphaFold 3 batch inference runs complete!")


# =============================================================================
# 2. mmCIF TO PDB CONVERSION
# =============================================================================

def convert_cif_to_pdb(cif_path, pdb_path):
    """Converts AF3 mmCIF into a clean standard PDB file with protein (chain A) and ligand (chain B)."""
    cif_path = Path(cif_path).resolve()
    pdb_path = Path(pdb_path).resolve()
    os.makedirs(pdb_path.parent, exist_ok=True)
    
    records = []
    atom_count = 1
    with open(cif_path, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                parts = line.split()
                if len(parts) >= 15:
                    rec_type = parts[0]
                    name = parts[3]
                    resname = "LIG" if "LIG" in parts[5].upper() else parts[5]
                    chain = "B" if "LIG" in resname else (parts[6] if len(parts[6]) == 1 else "A")
                    resseq = "1" if "LIG" in resname else (parts[8] if parts[8].isdigit() else "1")
                    try:
                        x, y, z = float(parts[10]), float(parts[11]), float(parts[12])
                        elem = parts[2] if len(parts[2]) <= 2 else parts[3][0]
                        record = f"{rec_type:<6}{atom_count:>5} {name:<4} {resname:<3} {chain}{int(resseq):>4}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 50.00          {elem:>2}\n"
                        records.append(record)
                        atom_count += 1
                    except Exception:
                        pass
    with open(pdb_path, "w") as f:
        f.writelines(records)
        f.write("END\n")
    return str(pdb_path)


# =============================================================================
# 3. POCKET EXTRACTION & xTB CALCULATION PER COMPLEX
# =============================================================================

def run_single_xtb_task(cif_file, qm_dir, cutoff=4.0):
    """Executes active site extraction and GFN2-xTB calculation for a single CIF."""
    cif_path = Path(cif_file).resolve()
    cif_name = cif_path.stem
    job_name = cif_path.parent.name
    
    target_out = Path(qm_dir) / f"run_{cif_name}"
    for d in ["01_prepare", "03_cluster", "04_xtb"]:
        os.makedirs(target_out / d, exist_ok=True)
        
    pdb_path = target_out / "01_prepare" / "complex_clean.pdb"
    convert_cif_to_pdb(cif_path, pdb_path)
    
    # Run pdb2pqr (fallback to clean pdb)
    pqr_path = target_out / "01_prepare" / "complex_parse.pqr"
    pdb2pqr_bin = "/verona/nobackup/rojashrj/miniconda3/envs/af3_312/bin/pdb2pqr"
    
    if os.path.exists(pdb2pqr_bin) and pdb_path.exists():
        try:
            cmd = [pdb2pqr_bin, "--ff=PARSE", "--titration-state-method=propka", "--with-ph=7.4", "complex_clean.pdb", "complex_parse.pqr"]
            subprocess.run(cmd, cwd=str(target_out / "01_prepare"), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=45)
            # Append HETATMs
            with open(pdb_path, "r") as pf:
                hetatms = [l for l in pf if l.startswith("HETATM")]
            with open(pqr_path, "a") as pqf:
                pqf.writelines(hetatms)
        except Exception:
            shutil.copy2(pdb_path, pqr_path)
    else:
        if pdb_path.exists():
            shutil.copy2(pdb_path, pqr_path)
            
    # Run take_interface.py
    pipeline_dir = Path(__file__).resolve().parent
    take_iface = (pipeline_dir / "take_interface.py").resolve()
    
    cmd_iface = [
        sys.executable, str(take_iface),
        str(pqr_path.resolve()), "LIG", str(cutoff), "cluster", "0"
    ]
    try:
        subprocess.run(cmd_iface, cwd=str(target_out / "03_cluster"), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    except Exception as e:
        print(f"Warning: take_interface failed for {cif_name}: {e}")
        return str(target_out)
        
    # Run xTB calculation
    xtb_bin = "/verona/nobackup/rojashrj/miniconda3/envs/af3_312/bin/xtb"
    cluster_xyz = (target_out / "03_cluster" / "cluster.xyz").resolve()
    charge_file = target_out / "03_cluster" / "cluster_charge.txt"
    
    chrg = "0"
    if charge_file.exists():
        try:
            with open(charge_file, "r") as f:
                chrg = f.read().strip() or "0"
        except Exception:
            chrg = "0"
            
    if cluster_xyz.exists() and os.path.exists(xtb_bin):
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = "4"
        env["MKL_NUM_THREADS"] = "4"
        out_log = target_out / "04_xtb" / "cluster_xtb.out"
        cmd_xtb = [
            xtb_bin, str(cluster_xyz),
            "--gfn", "2",
            "--opt", "loose",
            "--alpb", "water",
            "--chrg", chrg,
            "--wbo",
            "--dipole"
        ]
        try:
            with open(out_log, "w") as out_f:
                subprocess.run(cmd_xtb, cwd=str(target_out / "04_xtb"), stdout=out_f, stderr=subprocess.STDOUT, env=env, timeout=120)
        except Exception as e:
            print(f"Warning: xTB loose opt timed out or failed for {cif_name}: {e}")
            
        # If dipole or energy was not finalized, run fast property calculation
        need_sp = True
        if out_log.exists():
            txt = out_log.read_text()
            if "molecular dipole:" in txt and ("TOTAL ENERGY" in txt or "total energy" in txt):
                need_sp = False
        if need_sp:
            cmd_sp = [
                xtb_bin, str(cluster_xyz),
                "--gfn", "2",
                "--alpb", "water",
                "--chrg", chrg,
                "--wbo",
                "--dipole"
            ]
            try:
                with open(out_log, "a") as out_f:
                    subprocess.run(cmd_sp, cwd=str(target_out / "04_xtb"), stdout=out_f, stderr=subprocess.STDOUT, env=env, timeout=60)
            except Exception as e:
                print(f"Warning: xTB single-point fallback failed for {cif_name}: {e}")

        # Copy cluster.xyz to xtbopt.xyz if opt didn't generate it
        if not (target_out / "04_xtb" / "xtbopt.xyz").exists():
            shutil.copy2(cluster_xyz, target_out / "04_xtb" / "xtbopt.xyz")
            
    return str(target_out)


def run_xtb_batch_pipeline(
    af3_dir="TheraFlow/outputs/af3_predictions",
    qm_dir="TheraFlow/outputs/docking_md_runs",
    n_jobs=4,
    cutoff=4.0
):
    """Runs pocket extraction and xTB calculations across all AF3 predictions in parallel."""
    cif_files = sorted(glob.glob(f"{af3_dir}/docking_*/*_model.cif"))
    if not cif_files:
        cif_files = sorted(glob.glob(f"{af3_dir}/docking_*/*.cif"))
        
    print(f"Found {len(cif_files)} top-model CIF predictions for xTB QM analysis.")
    os.makedirs(qm_dir, exist_ok=True)
    
    print(f"Executing parallel active site extraction & GFN2-xTB calculations ({n_jobs} workers)...")
    results = Parallel(n_jobs=n_jobs)(
        delayed(run_single_xtb_task)(c, qm_dir, cutoff) for c in tqdm(cif_files, desc="xTB QM Processing")
    )
    print(f"Completed xTB QM calculations for {len(results)} complexes!")
    return results


# =============================================================================
# 4. KABSCH RMSD & QM METRIC PARSING
# =============================================================================

def parse_pdb_atoms(pdb_path):
    """Parses coordinate points and residue metadata from PDB file."""
    coords = []
    atom_types = []
    if not os.path.exists(pdb_path):
        return None, None
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    atom_name = line[12:16].strip()
                    resname = line[17:20].strip().upper()
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coords.append([x, y, z])
                    atom_types.append({"atom_name": atom_name, "resname": resname})
                except Exception:
                    pass
    return np.array(coords), atom_types


def read_xyz_coordinates(xyz_path):
    """Reads 3D coordinate array from XYZ file."""
    if not os.path.exists(xyz_path):
        return None
    coords = []
    with open(xyz_path, "r") as f:
        lines = f.readlines()
        if len(lines) < 3:
            return None
        for line in lines[2:]:
            parts = line.split()
            if len(parts) >= 4:
                try:
                    coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
                except Exception:
                    pass
    return np.array(coords)


def kabsch_alignment(P, Q):
    """Computes optimal rotation matrix R and translations via Kabsch algorithm (SVD)."""
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
    """Applies Kabsch translation and rotation to coordinate array."""
    return np.dot(coords - c_P, R.T) + c_Q


def calculate_rmsd(a, b):
    """Calculates root-mean-square deviation between two coordinate matrices."""
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def parse_xtb_output(out_path):
    """Parses Total SCF Energy, HOMO-LUMO Gap, and Dipole Moment from xTB output."""
    if not os.path.exists(out_path):
        return {}
    with open(out_path, "r") as f:
        content = f.read()
    energy_matches = re.findall(r"(?:TOTAL ENERGY|total energy)\s*[:]?\s*([\d.-]+)\s*Eh", content, re.IGNORECASE)
    gap_matches = re.findall(r"HOMO-LUMO\s+(?:GAP|gap)\s*[:]?\s*([\d.-]+)\s*eV", content, re.IGNORECASE)
    dipole_pos = content.rfind("molecular dipole:")
    dipole_moment = None
    if dipole_pos != -1:
        m = re.search(r"full:\s+[\d.-]+\s+[\d.-]+\s+[\d.-]+\s+([\d.-]+)", content[dipole_pos:dipole_pos+500])
        if m:
            dipole_moment = float(m.group(1))
    return {
        "Total_SCF_Energy_Eh": float(energy_matches[-1]) if energy_matches else None,
        "HOMO_LUMO_Gap_eV": float(gap_matches[-1]) if gap_matches else None,
        "Total_Dipole_Moment_Debye": dipole_moment
    }


def parse_wbo_min(wbo_path, ligand_indices, pocket_indices):
    """Parses minimum Wiberg bond order among interacting ligand-pocket atom pairs."""
    if not os.path.exists(wbo_path):
        return None
    wbo_vals = []
    with open(wbo_path, "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) == 3:
                try:
                    a1 = int(parts[0]) - 1
                    a2 = int(parts[1]) - 1
                    val = float(parts[2])
                    if (a1 in ligand_indices and a2 in pocket_indices) or (a2 in ligand_indices and a1 in pocket_indices):
                        wbo_vals.append(val)
                except Exception:
                    pass
    return min(wbo_vals) if wbo_vals else None


def extract_docking_and_qm_metrics(
    af3_dir="alphafold3/output_result5",
    qm_dir="TheraFlow/outputs/docking_md_runs",
    input_json_path="TheraFlow/outputs/af3_batch_input.json",
    output_csv="TheraFlow/outputs/docking_qm_metrics.csv"
):
    """
    Extracts unified AlphaFold 3 Docking & QM metrics:
    Directory, Method, Model_Index, Ranking_Score, Ligand_RMSD_Pocket_Aligned,
    HOMO_LUMO_Gap_eV, Total_Dipole_Moment_Debye, Total_SCF_Energy_Eh,
    iPTM, PTM, Pocket_Backbone_RMSD, Ligand_Self_Aligned_RMSD, Min_WBO,
    mean_protein_plddt, mean_ligand_plddt, SMILES, test_dataset.
    """
    af3_candidates = [Path(af3_dir), Path("alphafold3/output_result5"), Path("../alphafold3/output_result5"), Path("TheraFlow/outputs/af3_predictions")]
    af3_path = next((p.resolve() for p in af3_candidates if p.exists()), Path(af3_dir).resolve())

    qm_candidates = [Path(qm_dir), Path("TheraFlow/outputs/docking_md_runs"), Path("TheraFlow/outputs_new/docking_md_runs"), Path("outputs/docking_md_runs")]
    qm_path = next((p.resolve() for p in qm_candidates if p.exists()), Path(qm_dir).resolve())

    # Build job -> SMILES mapping from all available batch input JSONs
    job_meta_map = {}
    json_candidates = [
        Path(input_json_path),
        Path("TheraFlow/outputs/af3_batch_input.json"),
        Path("TheraFlow/outputs_new/af3_batch_input.json"),
        Path("alphafold3/af3_batch_input_result5.json"),
        Path("alphafold3/af3_batch_input_result5_transgem.json"),
        Path("alphafold3/af3_batch_input_supplement.json"),
    ]
    for j_path in json_candidates:
        if j_path.exists():
            try:
                with open(j_path, "r") as f:
                    jobs = json.load(f)
                for j in jobs:
                    name = j.get("name")
                    lig_smi = ""
                    for s in j.get("sequences", []):
                        if "ligand" in s:
                            lig_smi = s["ligand"].get("smiles", "")
                    if name and lig_smi:
                        job_meta_map[name] = lig_smi
                        job_meta_map[name.lower()] = lig_smi
                        job_meta_map[f"run_{name}_model"] = lig_smi
                        job_meta_map[f"run_{name.lower()}_model"] = lig_smi
            except Exception as e:
                pass

    # -------------------------------------------------------------------------
    # 1. Parse QM runs from qm_path (TheraFlow/outputs/docking_md_runs)
    # -------------------------------------------------------------------------
    qm_records = {}
    if qm_path.exists():
        qm_dirs = [d for d in qm_path.iterdir() if d.is_dir() and not d.name.startswith(("_", "."))]
        print(f"Parsing QM & Kabsch RMSD metrics across {len(qm_dirs)} simulation directories...")
        for qd in qm_dirs:
            raw_name = qd.name
            job_name = raw_name[4:] if raw_name.startswith("run_") else raw_name
            if job_name.endswith("_model"):
                job_name = job_name[:-6]

            method = "Unknown"
            dataset = "Unknown"
            model_idx = None

            for d in ["glaucoma_bulk", "glaucoma", "cardiomyopathy", "parkinsons", "hcc", "lincs1000_1", "lincs1000_2"]:
                if d in raw_name.lower():
                    dataset = d
                    break

            parts = raw_name.split("_")
            for p in parts:
                if p.isdigit():
                    model_idx = int(p)
                    break
            
            # Identify method name
            for m in ["biaee", "bicev", "fame", "gex2sgen", "gexmolgen", "gxrnn", "paccmann_rl", "smilesgen", "transgem", "triomphe"]:
                if m in raw_name.lower():
                    method = m
                    break

            # xTB properties
            xtb_out = qd / "04_xtb" / "cluster_xtb.out"
            xtb_props = parse_xtb_output(str(xtb_out))

            # RMSD & WBO
            pdb_path = qd / "03_cluster" / "cluster.pdb"
            xyz_path = qd / "04_xtb" / "xtbopt.xyz"
            wbo_path = qd / "04_xtb" / "wbo"

            pocket_bb_rmsd = None
            ligand_pocket_aligned_rmsd = None
            ligand_self_aligned_rmsd = None
            min_wbo = None

            if pdb_path.exists() and xyz_path.exists():
                try:
                    start_coords, atom_types = parse_pdb_atoms(str(pdb_path))
                    opt_coords = read_xyz_coordinates(str(xyz_path))
                    if start_coords is not None and opt_coords is not None and len(start_coords) == len(opt_coords):
                        lig_idx = [i for i, a in enumerate(atom_types) if a["resname"] == "LIG"]
                        pock_idx = [i for i, a in enumerate(atom_types) if a["resname"] != "LIG"]
                        pock_bb_idx = [i for i, a in enumerate(atom_types) if a["resname"] != "LIG" and a["atom_name"] in ["CA", "C", "N", "O"]]

                        if pock_bb_idx:
                            P_start = start_coords[pock_bb_idx]
                            P_opt = opt_coords[pock_bb_idx]
                            c_P, c_Q, R = kabsch_alignment(P_opt, P_start)
                            P_opt_aligned = apply_transform(P_opt, c_P, c_Q, R)
                            pocket_bb_rmsd = round(calculate_rmsd(P_opt_aligned, P_start), 3)

                            if lig_idx:
                                L_start = start_coords[lig_idx]
                                L_opt = opt_coords[lig_idx]
                                L_opt_aligned = apply_transform(L_opt, c_P, c_Q, R)
                                ligand_pocket_aligned_rmsd = round(calculate_rmsd(L_opt_aligned, L_start), 3)

                        if lig_idx:
                            L_start = start_coords[lig_idx]
                            L_opt = opt_coords[lig_idx]
                            c_L_opt, c_L_start, R_L = kabsch_alignment(L_opt, L_start)
                            L_opt_self_aligned = apply_transform(L_opt, c_L_opt, c_L_start, R_L)
                            ligand_self_aligned_rmsd = round(calculate_rmsd(L_opt_self_aligned, L_start), 3)

                        min_wbo = parse_wbo_min(str(wbo_path), lig_idx, pock_idx)
                        if min_wbo is not None:
                            min_wbo = round(min_wbo, 4)
                except Exception as e:
                    pass

            smi_val = job_meta_map.get(raw_name, job_meta_map.get(job_name, job_meta_map.get(f"docking_{job_name}", "Unknown")))

            qm_records[raw_name] = {
                "Directory": raw_name,
                "Method": method,
                "Model_Index": model_idx,
                "Ligand_RMSD_Pocket_Aligned": ligand_pocket_aligned_rmsd,
                "HOMO_LUMO_Gap_eV": xtb_props.get("HOMO_LUMO_Gap_eV"),
                "Total_Dipole_Moment_Debye": xtb_props.get("Total_Dipole_Moment_Debye"),
                "Total_SCF_Energy_Eh": xtb_props.get("Total_SCF_Energy_Eh"),
                "Pocket_Backbone_RMSD": pocket_bb_rmsd,
                "Ligand_Self_Aligned_RMSD": ligand_self_aligned_rmsd,
                "Min_WBO": min_wbo,
                "SMILES": smi_val,
                "test_dataset": dataset
            }

    # Build QM lookup maps by name and canonical SMILES
    def canon_smi(s):
        if not s or s == "Unknown":
            return None
        try:
            from rdkit import Chem
            m = Chem.MolFromSmiles(s)
            return Chem.MolToSmiles(m) if m else None
        except Exception:
            return None

    qm_by_name = {}
    qm_by_canon = {}
    matched_qm_names = set()

    for q_name, q_rec in qm_records.items():
        qm_by_name[q_name] = q_rec
        qm_by_name[q_name.lower()] = q_rec
        jn = q_name[4:] if q_name.startswith("run_") else q_name
        if jn.endswith("_model"):
            jn = jn[:-6]
        qm_by_name[jn] = q_rec
        qm_by_name[jn.lower()] = q_rec
        
        cs = canon_smi(q_rec.get("SMILES"))
        if cs:
            qm_by_canon[cs] = q_rec
            qm_by_canon[(q_rec["Method"].lower(), q_rec["test_dataset"].lower(), cs)] = q_rec

    # -------------------------------------------------------------------------
    # 2. Parse AF3 structural predictions from af3_path
    # -------------------------------------------------------------------------
    af3_records = []
    if af3_path.exists():
        af3_dirs = [d for d in af3_path.iterdir() if d.is_dir() and not d.name.startswith(("_", "."))]
        print(f"Extracting AF3 docking confidence metrics across {len(af3_dirs)} candidate predictions...")
        for td in af3_dirs:
            raw_name = td.name
            job_name = raw_name[4:] if raw_name.startswith("run_") else raw_name
            method = "Unknown"
            dataset = "Unknown"
            model_idx = None
            
            for d in ["glaucoma_bulk", "glaucoma", "cardiomyopathy", "parkinsons", "hcc", "lincs1000_1", "lincs1000_2"]:
                if d in job_name:
                    dataset = d
                    break

            parts = job_name.split("_")
            if len(parts) >= 2:
                try:
                    model_idx = int(parts[-1])
                except ValueError:
                    model_idx = None
                if parts[0] == "docking":
                    method = parts[1]
                else:
                    method = parts[0]
            
            # 1. AF3 Confidence JSON
            af3_json = td / f"{job_name}_summary_confidences.json"
            if not af3_json.exists():
                af3_json = td / f"{raw_name}_summary_confidences.json"
            iptm = None
            ptm = None
            ranking_score = None
            if af3_json.exists():
                try:
                    with open(af3_json, "r") as f:
                        conf = json.load(f)
                        iptm = conf.get("iptm")
                        ptm = conf.get("ptm")
                        ranking_score = conf.get("ranking_score")
                except Exception:
                    pass

            # 2. CIF pLDDT scores
            cif_file = td / f"{job_name}_model.cif"
            if not cif_file.exists():
                cif_file = td / f"{raw_name}_model.cif"
            mean_prot_plddt = None
            mean_lig_plddt = None
            if cif_file.exists():
                prot_plddts = []
                lig_plddts = []
                with open(cif_file, "r") as cf:
                    for line in cf:
                        if line.startswith(("ATOM", "HETATM")):
                            p = line.split()
                            if len(p) >= 15:
                                try:
                                    res = p[5].upper()
                                    biso = float(p[14])
                                    if "LIG" in res:
                                        lig_plddts.append(biso)
                                    else:
                                        prot_plddts.append(biso)
                                except Exception:
                                    pass
                if prot_plddts:
                    mean_prot_plddt = round(sum(prot_plddts) / len(prot_plddts), 2)
                if lig_plddts:
                    mean_lig_plddt = round(sum(lig_plddts) / len(lig_plddts), 2)
            
            # 3. Extract SMILES from _data.json
            smiles_val = None
            for cand in [td / f"{raw_name}_data.json", td / f"{job_name}_data.json"]:
                if cand.exists():
                    try:
                        with open(cand, "r") as dj:
                            dj_data = json.load(dj)
                            for s in dj_data.get("sequences", []):
                                if "ligand" in s:
                                    smiles_val = s["ligand"].get("smiles")
                                    break
                    except Exception:
                        pass
                if smiles_val:
                    break

            if not smiles_val:
                smiles_val = job_meta_map.get(raw_name, job_meta_map.get(job_name, "Unknown"))
            
            # Match with QM data (by name or canonical SMILES)
            qm_match = qm_by_name.get(raw_name, qm_by_name.get(job_name, qm_by_name.get(f"run_{job_name}_model", qm_by_name.get(f"run_{raw_name}_model", {}))))
            if not qm_match and smiles_val and smiles_val != "Unknown":
                cs = canon_smi(smiles_val)
                if cs:
                    qm_match = qm_by_canon.get((method.lower(), dataset.lower(), cs), qm_by_canon.get(cs, {}))
            
            if qm_match:
                matched_qm_names.add(qm_match.get("Directory", ""))

            af3_records.append({
                "Directory": raw_name,
                "Method": method,
                "Model_Index": model_idx,
                "Ranking_Score": ranking_score,
                "Ligand_RMSD_Pocket_Aligned": qm_match.get("Ligand_RMSD_Pocket_Aligned"),
                "HOMO_LUMO_Gap_eV": qm_match.get("HOMO_LUMO_Gap_eV"),
                "Total_Dipole_Moment_Debye": qm_match.get("Total_Dipole_Moment_Debye"),
                "Total_SCF_Energy_Eh": qm_match.get("Total_SCF_Energy_Eh"),
                "iPTM": iptm,
                "PTM": ptm,
                "Pocket_Backbone_RMSD": qm_match.get("Pocket_Backbone_RMSD"),
                "Ligand_Self_Aligned_RMSD": qm_match.get("Ligand_Self_Aligned_RMSD"),
                "Min_WBO": qm_match.get("Min_WBO"),
                "mean_protein_plddt": mean_prot_plddt,
                "mean_ligand_plddt": mean_lig_plddt,
                "SMILES": smiles_val,
                "test_dataset": dataset
            })

    # Also add any QM records that weren't directly merged into af3_records
    for q_name, q_rec in qm_records.items():
        if q_name not in matched_qm_names:
            af3_records.append({
                "Directory": q_rec["Directory"],
                "Method": q_rec["Method"],
                "Model_Index": q_rec["Model_Index"],
                "Ranking_Score": None,
                "Ligand_RMSD_Pocket_Aligned": q_rec["Ligand_RMSD_Pocket_Aligned"],
                "HOMO_LUMO_Gap_eV": q_rec["HOMO_LUMO_Gap_eV"],
                "Total_Dipole_Moment_Debye": q_rec["Total_Dipole_Moment_Debye"],
                "Total_SCF_Energy_Eh": q_rec["Total_SCF_Energy_Eh"],
                "iPTM": None,
                "PTM": None,
                "Pocket_Backbone_RMSD": q_rec["Pocket_Backbone_RMSD"],
                "Ligand_Self_Aligned_RMSD": q_rec["Ligand_Self_Aligned_RMSD"],
                "Min_WBO": q_rec["Min_WBO"],
                "mean_protein_plddt": None,
                "mean_ligand_plddt": None,
                "SMILES": q_rec.get("SMILES"),
                "test_dataset": q_rec["test_dataset"]
            })

    df = pd.DataFrame(af3_records)
    if "Ranking_Score" in df.columns:
        df = df.sort_values(by="Ranking_Score", ascending=False)

    out_file = Path(output_csv).resolve()
    os.makedirs(out_file.parent, exist_ok=True)
    df.to_csv(out_file, index=False)
    print(f"Successfully saved {len(df)} docking & QM metrics to: {out_file}")
    return df


def main():
    parser = argparse.ArgumentParser(
        description="AlphaFold 3 docking & semi-empirical xTB quantum mechanical calculation pipeline."
    )
    parser.add_argument("--run_af3", action="store_true", help="Launch AlphaFold 3 parallel batch inference.")
    parser.add_argument("--run_xtb", action="store_true", help="Execute active site extraction and GFN2-xTB calculations.")
    parser.add_argument("--extract_metrics", action="store_true", help="Extract docking confidence and xTB descriptors.")
    parser.add_argument(
        "--af3_json", 
        type=str, 
        default="TheraFlow/outputs/af3_batch_input.json",
        help="Input JSON containing AF3 docking jobs."
    )
    parser.add_argument(
        "--af3_dir", 
        type=str, 
        default="/users/stgrad/rojashrj/Projects/omics2mol/alphafold3",
        help="Path to AlphaFold 3 repository installation."
    )
    parser.add_argument(
        "--af3_output", 
        type=str, 
        default="TheraFlow/outputs/af3_predictions",
        help="Output directory for AF3 predictions."
    )
    parser.add_argument(
        "--qm_output", 
        type=str, 
        default="TheraFlow/outputs/docking_md_runs",
        help="Base output directory for xTB QM calculations."
    )
    parser.add_argument(
        "--output_metrics", 
        type=str, 
        default="TheraFlow/outputs/docking_qm_metrics.csv",
        help="Output CSV for extracted docking & QM metrics."
    )
    parser.add_argument("--n_jobs_xtb", type=int, default=4, help="Parallel workers for xTB QM calculations.")
    parser.add_argument("--cutoff", type=float, default=4.0, help="Interface cutoff in Angstroms.")
    
    args = parser.parse_args()
    
    if args.run_af3:
        run_af3_batch_pipeline(
            input_json_path=args.af3_json,
            af3_dir=args.af3_dir,
            output_dir=args.af3_output
        )
        
    if args.run_xtb:
        run_xtb_batch_pipeline(
            af3_dir=args.af3_output,
            qm_dir=args.qm_output,
            n_jobs=args.n_jobs_xtb,
            cutoff=args.cutoff
        )
        
    if args.extract_metrics or args.run_xtb:
        extract_docking_and_qm_metrics(
            af3_dir=args.af3_output,
            qm_dir=args.qm_output,
            input_json_path=args.af3_json,
            output_csv=args.output_metrics
        )
        
    if not args.run_af3 and not args.run_xtb and not args.extract_metrics:
        # Default: run xTB batch and extract metrics
        run_xtb_batch_pipeline(
            af3_dir=args.af3_output,
            qm_dir=args.qm_output,
            n_jobs=args.n_jobs_xtb,
            cutoff=args.cutoff
        )
        extract_docking_and_qm_metrics(
            af3_dir=args.af3_output,
            qm_dir=args.qm_output,
            input_json_path=args.af3_json,
            output_csv=args.output_metrics
        )


if __name__ == "__main__":
    main()
