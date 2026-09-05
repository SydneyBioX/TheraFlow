#!/usr/bin/env python3
"""
Step 7: AlphaFold 3 Docking & QM/MM Molecular Dynamics Setup, Execution & Full Metric Extraction
Standalone pipeline for:
  1. AlphaFold 3 parallel batch inference (accelerated MSA caching + multi-GPU splitting).
  2. Active site pocket interface extraction (integrated take_interface logic).
  3. Semi-empirical xTB geometry optimization & ORCA DFT single-point setup & execution.
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
import shutil
import argparse
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd


# =============================================================================
# 1. ALPHAFOLD 3 BATCH INFERENCE RUNNER
# =============================================================================

def clean_job_for_af3_dialect(job):
    """Ensures AF3 JSON dialect compliance and removes deprecated fields."""
    job_copy = json.loads(json.dumps(job))
    job_copy["dialect"] = "alphafold3"
    job_copy["version"] = 1
    for seq in job_copy.get("sequences", []):
        if "protein" in seq and "count" in seq["protein"]:
            del seq["protein"]["count"]
        elif "ligand" in seq and "count" in seq["ligand"]:
            del seq["ligand"]["count"]
    return job_copy


def run_af3_batch_pipeline(
    input_json_path="TheraFlow/outputs/af3_batch_input.json",
    af3_dir="alphafold3",
    output_dir="TheraFlow/outputs/af3_predictions",
    af3_python="/verona/nobackup/rojashrj/miniconda3/envs/af3_312/bin/python",
    gpus=("0", "1", "2"),
    model_dir="weights",
    db_dir="database"
):
    """
    Executes high-throughput AlphaFold 3 docking:
    1. Extracts unique protein sequences across all candidate jobs.
    2. Runs database search (MSA/templates) ONCE per unique protein.
    3. Caches MSAs and injects them into all candidate ligand docking jobs.
    4. Distributes inference across GPUs in parallel.
    """
    input_json = Path(input_json_path).resolve()
    if not input_json.exists():
        raise FileNotFoundError(f"AF3 input JSON not found at: {input_json}")
        
    af3_path = Path(af3_dir).resolve()
    if not af3_path.exists():
        raise FileNotFoundError(f"AlphaFold 3 directory not found at: {af3_path}")
        
    out_path = Path(output_dir).resolve()
    os.makedirs(out_path, exist_ok=True)
    
    if not os.path.exists(af3_python):
        af3_python = sys.executable
        
    with open(input_json, "r") as f:
        jobs = json.load(f)
        
    print(f"Loaded {len(jobs)} AF3 docking jobs from {input_json}")
    
    # Identify unique protein sequences
    unique_seqs = {}
    for job in jobs:
        for seq in job.get("sequences", []):
            if "protein" in seq:
                prot_seq = seq["protein"]["sequence"]
                if prot_seq not in unique_seqs:
                    job_cleaned = clean_job_for_af3_dialect(job)
                    job_cleaned["modelSeeds"] = [1]
                    unique_seqs[prot_seq] = job_cleaned
                break
                
    representative_jobs = list(unique_seqs.values())
    print(f"Identified {len(representative_jobs)} unique protein targets for MSA precomputation.")
    
    work_dir = out_path / "_af3_workdir"
    precompute_in = work_dir / "precompute_input_batch"
    precompute_out = work_dir / "output_precompute"
    os.makedirs(work_dir, exist_ok=True)
    
    # Phase 1 & 2: Check for existing cached MSAs first across output_precompute directories
    msa_cache = {}
    search_dirs = [
        af3_path / "output_precompute",
        work_dir / "output_precompute",
        Path("alphafold3/output_precompute").resolve()
    ]
    for s_dir in search_dirs:
        if s_dir.exists():
            for dj in s_dir.glob("*/*_data.json"):
                try:
                    with open(dj, "r") as f:
                        data = json.load(f)
                    for seq in data.get("sequences", []):
                        if "protein" in seq and seq["protein"].get("sequence"):
                            prot_seq = seq["protein"]["sequence"]
                            if prot_seq not in msa_cache and seq["protein"].get("unpairedMsa"):
                                msa_cache[prot_seq] = {
                                    "unpairedMsa": seq["protein"].get("unpairedMsa"),
                                    "pairedMsa": seq["protein"].get("pairedMsa"),
                                    "templates": seq["protein"].get("templates", [])
                                }
                except Exception:
                    pass
                    
    missing_jobs = [j for j in representative_jobs if not any(s.get("protein", {}).get("sequence") in msa_cache for s in j.get("sequences", []))]
    
    af3_env_bin = str(Path(af3_python).parent)
    env = os.environ.copy()
    env["PATH"] = f"{af3_env_bin}:{env.get('PATH', '')}"

    if missing_jobs:
        print(f"\n--- Phase 1: Running MSA precomputation on {len(missing_jobs)} missing targets ---")
        for d in [precompute_in, precompute_out]:
            if d.exists():
                shutil.rmtree(d)
            os.makedirs(d, exist_ok=True)
            
        for job in missing_jobs:
            with open(precompute_in / f"{job['name']}.json", "w") as f:
                json.dump(job, f, indent=2)
                
        precompute_cmd = [
            af3_python,
            str(af3_path / "run_alphafold.py"),
            f"--input_dir={precompute_in}",
            "--run_inference=False",
            f"--model_dir={af3_path / model_dir}",
            f"--db_dir={af3_path / db_dir}",
            f"--output_dir={precompute_out}",
            "--force_output_dir=True"
        ]
        
        t0 = time.time()
        subprocess.run(precompute_cmd, check=True, env=env, cwd=str(af3_path))
        print(f"MSA Precomputation completed in {time.time() - t0:.2f} seconds.")
        
        for job in missing_jobs:
            job_name = job["name"]
            data_json = precompute_out / job_name / f"{job_name}_data.json"
            if data_json.exists():
                with open(data_json, "r") as f:
                    data = json.load(f)
                for seq in data.get("sequences", []):
                    if "protein" in seq:
                        prot_seq = seq["protein"]["sequence"]
                        msa_cache[prot_seq] = {
                            "unpairedMsa": seq["protein"].get("unpairedMsa"),
                            "pairedMsa": seq["protein"].get("pairedMsa"),
                            "templates": seq["protein"].get("templates", [])
                        }
                        break
    else:
        print(f"Found precomputed MSAs in cache for all {len(representative_jobs)} target sequences! Skipping database search.")
        
    print(f"Successfully cached MSAs for {len(msa_cache)} unique sequences.")

    
    # Phase 3: Populate all candidate jobs with cached MSAs
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
        
    # Phase 4: Split into GPU chunks and launch parallel inference
    chunk_size = (len(populated_jobs) + len(gpus) - 1) // len(gpus)
    subprocs = []
    
    print(f"\n--- Phase 2: Launching AF3 parallel inference across {len(gpus)} GPUs ---")
    for idx, gpu in enumerate(gpus):
        chunk = populated_jobs[idx * chunk_size : (idx + 1) * chunk_size]
        if not chunk:
            continue
        gpu_dir = work_dir / f"populated_input_batch_gpu{gpu}"
        os.makedirs(gpu_dir, exist_ok=True)
        for j in chunk:
            with open(gpu_dir / f"{j['name']}.json", "w") as f:
                json.dump(j, f, indent=2)
                
        run_cmd = [
            af3_python,
            str(af3_path / "run_alphafold.py"),
            f"--input_dir={gpu_dir}",
            "--run_data_pipeline=False",
            f"--model_dir={af3_path / model_dir}",
            f"--db_dir={af3_path / db_dir}",
            f"--output_dir={out_path}",
            "--force_output_dir=True"
        ]
        
        log_file = open(out_path / f"run_af3_gpu{gpu}.log", "w")
        gpu_env = env.copy()
        gpu_env["CUDA_VISIBLE_DEVICES"] = str(gpu)
        
        print(f"  GPU {gpu}: Running {len(chunk)} jobs (log: {out_path / f'run_af3_gpu{gpu}.log'})")
        p = subprocess.Popen(run_cmd, stdout=log_file, stderr=subprocess.STDOUT, env=gpu_env, cwd=str(af3_path))
        subprocs.append((p, log_file, gpu))
        
    for p, log, gpu in subprocs:
        p.wait()
        log.close()
        print(f"  GPU {gpu} finished with return code {p.returncode}.")
        
    print(f"\nAll AF3 structure predictions saved to: {out_path}")
    return str(out_path)


# =============================================================================
# 2. CIF TO PDB CONVERSION & INTERFACE EXTRACTION
# =============================================================================

def convert_cif_to_pdb(cif_path, pdb_path):
    """Parses mmCIF file and writes clean standard PDB format."""
    cif_path = Path(cif_path)
    pdb_path = Path(pdb_path)
    
    with open(cif_path, "r") as f:
        lines = f.readlines()
        
    atom_site_headers = []
    atom_records = []
    in_loop = False
    in_atom_site = False
    
    for line in lines:
        stripped = line.strip()
        if stripped == "loop_":
            in_loop = True
            atom_site_headers = []
            continue
        if in_loop:
            if stripped.startswith("_atom_site."):
                in_atom_site = True
                atom_site_headers.append(stripped)
                continue
            elif in_atom_site:
                if stripped.startswith("_") or stripped.startswith("#"):
                    in_loop = False
                    in_atom_site = False
                elif stripped:
                    atom_records.append(stripped)
                else:
                    in_loop = False
                    in_atom_site = False
                    
    if not atom_site_headers:
        raise ValueError(f"Could not find _atom_site loop in {cif_path}")
        
    header_map = {h.replace("_atom_site.", ""): idx for idx, h in enumerate(atom_site_headers)}
    
    os.makedirs(pdb_path.parent, exist_ok=True)
    with open(pdb_path, "w") as out:
        for rec in atom_records:
            parts = rec.split()
            if not parts:
                continue
            group = parts[header_map.get("group_PDB", 0)]
            if group not in ("ATOM", "HETATM"):
                continue
            serial = int(parts[header_map.get("id", 1)])
            element = parts[header_map.get("type_symbol", 2)]
            atom_name = parts[header_map.get("label_atom_id", 3)]
            res_name = parts[header_map.get("label_comp_id", 5)]
            chain = parts[header_map.get("label_asym_id", 6)]
            res_seq = parts[header_map.get("label_seq_id", 8)]
            
            try:
                res_seq_num = int(res_seq)
            except ValueError:
                res_seq_num = 1
                
            x = float(parts[header_map.get("Cartn_x", 10)])
            y = float(parts[header_map.get("Cartn_y", 11)])
            z = float(parts[header_map.get("Cartn_z", 12)])
            occ = float(parts[header_map.get("occupancy", 13)]) if "occupancy" in header_map else 1.00
            biso = float(parts[header_map.get("B_iso_or_equiv", 14)]) if "B_iso_or_equiv" in header_map else 0.00
            
            if "LIG" in res_name.upper():
                res_name = "LIG"
            res_name = res_name[:3]
            
            atom_name_fmt = f" {atom_name:<3}" if len(atom_name) < 4 and not atom_name.startswith("H") else f"{atom_name:<4}"
            
            out.write(
                f"{group:<6}{serial:>5} {atom_name_fmt} {res_name:>3} {chain:1}{res_seq_num:>4}    "
                f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{occ:>6.2f}{biso:>6.2f}          {element:>2}\n"
            )
        out.write("END\n")
    return pdb_path


# =============================================================================
# 3. ORCA & xTB QM/MM DIRECTORY SETUP & EXECUTION
# =============================================================================

def setup_single_docking_md(
    cif_file,
    outdir,
    orca_path="/users/stgrad/rojashrj/software/orca/orca_6_1_0_linux_x86-64_shared_openmpi418_avx2",
    ligand="LIG",
    cutoff=4.0,
    charge=0,
    nproc=16
):
    """Sets up a complete ORCA and xTB calculation workspace for a given AF3 mmCIF prediction."""
    cif_path = Path(cif_file).resolve()
    cif_name = cif_path.stem
    outdir_path = Path(outdir).resolve()
    
    for d in ["01_prepare", "03_cluster", "04_xtb", "05_orca"]:
        os.makedirs(outdir_path / d, exist_ok=True)
        
    # 1. Convert mmCIF to PDB
    pdb_path = outdir_path / "01_prepare" / "complex_clean.pdb"
    convert_cif_to_pdb(cif_path, pdb_path)
    
    # 2. Copy take_interface.py helper if available
    pipeline_take_iface = Path(__file__).parent / "take_interface.py"
    if pipeline_take_iface.exists():
        shutil.copy2(pipeline_take_iface, outdir_path / "03_cluster" / "take_interface.py")
        
    # 3. Generate master run_pipeline.sh
    pipeline_sh = outdir_path / "run_pipeline.sh"
    orca_abs = os.path.abspath(orca_path)
    
    script = f"""#!/bin/bash
set -e
export PATH="{orca_abs}:${{PATH}}"
export LD_LIBRARY_PATH="{orca_abs}:${{LD_LIBRARY_PATH}}"
export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2

echo "=========================================================="
echo "Starting ORCA & xTB Analysis Pipeline for {cif_name}"
echo "=========================================================="

echo "--- Step 1: Running pdb2pqr ---"
cd 01_prepare
if command -v pdb2pqr &> /dev/null; then
    pdb2pqr --ff=PARSE --titration-state-method=propka --with-ph=7.4 complex_clean.pdb complex_parse.pqr
    grep "^HETATM" complex_clean.pdb >> complex_parse.pqr || true
else
    echo "[Warning] pdb2pqr not found; using complex_clean.pdb as fallback"
    cp complex_clean.pdb complex_parse.pqr
fi
cd ..

echo "--- Step 2: Extracting active site cluster ({cutoff} A) ---"
cd 03_cluster
if [ -f "take_interface.py" ]; then
    python take_interface.py ../01_prepare/complex_parse.pqr {ligand.upper()} {cutoff:.1f} cluster {charge}
fi
cd ..

echo "--- Step 3: Running xTB Geometry Optimization ---"
cd 04_xtb
Q=$(cat ../03_cluster/cluster_charge.txt 2>/dev/null || echo "{charge}")
if command -v xtb &> /dev/null; then
    xtb ../03_cluster/cluster.xyz --gfn 2 --opt loose --alpb water --chrg ${{Q}} --uhf 0 > cluster_xtb.out
    cp xtbopt.xyz ../05_orca/cluster_opt.xyz
else
    echo "[Warning] xtb not found in PATH; copying raw cluster coordinates"
    cp ../03_cluster/cluster.xyz ../05_orca/cluster_opt.xyz
fi
cd ..

echo "--- Step 4: Generating ORCA Input ---"
Q=$(cat 03_cluster/cluster_charge.txt 2>/dev/null || echo "{charge}")
MULT=$(python3 -c "
import sys
z_map = {{H:1,HE:2,LI:3,BE:4,B:5,C:6,N:7,O:8,F:9,NE:10,NA:11,MG:12,AL:13,SI:14,P:15,S:16,CL:17,AR:18,K:19,CA:20,FE:26,ZN:30,BR:35,I:53}}
atoms = []
try:
    with open('05_orca/cluster_opt.xyz') as f:
        for line in list(f)[2:]:
            p = line.split()
            if p: atoms.append(p[0].upper())
    total_z = sum(z_map.get(a, 6) for a in atoms)
    q = int(sys.argv[1])
    print(1 if (total_z - q) % 2 == 0 else 2)
except Exception:
    print(1)
" ${{Q}})

cat << EOF_ORCA > 05_orca/cluster_sp.inp
! wB97X-V def2-SVP def2/J RIJCOSX TightSCF CPCM(Water) HIRSHFELD CHELPG

%pal
  nprocs {nproc}
end

%chelpg
  GRID 0.3
  RMAX 2.8
  DIPOLE FALSE
end

%output
  Print[P_Hirshfeld] 1
  Print[P_BondOrder_M] 1
end

* xyzfile ${{Q}} ${{MULT}} cluster_opt.xyz
EOF_ORCA

echo "--- Step 5: Ready for ORCA calculation ---"
echo "To run ORCA: cd 05_orca && orca cluster_sp.inp > cluster_sp.out"
echo "=== Setup complete in {outdir_path} ==="
"""
    with open(pipeline_sh, "w") as f:
        f.write(script)
    os.chmod(pipeline_sh, 0o755)
    return str(outdir_path)


# =============================================================================
# 4. KABSCH RMSD & QM (xTB / ORCA) METRIC PARSING
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
    energy_matches = re.findall(r"TOTAL ENERGY\s+([\d.-]+)\s+Eh", content)
    gap_matches = re.findall(r"HOMO-LUMO GAP\s+([\d.-]+)\s+eV", content)
    dipole_pos = content.rfind("molecular dipole:")
    dipole_moment = None
    if dipole_pos != -1:
        dipole_block = content[dipole_pos:dipole_pos+500]
        m = re.search(r"full:\s+[\d.-]+\s+[\d.-]+\s+[\d.-]+\s+([\d.-]+)", dipole_block)
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
    af3_dir="TheraFlow/outputs/af3_predictions",
    qm_dir="TheraFlow/outputs/docking_md_runs",
    input_json_path="TheraFlow/outputs/af3_batch_input.json",
    output_csv="TheraFlow/outputs/docking_qm_metrics.csv"
):
    """
    Extracts unified AlphaFold 3 Docking & QM/MM metrics matching:
    Directory, Method, Model_Index, Ranking_Score, Ligand_RMSD_Pocket_Aligned,
    HOMO_LUMO_Gap_eV, Total_Dipole_Moment_Debye, Total_SCF_Energy_Eh,
    iPTM, PTM, Pocket_Backbone_RMSD, Ligand_Self_Aligned_RMSD, Min_WBO,
    mean_protein_plddt, mean_ligand_plddt, SMILES, test_dataset.
    """
    af3_path = Path(af3_dir).resolve()
    qm_path = Path(qm_dir).resolve()
    
    # Map input JSON metadata by job name
    job_meta_map = {}
    input_json = Path(input_json_path)
    if input_json.exists():
        try:
            with open(input_json, "r") as f:
                jobs = json.load(f)
            for j in jobs:
                name = j.get("name")
                lig_smi = ""
                for s in j.get("sequences", []):
                    if "ligand" in s:
                        lig_smi = s["ligand"].get("smiles", "")
                job_meta_map[name] = lig_smi
        except Exception as e:
            print(f"Note: Could not parse metadata from input JSON: {e}")
            
    # Find all job directories (from qm_dir or af3_dir)
    target_dirs = []
    if qm_path.exists():
        target_dirs.extend([d for d in qm_path.iterdir() if d.is_dir() and not d.name.startswith(("_", "."))])
    if not target_dirs and af3_path.exists():
        target_dirs.extend([d for d in af3_path.iterdir() if d.is_dir() and not d.name.startswith(("_", "."))])
        
    print(f"Extracting docking & QM metrics across {len(target_dirs)} candidate directories...")
    records = []
    
    for td in target_dirs:
        raw_name = td.name
        # Clean job name (e.g. run_docking_BiAAE_cardiomyopathy_1 -> docking_BiAAE_cardiomyopathy_1)
        job_name = raw_name[4:] if raw_name.startswith("run_") else raw_name
        
        # Parse method, dataset, index
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
        af3_json = af3_path / job_name / f"{job_name}_summary_confidences.json"
        if not af3_json.exists() and (af3_path / raw_name / f"{raw_name}_summary_confidences.json").exists():
            af3_json = af3_path / raw_name / f"{raw_name}_summary_confidences.json"
            
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
        cif_file = af3_path / job_name / f"{job_name}_model.cif"
        if not cif_file.exists() and (af3_path / raw_name / f"{raw_name}_model.cif").exists():
            cif_file = af3_path / raw_name / f"{raw_name}_model.cif"
            
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
                
        # 3. xTB and Kabsch RMSD Calculations
        qm_job_dir = qm_path / f"run_{job_name}" if (qm_path / f"run_{job_name}").exists() else qm_path / job_name
        if not qm_job_dir.exists():
            qm_job_dir = td
            
        pdb_path = qm_job_dir / "03_cluster" / "cluster.pdb"
        xyz_path = qm_job_dir / "04_xtb" / "xtbopt.xyz"
        xtb_out = qm_job_dir / "04_xtb" / "cluster_xtb.out"
        wbo_path = qm_job_dir / "04_xtb" / "wbo"
        
        xtb_props = parse_xtb_output(str(xtb_out))
        
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
                        pocket_bb_rmsd = calculate_rmsd(P_opt_aligned, P_start)
                        
                        if lig_idx:
                            L_start = start_coords[lig_idx]
                            L_opt = opt_coords[lig_idx]
                            L_opt_aligned = apply_transform(L_opt, c_P, c_Q, R)
                            ligand_pocket_aligned_rmsd = calculate_rmsd(L_opt_aligned, L_start)
                            
                    if lig_idx:
                        L_start = start_coords[lig_idx]
                        L_opt = opt_coords[lig_idx]
                        c_L_opt, c_L_start, R_L = kabsch_alignment(L_opt, L_start)
                        L_opt_self_aligned = apply_transform(L_opt, c_L_opt, c_L_start, R_L)
                        ligand_self_aligned_rmsd = calculate_rmsd(L_opt_self_aligned, L_start)
                        
                    min_wbo = parse_wbo_min(str(wbo_path), lig_idx, pock_idx)
            except Exception as e:
                print(f"Note: Could not calculate RMSD for {job_name}: {e}")
                
        smiles_val = job_meta_map.get(job_name, "Unknown")
        
        records.append({
            "Directory": raw_name,
            "Method": method,
            "Model_Index": model_idx,
            "Ranking_Score": ranking_score,
            "Ligand_RMSD_Pocket_Aligned": ligand_pocket_aligned_rmsd,
            "HOMO_LUMO_Gap_eV": xtb_props.get("HOMO_LUMO_Gap_eV"),
            "Total_Dipole_Moment_Debye": xtb_props.get("Total_Dipole_Moment_Debye"),
            "Total_SCF_Energy_Eh": xtb_props.get("Total_SCF_Energy_Eh"),
            "iPTM": iptm,
            "PTM": ptm,
            "Pocket_Backbone_RMSD": pocket_bb_rmsd,
            "Ligand_Self_Aligned_RMSD": ligand_self_aligned_rmsd,
            "Min_WBO": min_wbo,
            "mean_protein_plddt": mean_prot_plddt,
            "mean_ligand_plddt": mean_lig_plddt,
            "SMILES": smiles_val,
            "test_dataset": dataset
        })
        
    default_cols = ["name", "method", "test_dataset", "SMILES", "Ranking_Score", "iPTM", "PTM", "mean_protein_plddt", "mean_ligand_plddt", "Ligand_RMSD_Pocket_Aligned", "Pocket_Backbone_RMSD", "Ligand_Self_Aligned_RMSD", "HOMO_LUMO_Gap_eV", "Total_Dipole_Moment_Debye", "Total_SCF_Energy_Eh", "Min_WBO"]
    df = pd.DataFrame(records) if len(records) > 0 else pd.DataFrame(columns=default_cols)
    if len(df) > 0 and "Ranking_Score" in df.columns:
        df = df.sort_values(by="Ranking_Score", ascending=False).reset_index(drop=True)
        
    out_file = Path(output_csv).resolve()
    os.makedirs(out_file.parent, exist_ok=True)
    df.to_csv(out_file, index=False)
    print(f"Successfully saved {len(df)} docking & QM metrics to: {out_file}")
    return df


# =============================================================================
# 5. CLI INTERFACE
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="AlphaFold 3 parallel batch inference, docking metrics extraction, and QM/MM (xTB / ORCA) setup."
    )
    # AF3 Inference Options
    parser.add_argument("--run_af3", action="store_true", help="Launch AlphaFold 3 parallel batch inference.")
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
        "--af3_python", 
        type=str, 
        default="/verona/nobackup/rojashrj/miniconda3/envs/af3_312/bin/python",
        help="Python executable with AF3 and Jax installed."
    )
    parser.add_argument(
        "--af3_output", 
        type=str, 
        default="TheraFlow/outputs/af3_predictions",
        help="Output directory for AF3 predictions."
    )
    parser.add_argument(
        "--gpus", 
        nargs="+", 
        default=["0", "1", "2"],
        help="GPUs to split AF3 inference across."
    )
    
    # Metric Extraction Options
    parser.add_argument("--extract_metrics", action="store_true", help="Extract Ranking_Score, Ligand_RMSD_Pocket_Aligned, HOMO_LUMO_Gap, Dipole, and SCF Energy.")
    parser.add_argument(
        "--output_metrics", 
        type=str, 
        default="TheraFlow/outputs/docking_qm_metrics.csv",
        help="Output CSV for extracted docking & QM metrics."
    )
    
    # ORCA / xTB Setup Options
    parser.add_argument("--setup_qm", action="store_true", help="Prepare QM/MM cluster and ORCA/xTB scripts.")
    parser.add_argument("--cif", type=str, default=None, help="Single mmCIF prediction file.")
    parser.add_argument(
        "--cif_dir", 
        type=str, 
        default=None, 
        help="Directory of CIF files to batch setup QM/MM calculations for."
    )
    parser.add_argument(
        "--qm_output", 
        type=str, 
        default="TheraFlow/outputs/docking_md_runs",
        help="Base output directory for QM/MM setups."
    )
    parser.add_argument(
        "--orca_path", 
        type=str, 
        default="/users/stgrad/rojashrj/software/orca/orca_6_1_0_linux_x86-64_shared_openmpi418_avx2",
        help="Path to ORCA binary directory."
    )
    parser.add_argument("--cutoff", type=float, default=4.0, help="Interface cutoff in Angstroms.")
    parser.add_argument("--nproc", type=int, default=16, help="Processors for ORCA.")
    
    args = parser.parse_args()
    
    if args.run_af3:
        run_af3_batch_pipeline(
            input_json_path=args.af3_json,
            af3_dir=args.af3_dir,
            output_dir=args.af3_output,
            af3_python=args.af3_python,
            gpus=args.gpus
        )
        
    if args.extract_metrics and not args.run_af3:
        extract_docking_and_qm_metrics(
            af3_dir=args.af3_output,
            qm_dir=args.qm_output,
            input_json_path=args.af3_json,
            output_csv=args.output_metrics
        )
        
    if args.setup_qm or args.cif or args.cif_dir:
        cif_files = []
        if args.cif:
            cif_files.append(Path(args.cif))
        elif args.cif_dir:
            cif_files.extend(Path(args.cif_dir).rglob("*_model.cif"))
            if not cif_files:
                cif_files.extend(Path(args.cif_dir).rglob("*.cif"))
        else:
            def_af3 = Path(args.af3_output)
            if def_af3.exists():
                cif_files.extend(def_af3.rglob("*_model.cif"))
                if not cif_files:
                    cif_files.extend(def_af3.rglob("*.cif"))
                    
        if not cif_files:
            print("No mmCIF files found to setup QM/MM calculations.")
            return
            
        print(f"Setting up QM/MM calculations for {len(cif_files)} CIF predictions...")
        for c in cif_files:
            target_out = Path(args.qm_output) / f"run_{c.stem}"
            setup_single_docking_md(
                cif_file=str(c),
                outdir=str(target_out),
                orca_path=args.orca_path,
                cutoff=args.cutoff,
                nproc=args.nproc
            )
            print(f"  Prepared: {target_out} (Run: cd {target_out} && ./run_pipeline.sh)")
            
    if not args.run_af3 and not args.extract_metrics and not (args.setup_qm or args.cif or args.cif_dir):
        print("Please specify an action flag: --run_af3, --extract_metrics, or --setup_qm / --cif.")
        parser.print_help()


if __name__ == "__main__":
    main()
