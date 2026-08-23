#!/usr/bin/env python3
"""
Master Orchestration Script for Molecule Generation
===================================================
This script orchestrates the generation of molecules across:
  - 10 methods:
    1. SmilesGEN (VAE Generation)
    2. SmilesGEN (VAE Optimization)
    3. GExMolGen
    4. Triomphe
    5. TransGEM
    6. paccmann_rl
    7. GxRNN
    8. BiAAE
    9. BiCEV
    9. FAME
    10. GEx2sGen
  - 7 datasets:
    - Lincs1000_1 (MCF7 cell line, AKT-inhibitor-1-2)
    - Lincs1000_2 (A375 cell line, venlafaxine)
    - cardiomyopathy
    - glaucoma
    - glaucoma_bulk
    - hcc
    - parkinsons

Usage:
  python run_all_generation.py --method <method_name> --dataset <dataset_name>
  python run_all_generation.py --run-all
"""

import os
import sys
import subprocess
import argparse
import time

# Base path
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(BASE_DIR)
if not os.path.exists(os.path.join(PARENT_DIR, "SmilesGEN")) and os.path.exists(os.path.join(os.path.dirname(PARENT_DIR), "SmilesGEN")):
    MODEL_ROOT = os.path.dirname(PARENT_DIR)
else:
    MODEL_ROOT = PARENT_DIR

# Env path mapping
CONDA_ENV_DIR = "/verona/nobackup/rojashrj/miniconda3/envs"
ENV_INTERPRETERS = {
    "smilesgen": os.path.join(CONDA_ENV_DIR, "smilesgen/bin/python"),
    "triomphe": os.path.join(CONDA_ENV_DIR, "triomphe/bin/python"),
    "TransGEM": os.path.join(CONDA_ENV_DIR, "TransGEM/bin/python"),
    "paccmann_rl": os.path.join(CONDA_ENV_DIR, "paccmann_rl/bin/python"),
    "gxrnn": os.path.join(CONDA_ENV_DIR, "gxrnn/bin/python"),
    "biaee": os.path.join(CONDA_ENV_DIR, "biaee/bin/python"),
    "bicev": os.path.join(CONDA_ENV_DIR, "bicev/bin/python"),
}

# Cell line and perturbation mappings for the 7 datasets
DATASET_METADATA = {
    "Lincs1000_1": {
        "cell": "MCF7",
        "pert": "AKT-inhibitor-1-2",
        "keywords": ["AKT-inhibitor-1-2"],
        "gexmolgen_ctl": "Lincs1000_1_ctl.csv",
        "gexmolgen_pert": "Lincs1000_1_pert.csv",
        "smilesgen_file": "Lincs1000_1.csv",
        "smilesgen_mean_file": "Lincs1000_1_m.csv"
    },
    "Lincs1000_2": {
        "cell": "A375",
        "pert": "venlafaxine",
        "keywords": ["venlafaxine"],
        "gexmolgen_ctl": "Lincs1000_2_ctl.csv",
        "gexmolgen_pert": "Lincs1000_2_pert.csv",
        "smilesgen_file": "Lincs1000_2.csv",
        "smilesgen_mean_file": "Lincs1000_2_m.csv"
    },
    "cardiomyopathy": {
        "cell": "MCF7",
        "pert": "cardiomyopathy",
        "keywords": ["carvedilol", "metoprolol", "lisinopril", "enalapril", "valsartan", "losartan", "spironolactone"],
        "gexmolgen_ctl": "cardiomyopathy_ctl.csv",
        "gexmolgen_pert": "cardiomyopathy_disease.csv",
        "smilesgen_file": "cardiomyopathy.csv",
        "smilesgen_mean_file": "cardiomyopathy.csv"
    },
    "glaucoma": {
        "cell": "MCF7",
        "pert": "glaucoma",
        "keywords": ["timolol", "brimonidine", "latanoprost"],
        "gexmolgen_ctl": "glacoma_ctl.csv",
        "gexmolgen_pert": "glacoma_disease.csv",
        "smilesgen_file": "glaucoma.csv",
        "smilesgen_mean_file": "glaucoma.csv"
    },
    "glaucoma_bulk": {
        "cell": "MCF7",
        "pert": "glaucoma_bulk",
        "keywords": ["timolol", "brimonidine", "latanoprost"],
        "gexmolgen_ctl": "glaucoma_bulk_ctl.csv",
        "gexmolgen_pert": "glaucoma_bulk_disease.csv",
        "smilesgen_file": "glaucoma_bulk.csv",
        "smilesgen_mean_file": "glaucoma_bulk.csv"
    },
    "hcc": {
        "cell": "MCF7",
        "pert": "hcc",
        "keywords": ["sorafenib", "regorafenib", "lenvatinib"],
        "gexmolgen_ctl": "hcc_ctl.csv",
        "gexmolgen_pert": "hcc_disease.csv",
        "smilesgen_file": "hcc.csv",
        "smilesgen_mean_file": "hcc.csv"
    },
    "parkinsons": {
        "cell": "MCF7",
        "pert": "parkinsons",
        "keywords": ["levodopa", "pramipexole", "ropinirole"],
        "gexmolgen_ctl": "parkinsons_ctl.csv",
        "gexmolgen_pert": "parkinsons_disease.csv",
        "smilesgen_file": "parkinsons.csv",
        "smilesgen_mean_file": "parkinsons.csv"
    }
}

METHODS = [
    "smilesgen_vae_gen",
    "gexmolgen",
    "triomphe",
    "transgem",
    "paccmann_rl",
    "gxrnn",
    "biaee",
    "bicev",
    "fame",
    "gex2sgen"
]

def run_smilesgen_vae_gen(dataset):
    meta = DATASET_METADATA[dataset]
    python_bin = ENV_INTERPRETERS["smilesgen"]
    cwd = os.path.join(MODEL_ROOT, "SmilesGEN")
    
    protein_name = meta["pert"]
    if "Lincs" in dataset:
        protein_name = protein_name + "_mean"
        
    cmd = [python_bin, "main.py", "--generation", f"--cell_name={meta['cell']}", f"--protein_name={protein_name}"]
    print(f"Executing SmilesGEN VAE Gen: {' '.join(cmd)}")
    subprocess.run(cmd, cwd=cwd, check=True)

def run_gexmolgen(dataset):
    meta = DATASET_METADATA[dataset]
    python_bin = ENV_INTERPRETERS["smilesgen"]
    cwd = os.path.join(MODEL_ROOT, "GExMolGen")
    
    control_path = os.path.join(MODEL_ROOT, "data/test_datasets/gexmolgen", meta["gexmolgen_ctl"])
    pert_path = os.path.join(MODEL_ROOT, "data/test_datasets/gexmolgen", meta["gexmolgen_pert"])
    output_path = os.path.join(MODEL_ROOT, f"GExMolGen/results/{dataset}_gexmolgen.csv")
    
    cmd = [
        python_bin, os.path.join(BASE_DIR, "run_gexmolgen_cli.py"),
        "--control", control_path,
        "--pert", pert_path,
        "--output", output_path,
        "--beam_width", "100",
        "--seed", "0"
    ]
    print(f"Executing GExMolGen CLI: {' '.join(cmd)}")
    subprocess.run(cmd, cwd=cwd, check=True)

def run_dynamic_wrapper(method_dir, method_name, env_name, dataset, template_base_name):
    meta = DATASET_METADATA[dataset]
    python_bin = ENV_INTERPRETERS[env_name]
    cwd = os.path.join(MODEL_ROOT, method_dir)
    
    template_path = os.path.join(cwd, f"run_{template_base_name}_cardiomyopathy.py")
    temp_script = os.path.join(cwd, f"run_{template_base_name}_{dataset}_temp.py")
    
    if not os.path.exists(template_path):
        print(f"Error: Template script {template_path} not found.")
        return
        
    with open(template_path, 'r') as f:
        content = f.read()
        
    new_content = content.replace("cardiomyopathy", dataset)
    new_content = new_content.replace("Cardiomyopathy", dataset.capitalize())
    
    if method_name == "triomphe":
        old_kw_str = 'cardiomyopathy_keywords = ["carvedilol", "metoprolol", "lisinopril", "enalapril", "valsartan", "losartan", "spironolactone"]'
        new_kw_str = f"{dataset}_keywords = {str(meta['keywords'])}"
        new_content = new_content.replace(old_kw_str, new_kw_str)
        new_content = new_content.replace("sig_info['pert_iname'].str.lower().isin(cardiomyopathy_keywords)", f"sig_info['pert_iname'].str.lower().isin({dataset}_keywords)")
        
    with open(temp_script, 'w') as f:
        f.write(new_content)
        
    cmd = [python_bin, os.path.basename(temp_script)]
    print(f"Executing {method_name.upper()} on {dataset}: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, cwd=cwd, check=True)
    finally:
        if os.path.exists(temp_script):
            os.remove(temp_script)

def run_fame(dataset):
    run_dynamic_wrapper("reimplementation", "fame", "smilesgen", dataset, "fame")

def run_gex2sgen(dataset):
    run_dynamic_wrapper("reimplementation", "gex2sgen", "smilesgen", dataset, "gex2sgen")

def main():
    parser = argparse.ArgumentParser(description="Orchestrator to generate molecules across all methods and datasets")
    parser.add_argument("--method", choices=METHODS, help="Run a specific method")
    parser.add_argument("--dataset", choices=list(DATASET_METADATA.keys()), help="Run on a specific dataset")
    parser.add_argument("--run-all", action="store_true", help="Run all combinations sequentially")
    
    args = parser.parse_args()
    
    if not args.run_all and (not args.method or not args.dataset):
        parser.print_help()
        sys.exit(1)
        
    log_dir = os.path.join(BASE_DIR, "master_run_logs")
    os.makedirs(log_dir, exist_ok=True)
    
    combinations = []
    if args.run_all:
        for m in METHODS:
            for d in DATASET_METADATA.keys():
                combinations.append((m, d))
    else:
        combinations.append((args.method, args.dataset))
        
    print(f"Orchestrating {len(combinations)} generation tasks...")
    summary = []
    
    for idx, (method, dataset) in enumerate(combinations, 1):
        print(f"\n[{idx}/{len(combinations)}] Starting {method} on {dataset}...")
        start_time = time.time()
        status = "SUCCESS"
        
        try:
            if method == "smilesgen_vae_gen":
                run_smilesgen_vae_gen(dataset)
            elif method == "gexmolgen":
                run_gexmolgen(dataset)
            elif method == "triomphe":
                run_dynamic_wrapper("Triomphe/software/code/structure_generator", "triomphe", "triomphe", dataset, "triomphe")
            elif method == "transgem":
                run_dynamic_wrapper("TransGEM", "transgem", "TransGEM", dataset, "transgem")
            elif method == "paccmann_rl":
                run_dynamic_wrapper("paccmann_rl", "paccmann_rl", "paccmann_rl", dataset, "paccmann")
            elif method == "gxrnn":
                run_dynamic_wrapper("GxRNN/GxRNN", "gxrnn", "gxrnn", dataset, "gxrnn")
            elif method == "biaee":
                run_dynamic_wrapper("BiAAE/training_scripts", "biaee", "biaee", dataset, "biaee")
            elif method == "bicev":
                run_dynamic_wrapper("BiCEV", "bicev", "bicev", dataset, "bicev")
            elif method == "fame":
                run_fame(dataset)
            elif method == "gex2sgen":
                run_gex2sgen(dataset)
        except Exception as e:
            print(f"Error during execution: {e}")
            status = "FAILED"
            
        elapsed = time.time() - start_time
        print(f"Result: {status} (took {elapsed:.2f}s)")
        summary.append((method, dataset, status, elapsed))
        
    print("\n=======================================================")
    print("Execution Summary:")
    print("=======================================================")
    for method, dataset, status, elapsed in summary:
        print(f"  {method:<20} | {dataset:<15} : {status:<10} ({elapsed:.2f}s)")
    print("=======================================================")

if __name__ == "__main__":
    main()
