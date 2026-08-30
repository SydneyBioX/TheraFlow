# smilesgen repositroy should be cloned prior to running this script using gitclone  https://github.com/hliulab/SmilesGEN.git
# since paths to the data is hardcoded,  we need to place them in approprioate directores. (check setup_symlinks function )
# conda activate smilesgen (activate the environment which has all dependencies)
import os
import sys
import argparse
import subprocess
import shutil

# Paths setup
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
WORKSPACE_ROOT = os.path.dirname(os.path.dirname(BASE_DIR))
DATA_DIR = os.path.join(BASE_DIR, "data")
REST_METHODS_DIR = os.path.join(DATA_DIR, "test_datasets/rest_methods")
SMILESGEN_MODEL_DIR = os.path.join(WORKSPACE_ROOT, "SmilesGEN")
SMILESGEN_TEST_DIR = os.path.join(SMILESGEN_MODEL_DIR, "datasets/test_protein")
OUTPUTS_DIR = os.path.join(BASE_DIR, "outputs/smilesgen")
PYTHON_BIN = "/verona/nobackup/rojashrj/miniconda3/envs/smilesgen/bin/python"

DATASET_METADATA = {
    "Lincs1000_1": {
        "cell": "MCF7",
        "pert": "AKT-inhibitor-1-2",
        "src_file": "lincs1000_1.csv",
        "is_mean": False
    },
    "Lincs1000_2": {
        "cell": "A375",
        "pert": "venlafaxine",
        "src_file": "lincs1000_2.csv",
        "is_mean": False
    },
    "cardiomyopathy": {
        "cell": "MCF7",
        "pert": "cardiomyopathy",
        "src_file": "cardiomyopathy.csv",
        "is_mean": False
    },
    "parkinsons": {
        "cell": "MCF7",
        "pert": "parkinsons",
        "src_file": "parkinsons.csv",
        "is_mean": False
    },
    "hcc": {
        "cell": "MCF7",
        "pert": "hcc",
        "src_file": "hcc.csv",
        "is_mean": False
    },
    "glaucoma_bulk": {
        "cell": "MCF7",
        "pert": "glaucoma_bulk",
        "src_file": "glaucoma_bulk.csv",
        "is_mean": False
    },
    "glaucoma": {
        "cell": "MCF7",
        "pert": "glaucoma",
        "src_file": "glaucoma.csv",
        "is_mean": False
    }
}

def setup_symlinks(dataset):
    """
    Creates symlinks so that it matches the hardcoded directories 
    """
    meta = DATASET_METADATA[dataset]
    os.makedirs(SMILESGEN_TEST_DIR, exist_ok=True)
    
    src_path = os.path.join(REST_METHODS_DIR, meta["src_file"])
    if not os.path.exists(src_path):
        print(f"Error: Preprocessed file {src_path} does not exist. Run preprocessing first.")
        sys.exit(1)
        
    suffix = "_mean.csv" if meta["is_mean"] else ".csv"
    
    # Adjust pert name for SmilesGEN hardcoded loading
    pert_name = meta["pert"]
    if dataset == "Lincs1000_1":
        pert_name = "AKT1"
        
    tgt_name = f"{meta['cell']}_{pert_name}{suffix}"
    tgt_path = os.path.join(SMILESGEN_TEST_DIR, tgt_name)
    
    # Remove existing file/symlink if present
    if os.path.exists(tgt_path) or os.path.islink(tgt_path):
        os.remove(tgt_path)
    os.symlink(src_path, tgt_path)
    print(f"Created symlink: {tgt_path} -> {src_path}")
    
    # Special double mapping for AKT-inhibitor-1-2
    if meta["pert"] == "AKT-inhibitor-1-2":
        alt_name = f"{meta['cell']}_AKT-inhibitor-1-2{suffix}"
        alt_path = os.path.join(SMILESGEN_TEST_DIR, alt_name)
        if os.path.exists(alt_path) or os.path.islink(alt_path):
            os.remove(alt_path)
        os.symlink(src_path, alt_path)
        print(f"Created alternative symlink: {alt_path} -> {src_path}")

def run_pipeline(dataset):
    meta = DATASET_METADATA[dataset]
    
    protein_name = meta["pert"]
    if dataset == "Lincs1000_1":
        protein_name = "AKT1"
        
    # Command 1: VAE Generation
    cmd_gen = [
        PYTHON_BIN, 
        "main.py", 
        "--generation", 
        f"--cell_name={meta['cell']}", 
        f"--protein_name={protein_name}"
    ]
    
    print(f"\n[{dataset}] Step 1: Running SmilesGEN VAE Generation in {SMILESGEN_MODEL_DIR}:")
    print(" ".join(cmd_gen))
    subprocess.run(cmd_gen, cwd=SMILESGEN_MODEL_DIR, check=True)
    
    # Command 2: Molecular Optimization
    cmd_opt = [
        PYTHON_BIN, 
        "MolecularOptimization.py", 
        f"--cell_name={meta['cell']}", 
        f"--protein_name={protein_name}"
    ]
    
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = "0"
    
    print(f"\n[{dataset}] Step 2: Running SmilesGEN Molecular Optimization in {SMILESGEN_MODEL_DIR}:")
    print(" ".join(cmd_opt))
    subprocess.run(cmd_opt, cwd=SMILESGEN_MODEL_DIR, env=env, check=True)
    
    # Copy generated output to local outputs folder
    src_out_path = os.path.join(SMILESGEN_MODEL_DIR, "results", "MCF7", f"mcf7res-{protein_name}.csv")
    
    if not os.path.exists(src_out_path):
        print(f"Error: Expected output file {src_out_path} was not found.")
        return
        
    os.makedirs(OUTPUTS_DIR, exist_ok=True)
    
    # Renaming scheme
    if dataset == "Lincs1000_1":
        dest_name = "lincs1000_1.csv"
    elif dataset == "Lincs1000_2":
        dest_name = "lincs1000_2.csv"
    else:
        dest_name = f"{dataset}.csv"
        
    dest_path = os.path.join(OUTPUTS_DIR, dest_name)
    shutil.copy2(src_out_path, dest_path)
    print(f"Successfully saved optimized output to: {dest_path}")

def main():
    parser = argparse.ArgumentParser(description="Orchestrator to run SmilesGEN on test datasets")
    parser.add_argument("--dataset", choices=list(DATASET_METADATA.keys()), help="Run on a specific dataset")
    parser.add_argument("--run-all", action="store_true", help="Run all datasets sequentially")
    
    args = parser.parse_args()
    
    if not (args.dataset or args.run_all):
        parser.print_help()
        sys.exit(1)
        
    datasets = list(DATASET_METADATA.keys()) if args.run_all else [args.dataset]
    
    for d in datasets:
        setup_symlinks(d)
        run_pipeline(d)

if __name__ == "__main__":
    main()
