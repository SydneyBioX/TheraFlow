import os
import sys
import argparse
import subprocess
import shutil
import pandas as pd
import numpy as np

# Paths setup
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
WORKSPACE_ROOT = os.path.dirname(os.path.dirname(BASE_DIR))
DATA_DIR = os.path.join(BASE_DIR, "data")
REST_METHODS_DIR = os.path.join(DATA_DIR, "test_datasets/rest_methods")
OUTPUTS_DIR = os.path.join(BASE_DIR, "outputs/gxrnn")
GXRNN_ROOT = os.path.join(WORKSPACE_ROOT, "GxRNN/GxRNN")
PYTHON_BIN = "/verona/nobackup/rojashrj/miniconda3/envs/gxrnn/bin/python"

DATASETS = [
    "lincs1000_1",
    "lincs1000_2",
    "cardiomyopathy",
    "parkinsons",
    "hcc",
    "glaucoma",
    "glaucoma_bulk"
]

def main():
    parser = argparse.ArgumentParser(description="Orchestrator to run GxRNN on test datasets")
    parser.add_argument("--dataset", choices=DATASETS, help="Run on a specific dataset")
    parser.add_argument("--run-all", action="store_true", help="Run all datasets sequentially")
    
    args = parser.parse_args()
    if not (args.dataset or args.run_all):
        parser.print_help()
        sys.exit(1)
        
    datasets = DATASETS if args.run_all else [args.dataset]
    
    # Create outputs directory
    os.makedirs(OUTPUTS_DIR, exist_ok=True)
    
    for d in datasets:
        src_path = os.path.join(REST_METHODS_DIR, f"{d}.csv")
        if not os.path.exists(src_path):
            print(f"Error: Preprocessed file {src_path} does not exist. Run preprocessing first.")
            continue
            
        print(f"\n[{d}] Preparing GxRNN input signature...")
        df = pd.read_csv(src_path)
        
        # Load reference column mappings
        ref_path = os.path.join(GXRNN_ROOT, "datasets/test_protein/AKT2.csv")
        ref = pd.read_csv(ref_path)
        
        # Add hsa: prefix and align
        if "sig_id" in df.columns:
            df = df.drop(columns=['sig_id'])
            
        df.columns = ['hsa:' + c for c in df.columns]
        df = df[[c for c in ref.columns if c in df.columns]]
        
        # Save to datasets/test_protein/{d}.csv inside GxRNN
        os.makedirs(os.path.join(GXRNN_ROOT, "datasets/test_protein"), exist_ok=True)
        tgt_in_path = os.path.join(GXRNN_ROOT, "datasets/test_protein", f"{d}.csv")
        np.savetxt(tgt_in_path, df.values, delimiter=",", header=",".join(df.columns), comments="")
        print(f"[{d}] GxRNN input prepared at {tgt_in_path}")
        
        # Run GxRNN generation
        cmd = [
            PYTHON_BIN,
            "main.py",
            "--hidden_size", "1024",
            "--generation",
            "--saved_gxrnn", "results/saved_gxrnn_hidden1024_dropout03_500.pkl",
            "--candidate_num", "1000",
            "--protein_name", d
        ]
        print(f"Executing: {' '.join(cmd)}")
        
        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = "2"
        
        subprocess.run(cmd, cwd=GXRNN_ROOT, env=env, check=True)
        print(f"[{d}] GxRNN molecule generation complete!")
        
        # Copy output to local outputs directory
        src_out_path = os.path.join(GXRNN_ROOT, f"results/generation/res-{d}.csv")
        if os.path.exists(src_out_path):
            dest_path = os.path.join(OUTPUTS_DIR, f"{d}.csv")
            df_out = pd.read_csv(src_out_path)
            # Truncate at end-of-string token $ to remove trailing generator garbage
            df_out['SMILES'] = df_out['SMILES'].astype(str).apply(lambda x: x.split('$')[0])
            df_out.to_csv(dest_path, index=False)
            print(f"[{d}] Saved output to: {dest_path}")
        else:
            print(f"Error: Expected output file {src_out_path} was not found.")

if __name__ == "__main__":
    main()
