import os
import sys
import argparse
import subprocess
import shutil
import pandas as pd

# Paths setup
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
WORKSPACE_ROOT = os.path.dirname(os.path.dirname(BASE_DIR))
DATA_DIR = os.path.join(BASE_DIR, "data")
TRANSGEM_TEST_DIR = os.path.join(DATA_DIR, "test_datasets/transgem")
OUTPUTS_DIR = os.path.join(BASE_DIR, "outputs/transgem")
TRANSGEM_ROOT = os.path.join(WORKSPACE_ROOT, "TransGEM")
PYTHON_BIN = "/verona/nobackup/rojashrj/miniconda3/envs/TransGEM/bin/python"

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
    parser = argparse.ArgumentParser(description="Orchestrator to run TransGEM on test datasets")
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
        src_path = os.path.join(TRANSGEM_TEST_DIR, f"{d}.csv")
        if not os.path.exists(src_path):
            print(f"Error: Preprocessed file {src_path} does not exist. Run preprocessing first.")
            continue
            
        print(f"\n[{d}] Preparing TransGEM input dataset...")
        # Copy to TransGEM/data/test_{d}.csv
        tgt_in_path = os.path.join(TRANSGEM_ROOT, "data", f"test_{d}.csv")
        shutil.copy2(src_path, tgt_in_path)
        print(f"[{d}] TransGEM input prepared at {tgt_in_path}")
        
        # Run TransGEM generation
        cmd = [
            PYTHON_BIN,
            "app.py",
            "--data_path", "./data/",
            "--dataset", d,
            "--cell_line", "UNK",
            "--gene_encoder", "tenfold_binary",
            "--gpu", "cuda:0",
            "--seq_num", "1000"
        ]
        print(f"Executing: {' '.join(cmd)}")
        subprocess.run(cmd, cwd=TRANSGEM_ROOT, check=True)
        print(f"[{d}] TransGEM molecule generation complete!")
        
        # Format and save output
        src_out_path = os.path.join(TRANSGEM_ROOT, "result", f"result_{d}_tenfold_binary.csv")
        if os.path.exists(src_out_path):
            df_out = pd.read_csv(src_out_path)
            # Extract and rename SMILES column
            if "valid_smiles" in df_out.columns:
                df_clean = df_out[['valid_smiles']].rename(columns={'valid_smiles': 'SMILES'})
                dest_path = os.path.join(OUTPUTS_DIR, f"{d}.csv")
                df_clean.to_csv(dest_path, index=False)
                print(f"[{d}] Saved output to: {dest_path}")
            else:
                print(f"Error: Output file does not contain 'valid_smiles' column.")
        else:
            print(f"Error: Expected output file {src_out_path} was not found.")

if __name__ == "__main__":
    main()
