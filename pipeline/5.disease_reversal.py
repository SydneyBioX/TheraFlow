#!/usr/bin/env python3
"""
Step 5: Disease Reversal Scoring via DTSEA
Performs drug-target-set enrichment analysis against up/down disease gene signatures
and computes reversal scores scaled by drug action and IC50 potency.

Self-contained resources in TheraFlow/method_processing/data:
1. PPI Network: TheraFlow/method_processing/data/full_human_ppi.rds (auto-constructed from STRING DB v12 if missing)
2. LINCS Metadata: TheraFlow/method_processing/data/LINCS/lincs_gene_meta.csv (auto-downloaded from NCBI GEO GSE92742 if missing)
"""

import os
import urllib.request
import argparse
import subprocess
import pandas as pd
from pathlib import Path


def ensure_lincs_gene_meta(meta_path="TheraFlow/method_processing/data/LINCS/lincs_gene_meta.csv"):
    """
    Ensures that lincs_gene_meta.csv exists. If missing, downloads and parses 
    GSE92742_Broad_LINCS_gene_info.txt.gz directly from NCBI GEO.
    """
    meta_file = Path(meta_path)
    if meta_file.exists() and meta_file.stat().st_size > 1000:
        return str(meta_file.resolve())
        
    os.makedirs(meta_file.parent, exist_ok=True)
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_gene_info.txt.gz"
    print(f"LINCS metadata not found at {meta_file}. Downloading official Broad LINCS gene info from NCBI GEO:\n  {url}")
    try:
        temp_gz = str(meta_file) + ".txt.gz"
        urllib.request.urlretrieve(url, temp_gz)
        print("Parsing downloaded gene info...")
        df_info = pd.read_csv(temp_gz, sep="\t", compression="gzip", low_memory=False)
        cols_needed = ["pr_gene_id", "pr_gene_symbol", "pr_gene_title", "pr_is_lm", "pr_is_bing"]
        present_cols = [c for c in cols_needed if c in df_info.columns]
        df_info[present_cols].to_csv(meta_file, index=False)
        if os.path.exists(temp_gz):
            os.remove(temp_gz)
        print(f"Successfully generated {meta_file} ({len(df_info)} genes).")
        return str(meta_file.resolve())
    except Exception as e:
        print(f"Warning: Failed to download LINCS metadata from NCBI: {e}")
        return str(meta_file.resolve())


def run_disease_reversal_pipeline(
    targets_csv, 
    potency_csv, 
    output_reversal_csv, 
    ppi_rds="TheraFlow/method_processing/data/full_human_ppi.rds",
    lincs_meta="TheraFlow/method_processing/data/LINCS/lincs_gene_meta.csv",
    test_datasets_dir="TheraFlow/method_processing/data/test_datasets/rest_methods",
    workers=2
):
    """Executes disease_reversal.R to compute DTSEA enrichment scores across all test datasets."""
    resolved_meta = ensure_lincs_gene_meta(lincs_meta)
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    r_script_path = os.path.join(script_dir, "disease_reversal.R")
    
    if not os.path.exists(r_script_path):
        raise FileNotFoundError(f"R script not found at {r_script_path}")
        
    cmd = [
        "Rscript", r_script_path,
        "--targets_csv", str(targets_csv),
        "--potency_csv", str(potency_csv),
        "--output_csv", str(output_reversal_csv),
        "--ppi_rds", str(ppi_rds),
        "--lincs_meta", str(resolved_meta),
        "--test_datasets_dir", str(test_datasets_dir),
        "--workers", str(workers)
    ]
    
    cmd_str = " ".join(cmd)
    print(f"Launching DTSEA disease reversal R script:\n  {cmd_str}")
    result = subprocess.run(cmd, check=True)
    
    if result.returncode != 0:
        print("DTSEA execution failed!")
        print(result.stderr)
        raise RuntimeError(f"R script failed with exit code {result.returncode}")
        
    print(f"DTSEA disease reversal completed successfully: {output_reversal_csv}")
    return output_reversal_csv


def main():
    parser = argparse.ArgumentParser(description="Calculate DTSEA disease reversal scores in R.")
    parser.add_argument(
        "--targets_csv", 
        type=str, 
        default="TheraFlow/outputs/chembl_targets.csv",
        help="CSV of predicted targets."
    )
    parser.add_argument(
        "--potency_csv", 
        type=str, 
        default="TheraFlow/outputs/potency_scores.csv",
        help="CSV of potency and mechanism of action."
    )
    parser.add_argument(
        "--output", 
        type=str, 
        default="TheraFlow/outputs/reversal_res.csv",
        help="Output CSV for disease reversal results."
    )
    parser.add_argument(
        "--ppi_rds", 
        type=str, 
        default="TheraFlow/method_processing/data/full_human_ppi.rds",
        help="Path to PPI network RDS."
    )
    parser.add_argument(
        "--lincs_meta", 
        type=str, 
        default="TheraFlow/method_processing/data/LINCS/lincs_gene_meta.csv",
        help="Path to LINCS metadata."
    )
    parser.add_argument(
        "--test_datasets_dir", 
        type=str, 
        default="TheraFlow/method_processing/data/test_datasets/rest_methods",
        help="Directory of disease DEG signatures."
    )
    parser.add_argument(
        "--workers", 
        type=int, 
        default=2,
        help="Parallel R workers."
    )
    args = parser.parse_args()
    
    run_disease_reversal_pipeline(
        targets_csv=args.targets_csv,
        potency_csv=args.potency_csv,
        output_reversal_csv=args.output,
        ppi_rds=args.ppi_rds,
        lincs_meta=args.lincs_meta,
        test_datasets_dir=args.test_datasets_dir,
        workers=args.workers
    )


if __name__ == "__main__":
    main()
