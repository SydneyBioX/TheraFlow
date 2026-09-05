#!/usr/bin/env python3
import time
import os
import subprocess
import psutil
import pandas as pd
from pathlib import Path

def is_af3_running():
    for p in psutil.process_iter(['name', 'cmdline']):
        try:
            cmd = " ".join(p.info['cmdline'] or [])
            if "run_alphafold.py" in cmd and "af3_predictions" in cmd:
                return True
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
    return False

print("Monitoring AlphaFold 3 GPU workers...")
while is_af3_running():
    time.sleep(10)

print("All AlphaFold 3 GPU workers have finished!")

# Step 1: Extract Docking & QM Metrics
print("\nExtracting docking and confidence metrics...")
subprocess.run([
    "python", "TheraFlow/pipeline/7.docking_md_setup.py",
    "--extract_metrics",
    "--af3_output", "TheraFlow/outputs/af3_predictions",
    "--output_metrics", "TheraFlow/outputs/docking_qm_metrics.csv"
], check=True)

# Step 2: Final Safe Merge into processed_data.csv
print("\nAggregating final results back into TheraFlow/method_processing/data/processed_data.csv...")
base_csv = Path("TheraFlow/method_processing/data/processed_data.csv")
master_csv = Path("TheraFlow/outputs/merged_master_dataframe.csv")
consensus_csv = Path("TheraFlow/outputs/consensus_top.csv")
docking_csv = Path("TheraFlow/outputs/docking_qm_metrics.csv")
reversal_csv = Path("TheraFlow/outputs/reversal_res.csv")

if not base_csv.exists():
    raise FileNotFoundError(f"Base data not found: {base_csv}")

df_base = pd.read_csv(base_csv, low_memory=False)
initial_len = len(df_base)
print(f"Base processed_data.csv rows: {initial_len}")

# Merge master properties
if master_csv.exists():
    df_m = pd.read_csv(master_csv, low_memory=False)
    m_cols = [c for c in df_m.columns if c not in ["method", "test_dataset", "SMILES"]]
    if "method" in df_m.columns and "test_dataset" in df_m.columns:
        merge_keys = ["SMILES", "method", "test_dataset"]
        df_m_clean = df_m[merge_keys + m_cols].drop_duplicates(subset=merge_keys)
        cols_to_drop = [c for c in m_cols if c in df_base.columns]
        if cols_to_drop:
            df_base = df_base.drop(columns=cols_to_drop)
        df_base = df_base.merge(df_m_clean, on=merge_keys, how="left")
        print(f"Merged {len(m_cols)} columns from merged_master_dataframe.csv")

# Merge consensus ranking
if consensus_csv.exists():
    df_c = pd.read_csv(consensus_csv, low_memory=False)
    c_cols = [c for c in df_c.columns if c not in ["method", "test_dataset", "SMILES"]]
    if "method" in df_c.columns and "test_dataset" in df_c.columns:
        merge_keys = ["SMILES", "method", "test_dataset"]
        df_c_clean = df_c[merge_keys + c_cols].drop_duplicates(subset=merge_keys)
        cols_to_drop = [c for c in c_cols if c in df_base.columns]
        if cols_to_drop:
            df_base = df_base.drop(columns=cols_to_drop)
        df_base = df_base.merge(df_c_clean, on=merge_keys, how="left")
        print(f"Merged {len(c_cols)} columns from consensus_top.csv")

# Merge disease reversal
if reversal_csv.exists():
    df_r = pd.read_csv(reversal_csv, low_memory=False)
    r_cols = [c for c in df_r.columns if c not in ["method", "test_dataset", "SMILES", "method_dataset"]]
    if "method" in df_r.columns and "test_dataset" in df_r.columns:
        merge_keys = ["SMILES", "method", "test_dataset"]
        df_r_clean = df_r[merge_keys + r_cols].drop_duplicates(subset=merge_keys)
        cols_to_drop = [c for c in r_cols if c in df_base.columns]
        if cols_to_drop:
            df_base = df_base.drop(columns=cols_to_drop)
        df_base = df_base.merge(df_r_clean, on=merge_keys, how="left")
        print(f"Merged {len(r_cols)} columns from reversal_res.csv")

# Merge docking & QM metrics
if docking_csv.exists():
    df_d = pd.read_csv(docking_csv, low_memory=False)
    d_cols = [c for c in df_d.columns if c not in ["Method", "method", "test_dataset", "SMILES", "Directory"]]
    if "Method" in df_d.columns:
        df_d["method"] = df_d["Method"]
    if "method" in df_d.columns and "test_dataset" in df_d.columns:
        merge_keys = ["SMILES", "method", "test_dataset"]
        df_d_clean = df_d[merge_keys + d_cols].drop_duplicates(subset=merge_keys)
        cols_to_drop = [c for c in d_cols if c in df_base.columns]
        if cols_to_drop:
            df_base = df_base.drop(columns=cols_to_drop)
        df_base = df_base.merge(df_d_clean, on=merge_keys, how="left")
        print(f"Merged {len(d_cols)} columns from docking_qm_metrics.csv")

assert len(df_base) == initial_len, f"Row count changed from {initial_len} to {len(df_base)}!"
df_base.to_csv(base_csv, index=False)
print(f"\n=======================================================")
print(f"SUCCESS! Successfully updated {base_csv}")
print(f"Final shape: {df_base.shape} (exactly {len(df_base)} rows preserved)")
print(f"=======================================================")
