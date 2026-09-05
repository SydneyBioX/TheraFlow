#!/usr/bin/env python3
"""
Step 6: Master Metrics Merging & Consensus Ranking Execution
Merges all calculated property and test layers into a unified master dataframe,
and triggers consensus_ranking.R (from metrics_calculation/7.2. result5.qmd)
for Pareto optimization, RRA scoring, and AlphaFold 3 batch preparation.
"""

import os
import argparse
import subprocess
import pandas as pd
from pathlib import Path


def safe_merge(left_df, right_df, preferred_keys=None):
    """Merges right_df into left_df on preferred keys, avoiding duplicate columns and row expansions."""
    if right_df is None or len(right_df) == 0:
        return left_df
    if preferred_keys is None:
        preferred_keys = ["SMILES", "method", "test_dataset"]
    merge_keys = [k for k in preferred_keys if k in left_df.columns and k in right_df.columns]
    if not merge_keys and "SMILES" in left_df.columns and "SMILES" in right_df.columns:
        merge_keys = ["SMILES"]
    if not merge_keys:
        return left_df
    cols_to_add = [c for c in right_df.columns if c not in left_df.columns or c in merge_keys]
    right_dedup = right_df[cols_to_add].drop_duplicates(subset=merge_keys)
    return left_df.merge(right_dedup, on=merge_keys, how="left")


def run_merging_pipeline(
    properties_csv="TheraFlow/outputs/molecular_properties_metrics.csv",
    posebusters_csv="TheraFlow/outputs/molecular_properties_with_posebusters.csv",
    targets_csv="TheraFlow/outputs/chembl_targets.csv",
    potency_csv="TheraFlow/outputs/potency_scores.csv",
    reversal_csv="TheraFlow/outputs/reversal_res.csv",
    master_csv="TheraFlow/outputs/merged_master_dataframe.csv",
    output_consensus="TheraFlow/outputs/consensus_top.csv",
    output_af3_json="TheraFlow/outputs/af3_batch_input.json",
    top_n=10
):
    """Merges all metrics and triggers the consensus ranking R script."""
    print(f"Reading base properties from: {properties_csv}")
    if not os.path.exists(properties_csv):
        raise FileNotFoundError(f"Base properties CSV not found at: {properties_csv}")
        
    master_df = pd.read_csv(properties_csv)
    
    if os.path.exists(posebusters_csv):
        print(f"  Merging PoseBusters from {posebusters_csv}...")
        pb_df = pd.read_csv(posebusters_csv)
        master_df = safe_merge(master_df, pb_df)
    else:
        print(f"  [Note] PoseBusters CSV not found at {posebusters_csv}, skipping.")
        
    if os.path.exists(targets_csv):
        print(f"  Merging ChEMBL targets from {targets_csv}...")
        targets_df = pd.read_csv(targets_csv)
        master_df = safe_merge(master_df, targets_df, ["SMILES"])
    else:
        print(f"  [Note] Targets CSV not found at {targets_csv}, skipping.")
        
    if os.path.exists(potency_csv):
        print(f"  Merging potency from {potency_csv}...")
        pot_df = pd.read_csv(potency_csv)
        if "query_smiles" in pot_df.columns:
            pot_df = pot_df.rename(columns={"query_smiles": "SMILES"})
        master_df = safe_merge(master_df, pot_df, ["SMILES"])
    else:
        print(f"  [Note] Potency CSV not found at {potency_csv}, skipping.")
        
    if os.path.exists(reversal_csv):
        print(f"  Merging disease reversal from {reversal_csv}...")
        rev_df = pd.read_csv(reversal_csv)
        master_df = safe_merge(master_df, rev_df, ["SMILES", "test_dataset"])
    else:
        print(f"  [Note] Disease reversal CSV not found at {reversal_csv}, skipping.")
        
    out_master = Path(master_csv)
    os.makedirs(out_master.parent if str(out_master.parent) != "" else ".", exist_ok=True)
    master_df.to_csv(out_master, index=False)
    print(f"Successfully saved merged master dataframe ({len(master_df)} rows) to: {out_master.resolve()}")
    
    # Run consensus ranking in R
    script_dir = os.path.dirname(os.path.abspath(__file__))
    r_ranking_script = os.path.join(script_dir, "consensus_ranking.R")
    
    if os.path.exists(r_ranking_script):
        cmd = [
            "Rscript", r_ranking_script,
            "--master_csv", str(out_master),
            "--output_consensus", str(output_consensus),
            "--output_af3_json", str(output_af3_json),
            "--top_n", str(top_n)
        ]
        cmd_str = " ".join(cmd)
        print(f"Launching consensus ranking R script:\n  {cmd_str}")
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            print("Error running consensus ranking R script:")
            print(res.stderr)
            raise RuntimeError(f"consensus_ranking.R failed with exit code {res.returncode}")
        print(res.stdout)
        print("Consensus ranking pipeline completed successfully!")
    else:
        print(f"Warning: consensus_ranking.R not found at {r_ranking_script}")
        
    return str(out_master)


def main():
    parser = argparse.ArgumentParser(description="Merge calculated metrics and execute consensus ranking.")
    parser.add_argument(
        "--properties_csv", 
        type=str, 
        default="TheraFlow/outputs/molecular_properties_metrics.csv",
        help="Calculated molecular properties CSV."
    )
    parser.add_argument(
        "--posebusters_csv", 
        type=str, 
        default="TheraFlow/outputs/molecular_properties_with_posebusters.csv",
        help="PoseBusters test results CSV."
    )
    parser.add_argument(
        "--targets_csv", 
        type=str, 
        default="TheraFlow/outputs/chembl_targets.csv",
        help="ChEMBL predicted targets CSV."
    )
    parser.add_argument(
        "--potency_csv", 
        type=str, 
        default="TheraFlow/outputs/potency_scores.csv",
        help="Consensus potency CSV."
    )
    parser.add_argument(
        "--reversal_csv", 
        type=str, 
        default="TheraFlow/outputs/reversal_res.csv",
        help="DTSEA disease reversal CSV."
    )
    parser.add_argument(
        "--master_csv", 
        type=str, 
        default="TheraFlow/outputs/merged_master_dataframe.csv",
        help="Output merged master dataframe CSV."
    )
    parser.add_argument(
        "--output_consensus", 
        type=str, 
        default="TheraFlow/outputs/consensus_top.csv",
        help="Output top consensus candidates CSV."
    )
    parser.add_argument(
        "--output_af3_json", 
        type=str, 
        default="TheraFlow/outputs/af3_batch_input.json",
        help="Output AlphaFold 3 batch input JSON."
    )
    parser.add_argument(
        "--top_n", 
        type=int, 
        default=10,
        help="Number of top candidates to pick per group."
    )
    args = parser.parse_args()
    
    run_merging_pipeline(
        properties_csv=args.properties_csv,
        posebusters_csv=args.posebusters_csv,
        targets_csv=args.targets_csv,
        potency_csv=args.potency_csv,
        reversal_csv=args.reversal_csv,
        master_csv=args.master_csv,
        output_consensus=args.output_consensus,
        output_af3_json=args.output_af3_json,
        top_n=args.top_n
    )


if __name__ == "__main__":
    main()
