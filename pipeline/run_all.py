#!/usr/bin/env python3
"""
TheraFlow Master Pipeline Orchestrator (End-to-End Wrapper)
Sequentially runs all evaluation, property prediction, target profiling, 
disease reversal, consensus ranking, and docking/QM stages:

  Step 1: Generative Performance (Validity, Uniqueness, MOSES Novelty, Diversity)
  Step 2: 2D Molecular Properties, Safety Alerts (PAINS/BRENK/NIH) & Reference Similarities
  Step 3: InChIKey Deduplication, Lipinski/Complexity 3D ETKDGv3 Embedding & PoseBusters
  Step 4: ChEMBL Target Prediction, IC50/EC50 Potency & Mechanism of Action Retrieval
  Step 5: DTSEA Disease DEG Reversal Enrichment Analysis in R
  Step 6: Master Merging, 11D Pareto Multi-Objective Optimization & Robust Rank Aggregation
  Step 7: AlphaFold 3 Docking / QM Metric Extraction (Optional / If predictions exist)
  Final : Merges all calculated metrics directly back into processed_data.csv

Usage:
  python TheraFlow/run_pipeline.py --input TheraFlow/method_processing/data/processed_data.csv
"""

import os
import sys
import shutil
import argparse
import subprocess
import pandas as pd
from pathlib import Path


def log_header(title):
    print("\n" + "=" * 70)
    print(f"  {title}")
    print("=" * 70 + "\n")


def safe_merge(left_df, right_df, preferred_keys=None):
    """Robust dataframe merge avoiding column name collisions and row multiplications."""
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


def get_python_exe():
    candidates = [
        "/verona/nobackup/rojashrj/miniconda3/envs/smilesgen/bin/python",
        sys.executable
    ]
    for cand in candidates:
        if os.path.exists(cand):
            return cand
    return sys.executable

def run_command_stream(cmd, desc):
    print(f"--> Running {desc}:\n    " + " ".join(cmd) + "\n")
    res = subprocess.run(cmd)
    if res.returncode != 0:
        raise RuntimeError(f"Execution failed at stage: {desc} (exit code {res.returncode})")


def main():
    parser = argparse.ArgumentParser(
        description="TheraFlow End-to-End Master Pipeline Wrapper"
    )
    parser.add_argument(
        "--input",
        type=str,
        default="TheraFlow/method_processing/data/processed_data.csv",
        help="Path to raw candidate processed_data.csv."
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="TheraFlow/outputs",
        help="Directory to store intermediate and final metrics."
    )
    parser.add_argument(
        "--output_merged",
        type=str,
        default="TheraFlow/method_processing/data/processed_data.csv",
        help="Target path to save final merged dataset (defaults to updating processed_data.csv)."
    )
    parser.add_argument(
        "--backup",
        action="store_true",
        default=True,
        help="Create a .backup.csv copy before overwriting input dataset."
    )
    parser.add_argument(
        "--steps",
        type=str,
        default="1,2,3,4,5,6",
        help="Comma-separated list of pipeline steps to run (e.g. 1,2,3,4,5,6 or 1,2,3,4,5,6,7)."
    )
    
    # Execution options for sub-scripts
    parser.add_argument("--n_jobs_3d", type=int, default=4, help="Workers for 3D embedding in Step 3.")
    parser.add_argument("--chembl_workers", type=int, default=5, help="Workers for ChEMBL queries in Step 4.")
    parser.add_argument("--r_workers", type=int, default=2, help="Parallel workers for DTSEA in Step 5.")
    parser.add_argument("--top_n_consensus", type=int, default=10, help="Top candidates per group in Step 6.")
    parser.add_argument("--run_af3", action="store_true", help="Launch AlphaFold 3 batch inference in Step 7.")
    parser.add_argument("--extract_qm", action="store_true", help="Extract AF3 and QM metrics in Step 7.")
    
    args = parser.parse_args()
    
    input_file = Path(args.input).resolve()
    if not input_file.exists():
        raise FileNotFoundError(f"Input candidate CSV not found at: {input_file}")
        
    out_dir = Path(args.output_dir).resolve()
    os.makedirs(out_dir, exist_ok=True)
    
    pipeline_dir = Path(__file__).parent / "pipeline"
    if not pipeline_dir.exists():
        pipeline_dir = Path("TheraFlow/pipeline").resolve()
        
    active_steps = [int(s.strip()) for s in args.steps.split(",") if s.strip().isdigit()]
    print(f"Initialised TheraFlow Master Pipeline with {len(active_steps)} active steps: {active_steps}")
    print(f"Input Candidates: {input_file}")
    print(f"Outputs Directory: {out_dir}")
    
    # Output file paths
    set_props_csv = out_dir / "set_properties.csv"
    mol_props_csv = out_dir / "molecular_properties_metrics.csv"
    mol_props_filtered_csv = out_dir / "molecular_properties_filtered.csv"
    sdf_3d_file = out_dir / "filtered_3d.sdf"
    posebusters_csv = out_dir / "molecular_properties_with_posebusters.csv"
    chembl_targets_csv = out_dir / "chembl_targets.csv"
    potency_csv = out_dir / "potency_scores.csv"
    reversal_csv = out_dir / "reversal_res.csv"
    master_merged_csv = out_dir / "merged_master_dataframe.csv"
    consensus_csv = out_dir / "consensus_top.csv"
    af3_batch_json = out_dir / "af3_batch_input.json"
    docking_qm_csv = out_dir / "docking_qm_metrics.csv"
    
    # -------------------------------------------------------------------------
    # STEP 1: Generative Performance Macro Metrics
    # -------------------------------------------------------------------------
    if 1 in active_steps:
        log_header("Stage 1: Evaluating Generative Performance Macro-Metrics")
        cmd = [
            sys.executable, str(pipeline_dir / "1.generative_perf.py"),
            "--input", str(input_file),
            "--output", str(set_props_csv)
        ]
        run_command_stream(cmd, "Step 1: Generative Performance")
        
    # -------------------------------------------------------------------------
    # STEP 2: 2D Molecular Properties & Safety Alerts
    # -------------------------------------------------------------------------
    if 2 in active_steps:
        log_header("Stage 2: Calculating 2D Properties, Safety Filters & Similarities")
        cmd = [
            sys.executable, str(pipeline_dir / "2.molecular_properties.py"),
            "--input", str(input_file),
            "--output", str(mol_props_csv),
            "--output_filtered", str(mol_props_filtered_csv)
        ]
        run_command_stream(cmd, "Step 2: Molecular Properties")
        
    # -------------------------------------------------------------------------
    # STEP 3: InChIKey Deduplication, 3D Embedding & PoseBusters
    # -------------------------------------------------------------------------
    if 3 in active_steps:
        log_header("Stage 3: 3D Conformation Embedding & PoseBusters Validation")
        in_step3 = mol_props_csv if mol_props_csv.exists() else input_file
        cmd = [
            get_python_exe(), str(pipeline_dir / "3.posebuster_convert_sdf.py"),
            "--input", str(in_step3),
            "--output_sdf", str(sdf_3d_file),
            "--output_csv", str(posebusters_csv),
            "--n_jobs", str(args.n_jobs_3d)
        ]
        run_command_stream(cmd, "Step 3: PoseBusters & 3D Conformations")
        
    # -------------------------------------------------------------------------
    # STEP 4: ChEMBL Targets, Potency & Mechanism of Action
    # -------------------------------------------------------------------------
    if 4 in active_steps:
        log_header("Stage 4: Querying ChEMBL for Targets & Potency (IC50/EC50)")
        in_step4 = mol_props_filtered_csv if mol_props_filtered_csv.exists() else (mol_props_csv if mol_props_csv.exists() else input_file)
        cmd = [
            sys.executable, str(pipeline_dir / "4.potency_chembl.py"),
            "--input", str(in_step4),
            "--output_targets", str(chembl_targets_csv),
            "--output_potency", str(potency_csv),
            "--max_workers", str(args.chembl_workers)
        ]
        run_command_stream(cmd, "Step 4: ChEMBL Potency & Targets")
        
    # -------------------------------------------------------------------------
    # STEP 5: DTSEA Disease Reversal Scoring
    # -------------------------------------------------------------------------
    if 5 in active_steps:
        log_header("Stage 5: Calculating DTSEA Disease Gene Reversal in R")
        cmd = [
            sys.executable, str(pipeline_dir / "5.disease_reversal.py"),
            "--targets_csv", str(chembl_targets_csv),
            "--potency_csv", str(potency_csv),
            "--output", str(reversal_csv),
            "--workers", str(args.r_workers)
        ]
        run_command_stream(cmd, "Step 5: Disease Reversal")
        
    # -------------------------------------------------------------------------
    # STEP 6: Master Merging & Consensus Pareto / RRA Ranking
    # -------------------------------------------------------------------------
    if 6 in active_steps:
        log_header("Stage 6: Merging Metrics & Running Pareto / RRA Consensus Ranking")
        cmd = [
            sys.executable, str(pipeline_dir / "6.merging.py"),
            "--properties_csv", str(mol_props_csv if mol_props_csv.exists() else input_file),
            "--posebusters_csv", str(posebusters_csv),
            "--targets_csv", str(chembl_targets_csv),
            "--potency_csv", str(potency_csv),
            "--reversal_csv", str(reversal_csv),
            "--master_csv", str(master_merged_csv),
            "--output_consensus", str(consensus_csv),
            "--output_af3_json", str(af3_batch_json),
            "--top_n", str(args.top_n_consensus)
        ]
        run_command_stream(cmd, "Step 6: Master Merging & Consensus Ranking")
        
    # -------------------------------------------------------------------------
    # STEP 7: AlphaFold 3 Docking & QM Execution / Metric Extraction (Optional)
    # -------------------------------------------------------------------------
    if 7 in active_steps or args.run_af3 or args.extract_qm:
        log_header("Stage 7: AlphaFold 3 Docking & QM/MM Metric Extraction")
        cmd_7 = [sys.executable, str(pipeline_dir / "7.docking_md_setup.py")]
        if args.run_af3:
            cmd_7.extend(["--run_af3", "--af3_json", str(af3_batch_json)])
        if args.extract_qm or not args.run_af3:
            cmd_7.extend([
                "--extract_metrics",
                "--output_metrics", str(docking_qm_csv)
            ])
        run_command_stream(cmd_7, "Step 7: Docking & QM Metrics")
        
    # -------------------------------------------------------------------------
    # FINAL: Merge All Calculated Results Back into processed_data.csv
    # -------------------------------------------------------------------------
    log_header("Final Aggregation: Merging Annotated Metrics into Input Dataset")
    
    base_df = pd.read_csv(input_file)
    print(f"Loaded original candidates table ({len(base_df)} rows) from: {input_file}")
    
    if master_merged_csv.exists():
        print(f"Merging master metric layers from {master_merged_csv}...")
        master_df = pd.read_csv(master_merged_csv)
        final_df = safe_merge(base_df, master_df)
    else:
        final_df = base_df.copy()
        
    if consensus_csv.exists():
        print(f"Annotating top consensus ranked candidates from {consensus_csv}...")
        cons_df = pd.read_csv(consensus_csv)
        if "SMILES" in cons_df.columns:
            m_col = "method" if "method" in cons_df.columns else "Method"
            cons_keys = set(zip(cons_df["SMILES"], cons_df.get(m_col, ""), cons_df.get("test_dataset", "")))
            final_df["is_top_consensus"] = [
                (s, m, d) in cons_keys for s, m, d in zip(final_df["SMILES"], final_df.get("method", ""), final_df.get("test_dataset", ""))
            ]
            
    if docking_qm_csv.exists() and docking_qm_csv.stat().st_size > 0:
        print(f"Merging docking and QM descriptors from {docking_qm_csv}...")
        try:
            dqm_df = pd.read_csv(docking_qm_csv)
            if len(dqm_df) > 0:
                final_df = safe_merge(final_df, dqm_df, preferred_keys=["SMILES", "test_dataset"])
        except Exception as e:
            print(f"Warning: Could not merge docking QM CSV: {e}")
        
    target_out_path = Path(args.output_merged).resolve()
    
    if args.backup and target_out_path.exists():
        backup_path = target_out_path.with_name(f"{target_out_path.stem}.backup{target_out_path.suffix}")
        shutil.copy2(target_out_path, backup_path)
        print(f"Created backup of original dataset at: {backup_path}")
        
    os.makedirs(target_out_path.parent, exist_ok=True)
    final_df.to_csv(target_out_path, index=False)
    
    print(f"\n{"=" * 70}")
    print(f"SUCCESS: TheraFlow Pipeline Execution & Merging Complete!")
    print(f"Final fully annotated dataset ({len(final_df)} rows, {len(final_df.columns)} columns) saved to:")
    print(f"  --> {target_out_path}")
    print(f"{"=" * 70}\n")


if __name__ == "__main__":
    main()
