#!/usr/bin/env python3
"""
TheraFlow Master Pipeline Orchestrator
Executes the end-to-end evaluation and benchmarking suite across all tiers:
  Step 1: Generative Performance (Validity, Novelty, Uniqueness, IntDiv)
  Step 2: 2D Physicochemical Properties, Structural Alerts & Benchmark Filtering
  Step 3: InChIKey Deduplication, 3D Conformation Embedding (ETKDGv3) & PoseBusters Suite
  Step 4: ChEMBL Target Retrieval & Potency / Mechanism-of-Action Consensus Scoring
  Step 5: DTSEA Transcriptomic Disease Reversal Scoring (R)
  Step 6: AlphaFold 3 Docking & QM/MM Semi-empirical Metric Extraction
  Step 7: Final Master Dataset Merge & Method Standardization (processed_data.csv)
  Step 8: Multi-Tier Publication Benchmarking & Figure Generation (benchmark.R)
"""

import os
import sys
import time
import argparse
import subprocess
import pandas as pd
from pathlib import Path


def log_header(title):
    print("\n" + "=" * 80)
    print(f"  {title}")
    print("=" * 80)


def run_command_stream(cmd, step_name):
    print(f"\n[EXEC] {' '.join(cmd)}")
    start_time = time.time()
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        bufsize=1
    )
    for line in iter(process.stdout.readline, ''):
        print(f"[{step_name}] {line}", end='', flush=True)
    process.stdout.close()
    return_code = process.wait()
    elapsed = time.time() - start_time
    if return_code != 0:
        print(f"\n[ERROR] Step '{step_name}' failed with return code {return_code} after {elapsed:.1f}s.")
        raise RuntimeError(f"Step '{step_name}' failed!")
    print(f"[DONE] Step '{step_name}' completed successfully in {elapsed:.1f}s.")


def get_python_exe():
    return sys.executable


def main():
    parser = argparse.ArgumentParser(description="TheraFlow Master Benchmarking & Evaluation Pipeline Orchestrator.")
    parser.add_argument(
        "--input", 
        type=str, 
        default="TheraFlow/method_processing/data/processed_data.csv",
        help="Path to candidate molecules CSV dataset."
    )
    parser.add_argument(
        "--output_dir", 
        type=str, 
        default="TheraFlow/outputs",
        help="Directory to save pipeline outputs and intermediate files."
    )
    parser.add_argument(
        "--steps", 
        type=str, 
        default="1,2,3,4,5,6,7,8",
        help="Comma-separated list of steps to execute (1-8). Default: all."
    )
    parser.add_argument("--n_jobs_3d", type=int, default=4, help="Parallel workers for 3D ETKDGv3 embedding.")
    parser.add_argument("--chembl_workers", type=int, default=5, help="Concurrent threads for ChEMBL API.")
    parser.add_argument("--r_workers", type=int, default=7, help="Parallel R cores for DTSEA disease reversal.")
    parser.add_argument("--top_n_consensus", type=int, default=1, help="Top candidates per group for docking.")
    parser.add_argument("--run_af3", action="store_true", help="Launch AlphaFold 3 batch neural folding.")
    parser.add_argument("--extract_qm", action="store_true", help="Extract docking confidence and QM metrics.")
    args = parser.parse_args()

    active_steps = set(int(s.strip()) for s in args.steps.split(",") if s.strip().isdigit())
    pipeline_dir = Path(__file__).resolve().parent
    out_dir = Path(args.output_dir)
    os.makedirs(out_dir, exist_ok=True)
    input_file = Path(args.input)

    print(f"Starting TheraFlow Pipeline Orchestrator on: {input_file.resolve()}")
    print(f"Outputs will be stored in: {out_dir.resolve()}")
    print(f"Active stages: {sorted(list(active_steps))}")

    # Intermediate file definitions
    set_props_csv = out_dir / "set_properties.csv"
    mol_props_csv = out_dir / "molecular_properties_metrics.csv"
    mol_props_filtered_csv = out_dir / "molecular_properties_filtered.csv"
    sdf_3d_file = out_dir / "filtered_3d.sdf"
    posebusters_csv = out_dir / "molecular_properties_with_posebusters.csv"
    chembl_targets_csv = out_dir / "chembl_targets.csv"
    potency_csv = out_dir / "potency_scores.csv"
    reversal_csv = out_dir / "reversal_res.csv"
    consensus_csv = out_dir / "consensus_top.csv"
    af3_batch_json = out_dir / "af3_batch_input.json"
    docking_qm_csv = out_dir / "docking_qm_metrics.csv"
    master_merged_csv = out_dir / "merged_master_dataframe.csv"

    # -------------------------------------------------------------------------
    # STEP 1: Generative Performance
    # -------------------------------------------------------------------------
    if 1 in active_steps:
        log_header("Stage 1: Calculating Generative Performance (Set Properties)")
        cmd = [
            get_python_exe(), str(pipeline_dir / "1.generative_perf.py"),
            "--input", str(input_file),
            "--output", str(set_props_csv)
        ]
        run_command_stream(cmd, "Step 1: Generative Performance")

    # -------------------------------------------------------------------------
    # STEP 2: Molecular Properties & Safety Alerts
    # -------------------------------------------------------------------------
    if 2 in active_steps:
        log_header("Stage 2: Calculating 2D Physicochemical Properties & Safety Filtering")
        cmd = [
            get_python_exe(), str(pipeline_dir / "2.molecular_properties.py"),
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
    # STEP 4: ChEMBL Targets & Potency
    # -------------------------------------------------------------------------
    if 4 in active_steps:
        log_header("Stage 4: Querying ChEMBL for Targets & Potency (IC50/EC50)")
        in_step4 = mol_props_filtered_csv if mol_props_filtered_csv.exists() else (mol_props_csv if mol_props_csv.exists() else input_file)
        cmd = [
            get_python_exe(), str(pipeline_dir / "4.potency_chembl.py"),
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
    # STEP 6: AlphaFold 3 Docking & QM Execution / Metric Extraction
    # -------------------------------------------------------------------------
    if 6 in active_steps:
        log_header("Stage 6: AlphaFold 3 Docking & QM Metric Extraction")
        cmd_6 = [sys.executable, str(pipeline_dir / "6.docking_md.py")]
        if args.run_af3:
            cmd_6.extend(["--run_af3", "--af3_json", str(af3_batch_json)])
        if args.extract_qm or not args.run_af3:
            cmd_6.extend([
                "--extract_metrics",
                "--output_metrics", str(docking_qm_csv)
            ])
        run_command_stream(cmd_6, "Step 6: Docking & QM Metrics")

    # -------------------------------------------------------------------------
    # STEP 7: Final Master Merge & Method Standardization
    # -------------------------------------------------------------------------
    if 7 in active_steps:
        log_header("Stage 7: Final Master Dataset Merge & Method Standardization")
        cmd_7 = [
            sys.executable, str(pipeline_dir / "7.final_merge.py"),
            "--base_csv", str(input_file),
            "--properties_csv", str(mol_props_csv),
            "--posebusters_csv", str(posebusters_csv),
            "--targets_csv", str(chembl_targets_csv),
            "--potency_csv", str(potency_csv),
            "--reversal_csv", str(reversal_csv),
            "--consensus_csv", str(consensus_csv),
            "--docking_csv", str(docking_qm_csv),
            "--master_csv", str(master_merged_csv),
            "--output_processed_data", str(input_file)
        ]
        run_command_stream(cmd_7, "Step 7: Final Master Merge")

    # -------------------------------------------------------------------------
    # STEP 8: Multi-Tier Benchmarking & Figure Generation (benchmark.R)
    # -------------------------------------------------------------------------
    if 8 in active_steps:
        log_header("Stage 8: Multi-Tier Benchmarking & Figure Generation in R")
        figures_dir = out_dir / "figures"
        os.makedirs(figures_dir, exist_ok=True)
        cmd_8 = [
            "Rscript", str(pipeline_dir / "benchmark.R"),
            "--master_csv", str(master_merged_csv),
            "--set_properties_csv", str(set_props_csv),
            "--consensus_top_csv", str(consensus_csv),
            "--potency_scores_csv", str(potency_csv),
            "--reversal_res_csv", str(reversal_csv),
            "--docking_qm_metrics_csv", str(docking_qm_csv),
            "--architecture_csv", "TheraFlow/method_processing/data/architecture_complexity.csv",
            "--figures_dir", str(figures_dir)
        ]
        run_command_stream(cmd_8, "Step 8: Benchmarking & Figures")

    log_header("TheraFlow Pipeline Execution Finished Successfully!")


if __name__ == "__main__":
    main()
