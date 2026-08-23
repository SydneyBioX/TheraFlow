#!/usr/bin/env python
import os
import argparse
import json
import pandas as pd
from pathlib import Path

# Import pipeline steps
from utils import merge_outputs, filter_for_posebusters
from molecular_properties import calculate_properties_pipeline
from posebuster_analysis import run_posebusters
from potency_chembl import run_potency_chembl_pipeline
from disease_reversal import run_disease_reversal_pipeline
from consensus_ranking import run_consensus_ranking_pipeline


def main():
    parser = argparse.ArgumentParser(description="Omics2Drugs Metrics Calculation Pipeline")
    
    # Path configuration arguments
    parser.add_argument("--method_outputs_dir", type=str, default="../../results/method_outputs", help="Directory of method raw outputs")
    parser.add_argument("--output_dir", type=str, default="../../results/pipeline_runs", help="Output directory for pipeline files")
    
    # R script dependencies
    parser.add_argument("--ppi_rds", type=str, default="data/ppi.rds", help="Path to PPI rds network file (originally from STRING DB, filtered for experimental confidence score > 700)")
    parser.add_argument("--mapping_file", type=str, default="data/chembl_uniprot_mapping.txt", help="Path to ChEMBL-UniProt mapping")
    parser.add_argument("--lincs_meta", type=str, default="Omics2Drug/data/lincs_gene_meta.csv", help="Path to LINCS gene meta")
    parser.add_argument("--test_datasets_dir", type=str, default="Omics2Drug/data/test_datasets", help="Directory of disease signatures")
    
    # Pipeline execution arguments
    parser.add_argument("--steps", type=str, default="1,2,3,4,5,6,7,8", help="Steps to run (e.g. '1,2,5')")
    parser.add_argument("--n_jobs_3d", type=int, default=4, help="Parallel jobs for RDKit 3D conformation generation")
    parser.add_argument("--workers_r", type=int, default=4, help="Parallel workers for R DTSEA calculation")
    parser.add_argument("--workers_chembl", type=int, default=5, help="Thread workers for ChEMBL API client queries")
    
    args = parser.parse_args()
    
    # Normalize paths relative to workspace root
    base_dir = Path(os.getcwd())
    method_outputs_dir = base_dir / args.method_outputs_dir
    output_dir = base_dir / args.output_dir
    
    ppi_rds = base_dir / args.ppi_rds
    mapping_file = base_dir / args.mapping_file
    lincs_meta = base_dir / args.lincs_meta
    test_datasets_dir = base_dir / args.test_datasets_dir
    
    # Define filenames for each pipeline step
    step1_csv = output_dir / "step1_merged_outputs.csv"
    step2_csv = output_dir / "step2_molecular_properties.csv"
    step2_summary_csv = output_dir / "step2_generative_summary.csv"
    step3_sdf = output_dir / "step3_conformations.sdf"
    step4_csv = output_dir / "step4_posebusters.csv"
    step5_targets_csv = output_dir / "step5_chembl_targets.csv"
    step5_potency_csv = output_dir / "step5_chembl_potency.csv"
    step6_reversal_csv = output_dir / "step6_disease_reversal.csv"
    step7_master_csv = output_dir / "step7_merged_master_dataframe.csv"
    step8_consensus_csv = output_dir / "step8_consensus_top10.csv"
    
    steps_to_run = [int(s.strip()) for s in args.steps.split(",")]
    
    os.makedirs(output_dir, exist_ok=True)
    print(f"--- Launching Metrics calculation pipeline in: {output_dir} ---")
    
    # -------------------------------------------------------------------------
    # Step 1: Merge raw outputs
    # -------------------------------------------------------------------------
    if 1 in steps_to_run:
        print("\n=== STEP 1: Merging raw outputs ===")
        merge_outputs(method_outputs_dir, step1_csv)
    
    # -------------------------------------------------------------------------
    # Step 2: Calculate Molecular Properties
    # -------------------------------------------------------------------------
    if 2 in steps_to_run:
        print("\n=== STEP 2: Calculating molecular properties ===")
        if not step1_csv.exists():
            raise FileNotFoundError(f"Step 1 output missing at {step1_csv}. Run Step 1 first.")
        calculate_properties_pipeline(step1_csv, step2_csv, step2_summary_csv)
        
    # -------------------------------------------------------------------------
    # Step 3: Filter and embed 3D conformations
    # -------------------------------------------------------------------------
    if 3 in steps_to_run:
        print("\n=== STEP 3: Filtering & generating 3D conformations ===")
        if not step2_csv.exists():
            raise FileNotFoundError(f"Step 2 output missing at {step2_csv}. Run Step 2 first.")
        filter_for_posebusters(step2_csv, step3_sdf, n_jobs=args.n_jobs_3d)
        
    # -------------------------------------------------------------------------
    # Step 4: PoseBusters validation
    # -------------------------------------------------------------------------
    if 4 in steps_to_run:
        print("\n=== STEP 4: Running PoseBusters validation ===")
        if not step3_sdf.exists():
            raise FileNotFoundError(f"Step 3 output missing at {step3_sdf}. Run Step 3 first.")
        run_posebusters(step3_sdf, step4_csv)
        
    # -------------------------------------------------------------------------
    # Step 5: ChEMBL API potency lookup
    # -------------------------------------------------------------------------
    if 5 in steps_to_run:
        print("\n=== STEP 5: Querying ChEMBL targets and potency ===")
        if not step2_csv.exists():
            raise FileNotFoundError(f"Step 2 output missing at {step2_csv}. Run Step 2 first.")
        run_potency_chembl_pipeline(
            properties_csv=step2_csv, 
            output_targets_csv=step5_targets_csv, 
            output_potency_csv=step5_potency_csv, 
            max_workers=args.workers_chembl
        )
        
    # -------------------------------------------------------------------------
    # Step 6: Disease reversal in R (DTSEA)
    # -------------------------------------------------------------------------
    if 6 in steps_to_run:
        print("\n=== STEP 6: Executing DTSEA disease reversal in R ===")
        if not step5_targets_csv.exists() or not step5_potency_csv.exists():
            raise FileNotFoundError("Step 5 outputs missing. Run Step 5 first.")
            
        # Merge step2 properties metadata (method, test_dataset) with predicted targets for Step 6
        print("Merging molecular properties metadata with predicted targets for R DTSEA...")
        props_df = pd.read_csv(step2_csv)
        targets_df = pd.read_csv(step5_targets_csv)
        merged_targets = props_df[['SMILES', 'method', 'test_dataset']].merge(targets_df, on='SMILES', how='inner')
        step5_targets_merged_csv = output_dir / "step5_chembl_targets_merged.csv"
        merged_targets.to_csv(step5_targets_merged_csv, index=False)
        
        run_disease_reversal_pipeline(
            targets_csv=step5_targets_merged_csv, 
            potency_csv=step5_potency_csv, 
            output_reversal_csv=step6_reversal_csv, 
            ppi_rds=ppi_rds, 
            mapping_file=mapping_file, 
            lincs_meta=lincs_meta, 
            test_datasets_dir=test_datasets_dir, 
            workers=args.workers_r
        )
        
    # -------------------------------------------------------------------------
    # Step 7: Final merge of outputs
    # -------------------------------------------------------------------------
    if 7 in steps_to_run:
        print("\n=== STEP 7: Merging outputs into master dataframe ===")
        if not step2_csv.exists():
            raise FileNotFoundError(f"Step 2 molecular properties CSV missing at {step2_csv}.")
            
        print("Loading molecular properties base...")
        master_df = pd.read_csv(step2_csv)
        
        # Merge Step 4 PoseBusters
        if step4_csv.exists():
            print("Merging PoseBusters validation results...")
            pb_df = pd.read_csv(step4_csv)
            pb_cols = [c for c in pb_df.columns if c not in master_df.columns or c in ['SMILES', 'method', 'test_dataset']]
            master_df = master_df.merge(pb_df[pb_cols], on=['SMILES', 'method', 'test_dataset'], how='left')
            
        # Merge Step 5 targets
        if step5_targets_csv.exists():
            print("Merging ChEMBL predicted targets...")
            targets_df = pd.read_csv(step5_targets_csv)
            master_df = master_df.merge(targets_df, on='SMILES', how='left')
            
        # Merge Step 5 potency
        if step5_potency_csv.exists():
            print("Merging consensus potency scores...")
            potency_df = pd.read_csv(step5_potency_csv)
            if 'query_smiles' in potency_df.columns:
                potency_df = potency_df.rename(columns={'query_smiles': 'SMILES'})
            potency_cols = [c for c in potency_df.columns if c not in master_df.columns or c == 'SMILES']
            master_df = master_df.merge(potency_df[potency_cols], on='SMILES', how='left')
            
        # Merge Step 6 disease reversal
        if step6_reversal_csv.exists():
            print("Merging disease reversal DTSEA scores...")
            reversal_df = pd.read_csv(step6_reversal_csv)
            
            if 'method_dataset' in reversal_df.columns and 'method' not in reversal_df.columns:
                if 'method_dataset' not in master_df.columns:
                    master_df['method_dataset'] = master_df['method'] + '_' + master_df['test_dataset']
                reversal_cols = [c for c in reversal_df.columns if c not in master_df.columns or c in ['SMILES', 'method_dataset']]
                master_df = master_df.merge(reversal_df[reversal_cols], on=['SMILES', 'method_dataset'], how='left')
            else:
                reversal_cols = [c for c in reversal_df.columns if c not in master_df.columns or c in ['SMILES', 'method', 'test_dataset']]
                master_df = master_df.merge(reversal_df[reversal_cols], on=['SMILES', 'method', 'test_dataset'], how='left')
                
        # Write master dataframe
        master_df.to_csv(step7_master_csv, index=False)
        print(f"Master dataframe successfully generated: {step7_master_csv}")
        print(f"Total entries in final master: {len(master_df)}")
        
    # -------------------------------------------------------------------------
    # Step 8: Pareto and Robust Rank Aggregation (RRA) Consensus Ranking
    # -------------------------------------------------------------------------
    if 8 in steps_to_run:
        print("\n=== STEP 8: Performing Pareto & Robust Rank Aggregation (RRA) Consensus Ranking ===")
        if not step7_master_csv.exists():
            raise FileNotFoundError(f"Step 7 master dataframe missing at {step7_master_csv}. Run Step 7 first.")
        run_consensus_ranking_pipeline(
            input_csv=step7_master_csv, 
            output_csv=step8_consensus_csv
        )

    # -------------------------------------------------------------------------
    

if __name__ == "__main__":
    main()
