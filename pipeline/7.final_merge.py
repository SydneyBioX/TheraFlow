#!/usr/bin/env python3
"""
Step 7: Final Master Dataset Merge & Method Standardization
Aggregates all calculated metrics from Steps 1-6 (generative performance, 2D properties,
PoseBusters 3D structural validation via canonical SMILES, ChEMBL targets, potency,
DTSEA disease reversal, consensus rankings, and AF3 docking/xTB QM metrics) into a unified master dataset.
Also standardizes method nomenclature across all outputs.
"""

import os
import argparse
import pandas as pd
from pathlib import Path
from rdkit import Chem


# Method naming standardization mapping
METHOD_NAME_MAPPING = {
    "biaee": "BiAAE",
    "biaae": "BiAAE",
    "bicev": "BiCEV",
    "fame": "FAME",
    "gex2sgen": "Gex2SGen",
    "gxrnn": "GxRNN",
    "paccmann_rl": "Paccmann RL",
    "smilesgen": "SmilesGEN",
    "smilesgen (vae)": "SmilesGEN",
    "transgem": "TransGEM",
    "triomphe": "TRIOMPHE",
}

def get_canonical_smiles(s):
    """Safely converts SMILES to canonical form for robust joins."""
    if pd.isna(s) or not str(s).strip():
        return None
    s_clean = str(s).strip().replace("[C@@H1]", "[C@@H]").replace("[C@H1]", "[C@H]").replace("[C@@1]", "[C@@]").replace("[C@1]", "[C@]")
    try:
        mol = Chem.MolFromSmiles(s_clean)
        return Chem.MolToSmiles(mol, canonical=True) if mol else s_clean
    except Exception:
        return s_clean


def run_final_merge(
    base_csv="TheraFlow/method_processing/data/processed_data.csv",
    properties_csv="TheraFlow/outputs/molecular_properties_metrics.csv",
    posebusters_csv="TheraFlow/outputs/molecular_properties_with_posebusters.csv",
    targets_csv="TheraFlow/outputs/chembl_targets.csv",
    potency_csv="TheraFlow/outputs/potency_scores.csv",
    reversal_csv="TheraFlow/outputs/reversal_res.csv",
    consensus_csv="TheraFlow/outputs/consensus_top.csv",
    docking_csv="TheraFlow/outputs/docking_qm_metrics.csv",
    master_csv="TheraFlow/outputs/merged_master_dataframe.csv",
    output_processed_data="TheraFlow/method_processing/data/processed_data.csv"
):
    print(f"Reading base candidate dataset from: {base_csv}")
    base_path = Path(base_csv)
    if not base_path.exists():
        raise FileNotFoundError(f"Base candidate dataset not found: {base_csv}")

    df_base = pd.read_csv(base_path, low_memory=False)
    initial_len = len(df_base)
    print(f"Loaded {initial_len} total rows.")

    # Canonicalize SMILES for joins
    print("Computing canonical SMILES for exact cross-pipeline alignment...")
    df_base["canon_smi"] = df_base["SMILES"].apply(get_canonical_smiles)

    # 1. Merge Step 2 Molecular Properties
    if os.path.exists(properties_csv):
        print(f"1. Merging 2D Molecular Properties from {properties_csv}...")
        df_props = pd.read_csv(properties_csv, low_memory=False)
        df_props["canon_smi"] = df_props["SMILES"].apply(get_canonical_smiles)
        exclude = ["SMILES", "canon_smi", "method", "test_dataset", "Unnamed: 0"]
        prop_cols = [c for c in df_props.columns if c not in exclude]
        
        merge_keys = ["canon_smi"]
        if "test_dataset" in df_props.columns and "test_dataset" in df_base.columns:
            merge_keys = ["canon_smi", "test_dataset"]
            
        df_props_dedup = df_props.drop_duplicates(subset=merge_keys)[merge_keys + prop_cols]
        cols_to_drop = [c for c in prop_cols if c in df_base.columns]
        if cols_to_drop:
            df_base = df_base.drop(columns=cols_to_drop)
        df_base = df_base.merge(df_props_dedup, on=merge_keys, how="left")
        print(f"   -> Merged {len(prop_cols)} molecular property columns.")

    # 2. Merge Step 3 PoseBusters Validation (using canonical SMILES)
    if os.path.exists(posebusters_csv):
        print(f"2. Merging PoseBusters Structural Validation from {posebusters_csv}...")
        df_pb = pd.read_csv(posebusters_csv, low_memory=False)
        df_pb["canon_smi"] = df_pb["SMILES"].apply(get_canonical_smiles)
        exclude = ["mol_index", "SMILES", "canon_smi", "Unnamed: 0"]
        pb_cols = [c for c in df_pb.columns if c not in exclude]
        
        df_pb_dedup = df_pb.drop_duplicates(subset=["canon_smi"])[["canon_smi"] + pb_cols]
        cols_to_drop = [c for c in pb_cols if c in df_base.columns]
        if cols_to_drop:
            df_base = df_base.drop(columns=cols_to_drop)
        df_base = df_base.merge(df_pb_dedup, on="canon_smi", how="left")
        print(f"   -> Merged {len(pb_cols)} PoseBusters columns ({df_base['sanitization'].notna().sum()} rows matched).")

    # 3. Merge Step 4 ChEMBL Targets
    if os.path.exists(targets_csv):
        print(f"3. Merging ChEMBL Targets from {targets_csv}...")
        df_targets = pd.read_csv(targets_csv, low_memory=False)
        df_targets["canon_smi"] = df_targets["SMILES"].apply(get_canonical_smiles)
        exclude = ["SMILES", "canon_smi", "Unnamed: 0"]
        t_cols = [c for c in df_targets.columns if c not in exclude]
        
        df_t_dedup = df_targets.drop_duplicates(subset=["canon_smi"])[["canon_smi"] + t_cols]
        cols_to_drop = [c for c in t_cols if c in df_base.columns]
        if cols_to_drop:
            df_base = df_base.drop(columns=cols_to_drop)
        df_base = df_base.merge(df_t_dedup, on="canon_smi", how="left")
        print(f"   -> Merged {len(t_cols)} target columns.")

    # 4. Merge Step 4 Potency & Mechanism of Action
    if os.path.exists(potency_csv):
        print(f"4. Merging ChEMBL Potency Scores from {potency_csv}...")
        df_pot = pd.read_csv(potency_csv, low_memory=False)
        smi_col = "query_smiles" if "query_smiles" in df_pot.columns else "SMILES"
        if smi_col in df_pot.columns:
            df_pot["canon_smi"] = df_pot[smi_col].apply(get_canonical_smiles)
            exclude = [smi_col, "canon_smi", "Unnamed: 0"]
            pot_cols = [c for c in df_pot.columns if c not in exclude]
            
            df_pot_dedup = df_pot.drop_duplicates(subset=["canon_smi"])[["canon_smi"] + pot_cols]
            cols_to_drop = [c for c in pot_cols if c in df_base.columns]
            if cols_to_drop:
                df_base = df_base.drop(columns=cols_to_drop)
            df_base = df_base.merge(df_pot_dedup, on="canon_smi", how="left")
            print(f"   -> Merged {len(pot_cols)} potency columns.")

    # 5. Merge Step 5 Disease Reversal Results
    if os.path.exists(reversal_csv):
        print(f"5. Merging DTSEA Disease Reversal from {reversal_csv}...")
        df_rev = pd.read_csv(reversal_csv, low_memory=False)
        df_rev["canon_smi"] = df_rev["SMILES"].apply(get_canonical_smiles)
        exclude = ["SMILES", "canon_smi", "method", "test_dataset", "method_dataset", "Unnamed: 0"]
        rev_cols = [c for c in df_rev.columns if c not in exclude]
        
        merge_keys = ["canon_smi"]
        if "test_dataset" in df_rev.columns and "test_dataset" in df_base.columns:
            merge_keys = ["canon_smi", "test_dataset"]
            
        df_rev_dedup = df_rev.drop_duplicates(subset=merge_keys)[merge_keys + rev_cols]
        cols_to_drop = [c for c in rev_cols if c in df_base.columns]
        if cols_to_drop:
            df_base = df_base.drop(columns=cols_to_drop)
        df_base = df_base.merge(df_rev_dedup, on=merge_keys, how="left")
        print(f"   -> Merged {len(rev_cols)} disease reversal columns.")

    # 6. Merge Consensus Ranking Top Selections
    if os.path.exists(consensus_csv):
        print(f"6. Merging Consensus Ranking from {consensus_csv}...")
        df_cons = pd.read_csv(consensus_csv, low_memory=False)
        df_cons["canon_smi"] = df_cons["SMILES"].apply(get_canonical_smiles)
        exclude = ["SMILES", "canon_smi", "method", "test_dataset", "target_sequence", "Unnamed: 0"]
        cons_cols = [c for c in df_cons.columns if c not in exclude]
        
        merge_keys = ["canon_smi"]
        if "test_dataset" in df_cons.columns and "test_dataset" in df_base.columns:
            merge_keys = ["canon_smi", "test_dataset"]
            
        df_cons_dedup = df_cons.drop_duplicates(subset=merge_keys)[merge_keys + cons_cols]
        cols_to_drop = [c for c in cons_cols if c in df_base.columns]
        if cols_to_drop:
            df_base = df_base.drop(columns=cols_to_drop)
        df_base = df_base.merge(df_cons_dedup, on=merge_keys, how="left")
        print(f"   -> Merged {len(cons_cols)} consensus ranking columns.")

    # 7. Merge Step 6 AlphaFold 3 Docking & xTB QM Metrics
    if os.path.exists(docking_csv):
        print(f"7. Merging AF3 Docking & xTB QM Metrics from {docking_csv}...")
        df_dock = pd.read_csv(docking_csv, low_memory=False)
        if "SMILES" in df_dock.columns:
            df_dock["canon_smi"] = df_dock["SMILES"].apply(get_canonical_smiles)
            exclude = ["SMILES", "canon_smi", "method", "Method", "test_dataset", "Directory", "Model_Index", "Unnamed: 0"]
            dock_cols = [c for c in df_dock.columns if c not in exclude]
            
            merge_keys = ["canon_smi"]
            if "test_dataset" in df_dock.columns and "test_dataset" in df_base.columns:
                merge_keys = ["canon_smi", "test_dataset"]
                
            df_dock_dedup = df_dock.drop_duplicates(subset=merge_keys)[merge_keys + dock_cols]
            cols_to_drop = [c for c in dock_cols if c in df_base.columns]
            if cols_to_drop:
                df_base = df_base.drop(columns=cols_to_drop)
            df_base = df_base.merge(df_dock_dedup, on=merge_keys, how="left")
            print(f"   -> Merged {len(dock_cols)} docking and xTB QM metric columns.")

    # Drop temporary join key
    if "canon_smi" in df_base.columns:
        df_base = df_base.drop(columns=["canon_smi"])

    # -------------------------------------------------------------------------
    # 8. Standardize Method Names
    # -------------------------------------------------------------------------
    print("\nStandardizing method names across dataset...")
    if "method" in df_base.columns:
        df_base["method"] = df_base["method"].apply(lambda m: METHOD_NAME_MAPPING.get(str(m).strip().lower(), m) if pd.notna(m) else m)
        print("Updated method names distribution:")
        print(df_base["method"].value_counts())

    # -------------------------------------------------------------------------
    # Row Integrity Assertion & Saving
    # -------------------------------------------------------------------------
    assert len(df_base) == initial_len, f"Integrity error: row count changed from {initial_len} to {len(df_base)}!"
    
    out_proc = Path(output_processed_data)
    os.makedirs(out_proc.parent, exist_ok=True)
    df_base.to_csv(out_proc, index=False)
    print(f"\n[Saved] Processed data saved to: {out_proc.resolve()} (Shape: {df_base.shape})")

    out_master = Path(master_csv)
    os.makedirs(out_master.parent, exist_ok=True)
    df_base.to_csv(out_master, index=False)
    print(f"[Saved] Master merged dataset saved to: {out_master.resolve()} (Shape: {df_base.shape})")
    print(f"\n=======================================================")
    print(f"SUCCESS: Final merge complete! {len(df_base)} rows, {len(df_base.columns)} columns.")
    print(f"=======================================================")
    return str(out_proc)


def main():
    parser = argparse.ArgumentParser(description="Final master dataset merge and method name standardization.")
    parser.add_argument("--base_csv", type=str, default="TheraFlow/method_processing/data/processed_data.csv")
    parser.add_argument("--properties_csv", type=str, default="TheraFlow/outputs/molecular_properties_metrics.csv")
    parser.add_argument("--posebusters_csv", type=str, default="TheraFlow/outputs/molecular_properties_with_posebusters.csv")
    parser.add_argument("--targets_csv", type=str, default="TheraFlow/outputs/chembl_targets.csv")
    parser.add_argument("--potency_csv", type=str, default="TheraFlow/outputs/potency_scores.csv")
    parser.add_argument("--reversal_csv", type=str, default="TheraFlow/outputs/reversal_res.csv")
    parser.add_argument("--consensus_csv", type=str, default="TheraFlow/outputs/consensus_top.csv")
    parser.add_argument("--docking_csv", type=str, default="TheraFlow/outputs/docking_qm_metrics.csv")
    parser.add_argument("--master_csv", type=str, default="TheraFlow/outputs/merged_master_dataframe.csv")
    parser.add_argument("--output_processed_data", type=str, default="TheraFlow/method_processing/data/processed_data.csv")
    args = parser.parse_args()

    run_final_merge(
        base_csv=args.base_csv,
        properties_csv=args.properties_csv,
        posebusters_csv=args.posebusters_csv,
        targets_csv=args.targets_csv,
        potency_csv=args.potency_csv,
        reversal_csv=args.reversal_csv,
        consensus_csv=args.consensus_csv,
        docking_csv=args.docking_csv,
        master_csv=args.master_csv,
        output_processed_data=args.output_processed_data
    )


if __name__ == "__main__":
    main()
