#!/usr/bin/env python3
"""
Stage I: Generative Performance & Set Properties Evaluation
Extracts macro-metrics (Validity, Uniqueness, Novelty, and Internal Diversity)
from the combined method outputs for each (method, test_dataset) group.

Corresponds to Chunks 2-5 of metrics_calculation/7.1. results_5.ipynb.
"""

import os
import sys
import argparse
import urllib.request
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


def load_external_reference_set(cache_path="TheraFlow/method_processing/data/moses_test.csv"):
    """
    Downloads and caches the standardized external reference set (MOSES test set, ~176k molecules)
    used for drug-like novelty checks.
    """
    cache_file = Path(cache_path)
    if not cache_file.exists():
        url = "https://github.com/molecularsets/moses/raw/master/data/test.csv"
        print(f"Downloading MOSES external reference set from {url}...")
        os.makedirs(cache_file.parent, exist_ok=True)
        try:
            urllib.request.urlretrieve(url, cache_file)
            print(f"Downloaded and cached reference set to {cache_file}")
        except Exception as e:
            print(f"Warning: Failed to download reference set: {e}. Novelty will default to 0.")
            return set()
            
    print(f"Canonicalizing external reference set from {cache_file}...")
    ref_df = pd.read_csv(cache_file)
    external_ref_set = set()
    for s in tqdm(ref_df["SMILES"], desc="Loading Reference Molecules"):
        if isinstance(s, str) and s.strip():
            mol = Chem.MolFromSmiles(s)
            if mol:
                external_ref_set.add(Chem.MolToSmiles(mol, canonical=True))
                
    print(f"External Reference Set Loaded: {len(external_ref_set)} unique canonical molecules.")
    return external_ref_set


def calculate_macro_metrics(smiles_list, training_smiles_set=None):
    """
    Calculates set-based metrics for a list of SMILES:
    1. Validity: Fraction of valid chemical structures.
    2. Uniqueness: Fraction of unique valid structures.
    3. Novelty: Fraction of unique valid molecules not found in the external reference set.
    4. Internal Diversity (IntDiv): Average pairwise Tanimoto distance of Morgan fingerprints.
    """
    total = len(smiles_list)
    if total == 0:
        return {"Validity": 0.0, "Uniqueness": 0.0, "Novelty": 0.0, "IntDiv": 0.0}

    # 1. Validity
    mols = [Chem.MolFromSmiles(s) for s in smiles_list if isinstance(s, str) and s.strip()]
    valid_mols = [m for m in mols if m is not None]
    valid_smiles = [Chem.MolToSmiles(m, canonical=True) for m in valid_mols]
    valid_count = len(valid_mols)
    
    # 2. Uniqueness
    unique_smiles = set(valid_smiles)
    unique_count = len(unique_smiles)
    
    # 3. Novelty
    novel_count = 0
    if training_smiles_set is not None and unique_count > 0:
        novel_smiles = [s for s in unique_smiles if s not in training_smiles_set]
        novel_count = len(novel_smiles)
        
    # 4. Internal Diversity
    div_score = 0.0
    if unique_count > 1:
        try:
            mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
            fps = [mfpgen.GetFingerprint(Chem.MolFromSmiles(s)) for s in unique_smiles]
            
            # Subsample if group is large to avoid O(N^2) memory bottleneck
            if len(fps) > 500:
                idx = np.random.choice(len(fps), 500, replace=False)
                fps = [fps[i] for i in idx]
                
            distances = []
            for i in range(len(fps)):
                sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
                distances.extend([1.0 - x for x in sims])
            if len(distances) > 0:
                div_score = float(np.mean(distances))
        except Exception as e:
            print(f"Warning calculating Internal Diversity: {e}")

    return {
        "Validity": round(valid_count / total, 4),
        "Uniqueness": round(unique_count / valid_count, 4) if valid_count > 0 else 0.0,
        "Novelty": round(novel_count / unique_count, 4) if unique_count > 0 else 0.0,
        "IntDiv": round(div_score, 4)
    }


def main():
    parser = argparse.ArgumentParser(description="Calculate generative performance set properties across methods and test datasets.")
    parser.add_argument(
        "--input", 
        type=str, 
        default="/users/stgrad/rojashrj/Projects/omics2mol/TheraFlow/method_processing/data/processed_data.csv",
        help="Path to processed_data.csv containing SMILES, method, and test_dataset."
    )
    parser.add_argument(
        "--output", 
        type=str, 
        default="TheraFlow/outputs/set_properties.csv",
        help="Path for output set_properties.csv file."
    )
    parser.add_argument(
        "--ref_cache",
        type=str,
        default="TheraFlow/method_processing/data/moses_test.csv",
        help="Path to local cache for MOSES test dataset."
    )
    args = parser.parse_args()
    
    input_path = Path(args.input)
    if not input_path.exists():
        # Check relative fallback
        rel_path = Path("TheraFlow/method_processing/data/processed_data.csv")
        if rel_path.exists():
            input_path = rel_path
        else:
            raise FileNotFoundError(f"Input file not found: {args.input}")
            
    print(f"Reading input dataset from: {input_path}")
    df = pd.read_csv(input_path)
    print(f"Total rows loaded: {len(df)}")
    
    # Verify expected columns
    if "SMILES" not in df.columns or "method" not in df.columns or "test_dataset" not in df.columns:
        # Check case-insensitive
        df.columns = [c.upper() if c.lower() == "smiles" else c for c in df.columns]
        
    external_ref = load_external_reference_set(args.ref_cache)
    
    macro_results = []
    grouped = df.groupby(["method", "test_dataset"])
    print(f"Computing set properties across {len(grouped)} (method, test_dataset) groups...")
    
    for (method_name, dataset_name), group_df in tqdm(grouped, desc="Evaluating Groups"):
        smiles_list = group_df["SMILES"].dropna().tolist()
        metrics = calculate_macro_metrics(smiles_list, training_smiles_set=external_ref)
        metrics["method"] = method_name
        metrics["test_dataset"] = dataset_name
        macro_results.append(metrics)
        
    summary_df = pd.DataFrame(macro_results)
    cols = ["method", "test_dataset", "Validity", "Uniqueness", "Novelty", "IntDiv"]
    summary_df = summary_df[cols]
    
    output_path = Path(args.output)
    os.makedirs(output_path.parent if str(output_path.parent) != "" else ".", exist_ok=True)
    summary_df.to_csv(output_path, index=False)
    
    print("\n=== Set Properties Summary ===")
    print(summary_df.to_string(index=False))
    print(f"\nSuccessfully wrote set properties to: {output_path.resolve()}")


if __name__ == "__main__":
    main()
