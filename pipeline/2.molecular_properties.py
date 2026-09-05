#!/usr/bin/env python3
"""
Stage II: Molecular Properties & Reference Similarity Metrics
Calculates 2D physicochemical descriptors (MW, LogP, TPSA, QED, SA_Score, F_Csp3,
HBA, HBD, Rotatable Bonds, Rings, Aromatic Rings, Stereo Centers, PAINS/BRENK/NIH Safety Alerts)
and similarity suite against the dataset-specific reference inhibitor.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm
from rdkit import Chem, DataStructs
from rdkit.Chem import (
    AllChem, Descriptors, rdMolDescriptors,
    FilterCatalog, QED, RDConfig, MACCSkeys, Crippen, Lipinski
)
from rdkit.Chem.Scaffolds import MurckoScaffold

# SA Score setup
try:
    sa_path = os.path.join(RDConfig.RDContribDir, "SA_Score")
    if sa_path not in sys.path:
        sys.path.append(sa_path)
    import sascorer
except Exception as e:
    print(f"Warning: Could not configure sascorer path: {e}")
    sascorer = None

# Reference inhibitor SMILES dictionary for all test datasets
INHIBITOR_LOOKUP = {
    "glaucoma": "CC(C)OC(=O)CCC/C=C\C[C@H]1[C@H](C[C@H]([C@@H]1CC[C@H](CCC2=CC=CC=C2)O)O)O",
    "glaucoma_bulk": "CC(C)OC(=O)CCC/C=C\C[C@H]1[C@H](C[C@H]([C@@H]1CC[C@H](CCC2=CC=CC=C2)O)O)O",
    "cardiomyopathy": "COc1ccccc1OCCNCC(O)COc1cccc2[nH]c3ccccc3c12",
    "parkinsons": "C#CCN(C)[C@H](C)Cc1ccccc1",
    "hcc": "CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(C(F)(F)F)c3)cc2)ccn1",
    "lincs1000_1": "C1CN(CCC1N2C3=CC=CC=C3NC2=O)CC4=CC=C(C=C4)C5=C(N=C6C=C7C(=NC=N7)C=C6N5)C8=CC=CC=C8",
    "lincs1000_2": "CN(C)CC(C1=CC=C(C=C1)OC)C2(CCCCC2)O"
}


def get_safety_catalog():
    """Initializes PAINS, BRENK, and NIH filter catalog."""
    params = FilterCatalog.FilterCatalogParams()
    for cat in [
        FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS, 
        FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK, 
        FilterCatalog.FilterCatalogParams.FilterCatalogs.NIH
    ]:
        params.AddCatalog(cat)
    return FilterCatalog.FilterCatalog(params)


def prepare_inhibitor_data(smiles):
    """Pre-calculates Morgan, MACCS, and Murcko Scaffold reference fingerprints."""
    if not isinstance(smiles, str) or not smiles.strip() or smiles == "Unknown":
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None
        try:
            scaf = MurckoScaffold.GetScaffoldForMol(mol)
            scaf_fp = AllChem.GetMorganFingerprintAsBitVect(scaf, 2) if scaf else None
        except:
            scaf_fp = None
            
        return {
            "morg": AllChem.GetMorganFingerprintAsBitVect(mol, 2),
            "maccs": MACCSkeys.GenMACCSKeys(mol),
            "scaf": scaf_fp
        }
    except:
        return None


def get_sim_suite(ref_fp, gen_fp, label):
    """Computes Tanimoto, Dice, Tversky, and Cosine similarity."""
    if gen_fp is None or ref_fp is None:
        return {
            f"Sim_{label}_Tani": 0.0,
            f"Sim_{label}_Dice": 0.0,
            f"Sim_{label}_Tversky": 0.0,
            f"Sim_{label}_Cosine": 0.0
        }
    return {
        f"Sim_{label}_Tani": round(DataStructs.TanimotoSimilarity(ref_fp, gen_fp), 3),
        f"Sim_{label}_Dice": round(DataStructs.DiceSimilarity(ref_fp, gen_fp), 3),
        f"Sim_{label}_Tversky": round(DataStructs.TverskySimilarity(ref_fp, gen_fp, 0.1, 0.9), 3),
        f"Sim_{label}_Cosine": round(DataStructs.CosineSimilarity(ref_fp, gen_fp), 3)
    }


def calculate_row_properties(row, catalog, ref_lookups):
    """Calculates all 2D properties and reference similarities for a single row."""
    smiles = row.get("SMILES")
    dataset_name = row.get("test_dataset")
    
    if not isinstance(smiles, str) or not smiles.strip():
        return None
        
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
        
    # Molecule fingerprints
    fp_morg = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
    fp_maccs = MACCSkeys.GenMACCSKeys(mol)
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        fp_scaf = AllChem.GetMorganFingerprintAsBitVect(scaffold, 2) if scaffold else None
    except:
        fp_scaf = None
        
    # Reference lookup
    ref = ref_lookups.get(dataset_name)
    morg_sims = get_sim_suite(ref.get("morg") if ref else None, fp_morg, "Morg")
    maccs_sims = get_sim_suite(ref.get("maccs") if ref else None, fp_maccs, "MACCS")
    scaf_sims = get_sim_suite(ref.get("scaf") if ref else None, fp_scaf, "Scaffold")
    
    sa_score = round(sascorer.calculateScore(mol), 3) if sascorer else np.nan
    
    metrics = {
        "QED": round(QED.qed(mol), 3),
        "SA_Score": sa_score,
        "MW": round(Descriptors.MolWt(mol), 2),
        "LogP": round(Descriptors.MolLogP(mol), 2),
        "TPSA": round(Descriptors.TPSA(mol), 2),
        "F_Csp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 2),
        "Passes_Safety": not catalog.HasMatch(mol),
        "HBA": Descriptors.NumHAcceptors(mol),
        "HBD": Descriptors.NumHDonors(mol),
        "Rotatable_Bonds": Descriptors.NumRotatableBonds(mol),
        "Ring_Count": rdMolDescriptors.CalcNumRings(mol),
        "Aromatic_Rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "Stereo_Centers": len(Chem.FindMolChiralCenters(mol, includeUnassigned=True)),
    }
    metrics.update(morg_sims)
    metrics.update(maccs_sims)
    metrics.update(scaf_sims)
    return metrics


def apply_benchmark_filters(df):
    """Filters molecules based on safety alerts, SA score <= 5, MW <= 600, -1 <= LogP <= 6, and Lipinski <= 1."""
    filtered_df = df[
        (df["Passes_Safety"] == True) &
        (df["SA_Score"] <= 5.0) &
        (df["MW"] <= 600) &
        (df["LogP"] >= -1.0) & (df["LogP"] <= 6.0)
    ].copy()
    
    violations = (
        (filtered_df["MW"] > 500).astype(int) +
        (filtered_df["LogP"] > 5.0).astype(int) +
        (filtered_df["HBD"] > 5).astype(int) +
        (filtered_df["HBA"] > 10).astype(int)
    )
    filtered_df["Lipinski_Violations"] = violations
    return filtered_df[filtered_df["Lipinski_Violations"] <= 1].copy()


def main():
    parser = argparse.ArgumentParser(description="Calculate 2D molecular properties and reference similarity metrics.")
    parser.add_argument(
        "--input", 
        type=str, 
        default="/users/stgrad/rojashrj/Projects/omics2mol/TheraFlow/method_processing/data/processed_data.csv",
        help="Input CSV containing SMILES, method, and test_dataset."
    )
    parser.add_argument(
        "--output", 
        type=str, 
        default="TheraFlow/outputs/molecular_properties_metrics.csv",
        help="Output CSV for full molecular properties."
    )
    parser.add_argument(
        "--output_filtered", 
        type=str, 
        default="TheraFlow/outputs/molecular_properties_filtered.csv",
        help="Output CSV for benchmark-filtered molecules."
    )
    args = parser.parse_args()
    
    input_path = Path(args.input)
    if not input_path.exists():
        rel_path = Path("TheraFlow/method_processing/data/processed_data.csv")
        if rel_path.exists():
            input_path = rel_path
        else:
            raise FileNotFoundError(f"Input file not found: {args.input}")
            
    print(f"Reading input SMILES from {input_path}...")
    df = pd.read_csv(input_path)
    print(f"Loaded {len(df)} candidate molecules.")
    
    # Pre-calculate reference fingerprints
    print("Pre-calculating reference inhibitor fingerprints for benchmark datasets...")
    ref_lookups = {}
    for ds, smi in INHIBITOR_LOOKUP.items():
        ref_data = prepare_inhibitor_data(smi)
        if ref_data:
            ref_lookups[ds] = ref_data
            
    catalog = get_safety_catalog()
    
    print("Computing 2D physicochemical properties, structural alerts, and reference similarities...")
    metrics_list = []
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Calculating Properties"):
        metrics_list.append(calculate_row_properties(row, catalog, ref_lookups))
        
    valid_mask = [m is not None for m in metrics_list]
    df_valid = df[valid_mask].reset_index(drop=True)
    metrics_df = pd.DataFrame([m for m in metrics_list if m is not None])
    
    combined = pd.concat([df_valid, metrics_df], axis=1)
    
    out_file = Path(args.output)
    os.makedirs(out_file.parent if str(out_file.parent) != "" else ".", exist_ok=True)
    combined.to_csv(out_file, index=False)
    print(f"Successfully saved {len(combined)} properties to {out_file.resolve()}")
    
    if args.output_filtered:
        filtered_df = apply_benchmark_filters(combined)
        out_filt = Path(args.output_filtered)
        os.makedirs(out_filt.parent if str(out_filt.parent) != "" else ".", exist_ok=True)
        filtered_df.to_csv(out_filt, index=False)
        print(f"Successfully saved {len(filtered_df)} / {len(combined)} benchmark-filtered properties to {out_filt.resolve()}")


if __name__ == "__main__":
    main()
