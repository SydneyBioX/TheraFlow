#!/usr/bin/env python3
"""
Step 3: 2D Filtering, InChIKey Deduplication, 3D Embedding, and PoseBusters Validation
Applies InChIKey-based deduplication, Lipinski (MW<=500, LogP<=5, HBD<=5, HBA<=10) &
Complexity (RotatableBonds<=18, HeavyAtoms<=60) filtering, generates 3D conformations via ETKDGv3,
and runs PoseBusters chemical validity tests.

Corresponds to 2.posebusters_run.ipynb and metrics_calculation/4.posebuster_convert_sdf.ipynb.
"""

import os
import argparse
import pandas as pd
from pathlib import Path
from tqdm import tqdm
from joblib import Parallel, delayed
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, AllChem, inchi
from posebusters import PoseBusters


def get_inchikey_and_canonical_smiles(smiles):
    """Generates canonical SMILES and InChIKey for a given SMILES string."""
    if not isinstance(smiles, str) or not smiles.strip():
        return None, None, None
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None, None, None
    try:
        canon_smiles = Chem.MolToSmiles(mol, canonical=True)
        inchikey = inchi.MolToInchiKey(mol)
        return mol, canon_smiles, inchikey
    except Exception:
        return mol, Chem.MolToSmiles(mol, canonical=True), None


def process_single_molecule_3d(row_dict):
    """Applies 2D Lipinski & complexity filters and embeds ETKDGv3 3D coordinates."""
    smiles = row_dict.get("SMILES")
    mol, canon_smiles, inchikey = get_inchikey_and_canonical_smiles(smiles)
    if not mol:
        return None
        
    try:
        # 1. 2D Lipinski filter
        if not (Descriptors.MolWt(mol) <= 500 and 
                Crippen.MolLogP(mol) <= 5 and 
                Lipinski.NumHDonors(mol) <= 5 and 
                Lipinski.NumHAcceptors(mol) <= 10):
            return None
            
        # 2. Complexity filter
        if not (Lipinski.NumRotatableBonds(mol) <= 18 and 
                mol.GetNumHeavyAtoms() <= 60):
            return None
            
        # 3. Add Hydrogens and generate 3D coordinates via ETKDGv3
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        if AllChem.EmbedMolecule(mol, params) == -1:
            return None
            
        # Attach InChIKey and original dataframe properties as SD tags
        if inchikey:
            mol.SetProp("InChIKey", str(inchikey))
        mol.SetProp("SMILES", str(canon_smiles))
        
        for key, val in row_dict.items():
            if pd.notna(val) and key not in ("InChIKey", "SMILES"):
                mol.SetProp(str(key), str(val))
                
        return mol
    except Exception:
        return None


def extract_sdf_tags(sdf_path):
    """Extracts metadata properties from SDF."""
    suppl = Chem.SDMolSupplier(sdf_path, sanitize=False)
    records = []
    for i, mol in enumerate(suppl):
        if mol is None:
            records.append({"SMILES": "Unreadable", "mol_index": i})
            continue
        props = {k: mol.GetProp(k) for k in mol.GetPropNames()}
        props["mol_index"] = i
        if "SMILES" not in props:
            try:
                mol_no_h = Chem.RemoveHs(mol)
                props["SMILES"] = Chem.MolToSmiles(mol_no_h, canonical=True)
            except:
                props["SMILES"] = "Unknown"
        records.append(props)
    return pd.DataFrame(records)


def main():
    parser = argparse.ArgumentParser(description="InChIKey deduplication, 3D SDF conversion, and PoseBusters validation.")
    parser.add_argument(
        "--input", 
        type=str, 
        default="TheraFlow/outputs/molecular_properties_metrics.csv",
        help="Input CSV properties file or processed_data.csv."
    )
    parser.add_argument(
        "--output_sdf", 
        type=str, 
        default="TheraFlow/outputs/filtered_3d.sdf",
        help="Path for generated 3D SDF file."
    )
    parser.add_argument(
        "--output_csv", 
        type=str, 
        default="TheraFlow/outputs/molecular_properties_with_posebusters.csv",
        help="Path for PoseBusters test results CSV."
    )
    parser.add_argument(
        "--n_jobs", 
        type=int, 
        default=4,
        help="Number of parallel workers for 3D embedding."
    )
    args = parser.parse_args()
    
    input_path = Path(args.input)
    if not input_path.exists():
        rel_fallback = Path("TheraFlow/method_processing/data/processed_data.csv")
        if rel_fallback.exists():
            input_path = rel_fallback
        else:
            raise FileNotFoundError(f"Input file not found: {args.input}")
            
    print(f"Reading candidate molecules from {input_path}...")
    df = pd.read_csv(input_path)
    initial_count = len(df)
    
    # -------------------------------------------------------------------------
    # InChIKey & Canonical SMILES Deduplication (as in 2.posebusters_run.ipynb)
    # -------------------------------------------------------------------------
    print("Computing InChIKeys and canonical SMILES for exact chemical deduplication...")
    inchikeys = []
    canon_smiles_list = []
    valid_mask = []
    
    for s in tqdm(df["SMILES"], desc="Generating InChIKeys"):
        mol, c_smi, ikey = get_inchikey_and_canonical_smiles(s)
        if mol is not None:
            valid_mask.append(True)
            inchikeys.append(ikey if ikey else c_smi)
            canon_smiles_list.append(c_smi)
        else:
            valid_mask.append(False)
            inchikeys.append(None)
            canon_smiles_list.append(None)
            
    df = df[valid_mask].copy()
    df["InChIKey"] = [k for k in inchikeys if k is not None]
    df["SMILES"] = [s for s in canon_smiles_list if s is not None]
    
    # Deduplicate primarily on InChIKey (fallback to canonical SMILES)
    df = df.drop_duplicates(subset=["InChIKey"]).copy()
    print(f"Deduplication (via InChIKey): {initial_count} rows -> {len(df)} unique chemical structures.")
    
    # -------------------------------------------------------------------------
    # Parallel 2D Filtering & 3D Conformation Embedding (ETKDGv3)
    # -------------------------------------------------------------------------
    records = df.to_dict("records")
    print(f"Generating 3D conformations (Lipinski + Complexity filters, ETKDGv3) using {args.n_jobs} jobs...")
    
    embedded_mols = Parallel(n_jobs=args.n_jobs)(
        delayed(process_single_molecule_3d)(r) for r in tqdm(records, desc="Embedding 3D")
    )
    
    valid_3d = [m for m in embedded_mols if m is not None]
    print(f"Valid 3D structures generated: {len(valid_3d)} / {len(df)}")
    
    out_sdf = Path(args.output_sdf)
    os.makedirs(out_sdf.parent if str(out_sdf.parent) != "" else ".", exist_ok=True)
    
    writer = Chem.SDWriter(str(out_sdf))
    for mol in valid_3d:
        writer.write(mol)
    writer.close()
    print(f"Saved 3D conformations to {out_sdf.resolve()}")
    
    # -------------------------------------------------------------------------
    # PoseBusters Validation
    # -------------------------------------------------------------------------
    print(f"\nRunning PoseBusters chemical validity check (config='mol') on {out_sdf}...")
    buster = PoseBusters(config="mol")
    pb_df = buster.bust(str(out_sdf), full_report=True)
    pb_df_reset = pb_df.reset_index(drop=True)
    
    props_df = extract_sdf_tags(str(out_sdf))
    combined_pb = pd.concat([props_df, pb_df_reset], axis=1)
    
    out_csv = Path(args.output_csv)
    os.makedirs(out_csv.parent if str(out_csv.parent) != "" else ".", exist_ok=True)
    combined_pb.to_csv(out_csv, index=False)
    print(f"Successfully saved PoseBusters results ({len(combined_pb)} rows) to: {out_csv.resolve()}")


if __name__ == "__main__":
    main()
