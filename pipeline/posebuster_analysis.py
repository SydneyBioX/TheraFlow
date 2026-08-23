import os
import pandas as pd
from rdkit import Chem
from posebusters import PoseBusters


def extract_sdf_properties(sdf_path):
    """Reads all metadata tags (e.g. SMILES, method, test_dataset) stored in the SDF molecules."""
    suppl = Chem.SDMolSupplier(sdf_path, sanitize=False)
    props_list = []
    
    for i, mol in enumerate(suppl):
        if mol is None:
            # Add a placeholder if the molecule is unreadable
            props_list.append({"SMILES": "Unreadable", "mol_index": i})
            continue
            
        props = {k: mol.GetProp(k) for k in mol.GetPropNames()}
        props["mol_index"] = i
        
        # If SMILES is missing from tags, try generating it from the structure
        if "SMILES" not in props:
            try:
                mol_no_h = Chem.RemoveHs(mol)
                props["SMILES"] = Chem.MolToSmiles(mol_no_h, canonical=True)
            except:
                props["SMILES"] = "Unknown"
                
        props_list.append(props)
        
    return pd.DataFrame(props_list)


def run_posebusters(sdf_path, output_csv_path):
    """Runs PoseBusters on the SDF file and merges the validation results with the original SDF tags."""
    if not os.path.exists(sdf_path):
        raise FileNotFoundError(f"SDF file not found: {sdf_path}")
        
    print(f"Running PoseBusters check (config='mol') on: {sdf_path}...")
    
    # Initialize PoseBusters with the 'mol' chemical validity config
    buster = PoseBusters(config="mol")
    
    # Run tests on the SDF file with full report enabled (all columns)
    pb_df = buster.bust(sdf_path, full_report=True)
    pb_df_reset = pb_df.reset_index(drop=True)
    
    print("Extracting metadata tags from SDF file...")
    props_df = extract_sdf_properties(sdf_path)
    
    # Make sure we have the same row count
    if len(props_df) != len(pb_df_reset):
        print(f"Warning: SDF metadata row count ({len(props_df)}) does not match PoseBusters row count ({len(pb_df_reset)}). Concat might be misaligned.")
        
    # Concatenate properties and posebusters validation results
    combined = pd.concat([props_df, pb_df_reset], axis=1)
    
    # Save output to CSV
    os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
    combined.to_csv(output_csv_path, index=False)
    print(f"Successfully ran PoseBusters. Results saved to {output_csv_path}")
    
    return combined
