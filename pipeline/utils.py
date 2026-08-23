import os
import pandas as pd
from pathlib import Path
from tqdm import tqdm
from joblib import Parallel, delayed
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, AllChem

# Consolidated lookup mapping test datasets to their reference inhibitor SMILES
INHIBITOR_MAPPING = {
    # AKT1 inhibitors
    "lincs1000_1": "C1CN(CCC1N2C3=CC=CC=C3NC2=O)CC4=CC=C(C=C4)C5=C(N=C6C=C7C(=NC=N7)C=C6N5)C8=CC=CC=C8",
    "lincs1000_1_m": "C1CN(CCC1N2C3=CC=CC=C3NC2=O)CC4=CC=C(C=C4)C5=C(N=C6C=C7C(=NC=N7)C=C6N5)C8=CC=CC=C8",
    "smilesgen_1": "C1CN(CCC1N2C3=CC=CC=C3NC2=O)CC4=CC=C(C=C4)C5=C(N=C6C=C7C(=NC=N7)C=C6N5)C8=CC=CC=C8",
    "smilesgen_1_m": "C1CN(CCC1N2C3=CC=CC=C3NC2=O)CC4=CC=C(C=C4)C5=C(N=C6C=C7C(=NC=N7)C=C6N5)C8=CC=CC=C8",
    "excape_1": "C1CN(CCC1N2C3=CC=CC=C3NC2=O)CC4=CC=C(C=C4)C5=C(N=C6C=C7C(=NC=N7)C=C6N5)C8=CC=CC=C8",
    "gxrnn_1": "C1CN(CCC1N2C3=CC=CC=C3NC2=O)CC4=CC=C(C=C4)C5=C(N=C6C=C7C(=NC=N7)C=C6N5)C8=CC=CC=C8",
    # Venlafaxine (depression)
    "lincs1000_2": "CN(C)CC(C1=CC=C(C=C1)OC)C2(CCCCC2)O",
    "lincs1000_2_m": "CN(C)CC(C1=CC=C(C=C1)OC)C2(CCCCC2)O",
    # Latanoprost (Glaucoma)
    "glaucoma": "CC(C)OC(=O)CCC/C=C\\C[C@H]1[C@H](C[C@H]([C@@H]1CC[C@H](CCC2=CC=CC=C2)O)O)O",
    "glaucoma_bulk": "CC(C)OC(=O)CCC/C=C\\C[C@H]1[C@H](C[C@H]([C@@H]1CC[C@H](CCC2=CC=CC=C2)O)O)O",
    # Olaparib (Cancer)
    "breast_cancer": "C1CC1C(=O)N2CCN(CC2)C(=O)C3=C(C=CC(=C3)CC4=NNC(=O)C5=CC=CC=C54)F",
    "prostate_cancer": "C1CC1C(=O)N2CCN(CC2)C(=O)C3=C(C=CC(=C3)CC4=NNC(=O)C5=CC=CC=C54)F",
    # Ethanol
    "triomphe_1": "CCO",
    # Carvedilol (Cardiomyopathy)
    "cardiomyopathy": "COc1ccccc1OCCNCC(O)COc1cccc2[nH]c3ccccc3c12",
    # Selegiline (Parkinson's)
    "parkinsons": "C#CCN(C)[C@H](C)Cc1ccccc1",
    # Sorafenib (HCC)
    "hcc": "CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(C(F)(F)F)c3)cc2)ccn1"
}

# Mapping methods to their relative subfolder names
METHOD_FOLDERS = {
    "smilesgen": "smilesgen",
    "gexmolgen": "gexmolgen",
    "transgem": "transgem",
    "bicev": "Bicev",
    "triomphe": "triomphe",
    "paccmann_rl": "paccmannrl",
    "biaee": "biaee",
    "fame": "fame",
    "gex2sgen": "gex2sgen",
    "gxrnn": "gxrnn"
}


def process_method_file(file_path, method):
    """Processes a single method output file and standardizes column names and metadata."""
    file_path = Path(file_path)
    test_dataset = file_path.stem
    
    # Read file based on type
    if method == "triomphe":
        df = pd.read_csv(file_path, sep="\t")
    else:
        df = pd.read_csv(file_path)
        
    if len(df) == 0:
        return pd.DataFrame(columns=["SMILES", "method", "test_dataset", "inhibitor"])

    # Standardize column selection and rename to 'SMILES'
    if method == "gexmolgen":
        # Flatten from index 1 onwards
        df = pd.DataFrame(df.iloc[:, 1:].values.flatten(), columns=["SMILES"])
    elif method in ["transgem", "triomphe"]:
        df = df.iloc[:, [1]]
        df.columns = ["SMILES"]
    elif method == "paccmann_rl":
        df = df[["SMILES"]]
    elif method == "gex2sgen":
        df = df[["smiles"]]
        df.columns = ["SMILES"]
    else:
        # Check standard columns
        if "SMILES" in df.columns:
            df = df[["SMILES"]]
        elif "smiles" in df.columns:
            df = df[["smiles"]]
            df.columns = ["SMILES"]
        else:
            # Fallback: take first column
            df = df.iloc[:, [0]]
            df.columns = ["SMILES"]

    # Clean up SMILES strings
    df = df.dropna(subset=["SMILES"])
    df["SMILES"] = df["SMILES"].astype(str).str.strip()
    df = df[df["SMILES"] != ""]
    
    # Add metadata
    df["method"] = method
    df["test_dataset"] = test_dataset
    
    # Map inhibitor
    inhibitor = INHIBITOR_MAPPING.get(test_dataset, "Unknown")
    df["inhibitor"] = inhibitor
    
    return df


def merge_outputs(method_outputs_dir, output_csv_path):
    """Merges output files from all specified method folders into a single CSV file."""
    base_path = Path(method_outputs_dir)
    merged_dfs = []
    
    for method, subfolder in METHOD_FOLDERS.items():
        subfolder_path = base_path / subfolder
        if not subfolder_path.exists():
            print(f"Warning: Directory {subfolder_path} does not exist. Skipping method '{method}'.")
            continue
            
        print(f"Processing method directory: {subfolder_path}")
        # Match csv or txt files
        files = list(subfolder_path.glob("*.csv")) + list(subfolder_path.glob("*.txt"))
        for file in files:
            try:
                df = process_method_file(file, method)
                if len(df) > 0:
                    merged_dfs.append(df)
            except Exception as e:
                print(f"Error processing file {file} for method '{method}': {e}")
                
    if not merged_dfs:
        raise ValueError(f"No valid output files found to merge in {method_outputs_dir}")
        
    combined = pd.concat(merged_dfs, ignore_index=True)
    
    # Write to target CSV
    os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
    combined.to_csv(output_csv_path, index=False)
    print(f"Successfully combined {len(combined)} rows. Output saved to {output_csv_path}")
    return combined


def process_molecule_3d(row_dict):
    """Helper to perform 2D filtering and RDKit 3D conformation embedding on a single molecule dict."""
    smiles = row_dict.get("SMILES")
    if not isinstance(smiles, str) or smiles.strip() == "":
        return None
        
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
        
    try:
        # 1. Apply 2D Lipinski Filters
        if not (Descriptors.MolWt(mol) <= 500 and 
                Crippen.MolLogP(mol) <= 5 and 
                Lipinski.NumHDonors(mol) <= 5 and 
                Lipinski.NumHAcceptors(mol) <= 10):
            return None
            
        # 2. Apply Complexity Filters
        if not (Lipinski.NumRotatableBonds(mol) <= 18 and 
                mol.GetNumHeavyAtoms() <= 60):
            return None
            
        # 3. Add Hs and generate 3D embedding using ETKDGv3
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        if AllChem.EmbedMolecule(mol, params) == -1:
            return None
            
        # Attach all original metadata fields as SD properties
        for key, val in row_dict.items():
            if pd.notna(val):
                mol.SetProp(str(key), str(val))
                
        return mol
    except:
        return None


def filter_for_posebusters(properties_csv, output_sdf_path, n_jobs=4):
    """Filters molecules from a CSV properties file and writes their 3D structures to an SDF file."""
    df = pd.read_csv(properties_csv)
    
    # Deduplicate based on SMILES column
    initial_count = len(df)
    df = df.drop_duplicates(subset=["SMILES"])
    print(f"Deduplicated SMILES: {initial_count} rows -> {len(df)} unique rows.")
    
    # Prepare records for parallel processing
    tasks = df.to_dict("records")
    print(f"Generating 3D conformations using {n_jobs} parallel jobs...")
    
    results = Parallel(n_jobs=n_jobs)(
        delayed(process_molecule_3d)(row) for row in tqdm(tasks)
    )
    
    # Collect valid 3D structures and attach properties in the main thread
    valid_mols = []
    for row, mol in zip(tasks, results):
        if mol is not None:
            for key, val in row.items():
                if pd.notna(val):
                    mol.SetProp(str(key), str(val))
            valid_mols.append(mol)
            
    # Write to SDF
    os.makedirs(os.path.dirname(output_sdf_path), exist_ok=True)
    writer = Chem.SDWriter(output_sdf_path)
    for mol in valid_mols:
        writer.write(mol)
    writer.close()
    
    print(f"Success! Saved {len(valid_mols)} unique 3D structures to {output_sdf_path}")
    return len(valid_mols)
