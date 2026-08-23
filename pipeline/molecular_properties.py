import os
import sys
import urllib.request
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm
from rdkit import Chem, DataStructs
from rdkit.Chem import (
    AllChem, Descriptors, rdFingerprintGenerator, rdMolDescriptors,
    FilterCatalog, QED, RDConfig, MACCSkeys
)
from rdkit.Chem.Scaffolds import MurckoScaffold

# Setup RDKit SA Score path and import sascorer
try:
    sa_path = os.path.join(RDConfig.RDContribDir, "SA_Score")
    if sa_path not in sys.path:
        sys.path.append(sa_path)
    import sascorer
except Exception as e:
    print(f"Warning: Could not configure sascorer path: {e}")
    sascorer = None


def get_safety_catalog():
    """Returns the combined PAINS, BRENK, and NIH safety filter catalog."""
    params = FilterCatalog.FilterCatalogParams()
    for cat in [
        FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS, 
        FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK, 
        FilterCatalog.FilterCatalogParams.FilterCatalogs.NIH
    ]:
        params.AddCatalog(cat)
    return FilterCatalog.FilterCatalog(params)


def prepare_inhibitor_data(smiles):
    """Helper to pre-calculate all reference fingerprints for an inhibitor."""
    if not isinstance(smiles, str) or smiles.strip() == "" or smiles == "Unknown":
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None
        # Scaffold calculation
        try:
            scaf = MurckoScaffold.GetScaffoldForMol(mol)
            scaf_fp = AllChem.GetMorganFingerprintAsBitVect(scaf, 2)
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
    """Calculates Tanimoto, Dice, Tversky, and Cosine similarities between two fingerprints."""
    if gen_fp is None or ref_fp is None:
        return {
            f"Sim_{label}_Tani": 0.0,
            f"Sim_{label}_Dice": 0.0,
            f"Sim_{label}_Tversky": 0.0,
            f"Sim_{label}_Cosine": 0.0
        }
    try:
        return {
            f"Sim_{label}_Tani": round(DataStructs.TanimotoSimilarity(ref_fp, gen_fp), 3),
            f"Sim_{label}_Dice": round(DataStructs.DiceSimilarity(ref_fp, gen_fp), 3),
            f"Sim_{label}_Tversky": round(DataStructs.TverskySimilarity(ref_fp, gen_fp, 0.1, 0.9), 3),
            f"Sim_{label}_Cosine": round(DataStructs.CosineSimilarity(ref_fp, gen_fp), 3)
        }
    except:
        return {
            f"Sim_{label}_Tani": 0.0,
            f"Sim_{label}_Dice": 0.0,
            f"Sim_{label}_Tversky": 0.0,
            f"Sim_{label}_Cosine": 0.0
        }


def get_molecular_metrics_for_row(row, universal_catalog, inhibitor_lookups):
    """Calculates all 2D descriptors and similarity suites for a single dataframe row."""
    smiles = row.get("SMILES")
    dataset_name = row.get("test_dataset")
    inhibitor_smiles = row.get("inhibitor")
    
    if pd.isna(smiles) or not isinstance(smiles, str) or smiles.strip() == "":
        return None
        
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
        
    # Get reference inhibitor info (fallback to calculating on the fly if not in lookups)
    ref = inhibitor_lookups.get(dataset_name)
    if not ref and isinstance(inhibitor_smiles, str) and inhibitor_smiles != "Unknown":
        ref = prepare_inhibitor_data(inhibitor_smiles)
        
    # Calculate descriptors
    fp_morg = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
    fp_maccs = MACCSkeys.GenMACCSKeys(mol)
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        fp_scaf = AllChem.GetMorganFingerprintAsBitVect(scaffold, 2)
    except:
        fp_scaf = None
        
    morg_sims = get_sim_suite(ref["morg"] if ref else None, fp_morg, "Morg")
    maccs_sims = get_sim_suite(ref["maccs"] if ref else None, fp_maccs, "MACCS")
    scaf_sims = get_sim_suite(ref["scaf"] if ref else None, fp_scaf, "Scaffold")
    
    # SA Score calculation
    sa_score = 0.0
    if sascorer:
        try:
            sa_score = sascorer.calculateScore(mol)
        except:
            pass
            
    metrics = {
        "QED": round(QED.qed(mol), 3),
        "SA_Score": round(sa_score, 3),
        "MW": round(Descriptors.MolWt(mol), 2),
        "LogP": round(Descriptors.MolLogP(mol), 2),
        "TPSA": round(Descriptors.TPSA(mol), 2),
        "F_Csp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 2),
        "Passes_Safety": not universal_catalog.HasMatch(mol),
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


def rank_all_properties(df):
    """Calculates rank columns and a global consensus score/rank grouped by method and test_dataset."""
    group_cols = ["method", "test_dataset"]
    ranked_df = df.copy()
    
    higher_is_better = [
        "Sim_Morg_Tani", "Sim_Morg_Dice", "Sim_Morg_Tversky",
        "Sim_Morg_Cosine", "Sim_MACCS_Tani", "Sim_MACCS_Dice",
        "Sim_MACCS_Tversky", "Sim_MACCS_Cosine", "Sim_Scaffold_Tani",
        "Sim_Scaffold_Dice", "Sim_Scaffold_Tversky", "Sim_Scaffold_Cosine",
        "QED", "F_Csp3"
    ]
    lower_is_better = ["SA_Score", "MW", "LogP", "TPSA", "HBD", "HBA"]
    
    rank_cols = []
    for col in higher_is_better:
        if col in ranked_df.columns:
            col_rank = f"{col}_Rank"
            ranked_df[col_rank] = ranked_df.groupby(group_cols)[col].rank(ascending=False, method="min")
            rank_cols.append(col_rank)
            
    for col in lower_is_better:
        if col in ranked_df.columns:
            col_rank = f"{col}_Rank"
            ranked_df[col_rank] = ranked_df.groupby(group_cols)[col].rank(ascending=True, method="min")
            rank_cols.append(col_rank)
            
    if rank_cols:
        ranked_df["Global_Consensus_Score"] = ranked_df[rank_cols].mean(axis=1)
        ranked_df["Global_Consensus_Rank"] = ranked_df.groupby(group_cols)["Global_Consensus_Score"].rank(ascending=True, method="min")
        
    return ranked_df


def load_external_reference_set(local_cache_path="data/moses_test.csv"):
    """Loads external reference set (MOSES test set) from local cache or downloads it if missing."""
    external_ref_set = set()
    cache_file = Path(local_cache_path)
    
    if not cache_file.exists():
        url = "https://github.com/molecularsets/moses/raw/master/data/test.csv"
        print(f"MOSES reference set not found locally. Downloading from {url}...")
        try:
            os.makedirs(cache_file.parent, exist_ok=True)
            urllib.request.urlretrieve(url, cache_file)
            print("Download complete.")
        except Exception as e:
            print(f"Error downloading reference set: {e}. Novelty calculations will be skipped.")
            return external_ref_set
            
    try:
        print(f"Loading and canonicalizing external reference set from {cache_file}...")
        ref_df = pd.read_csv(cache_file)
        for s in tqdm(ref_df["SMILES"]):
            mol = Chem.MolFromSmiles(s)
            if mol:
                external_ref_set.add(Chem.MolToSmiles(mol, canonical=True))
        print(f"Loaded {len(external_ref_set)} unique reference molecules.")
    except Exception as e:
        print(f"Error loading reference set: {e}")
        
    return external_ref_set


def calculate_macro_metrics(smiles_list, training_smiles_set=None):
    """Calculates set-based Validity, Uniqueness, Novelty, and Internal Diversity for a list of SMILES."""
    total = len(smiles_list)
    if total == 0:
        return {"Validity": 0.0, "Uniqueness": 0.0, "Novelty": 0.0, "IntDiv": 0.0}

    # 1. Validity
    mols = [Chem.MolFromSmiles(s) for s in smiles_list]
    valid_mols = [m for m in mols if m is not None]
    valid_smiles = [Chem.MolToSmiles(m, canonical=True) for m in valid_mols]
    valid_count = len(valid_mols)
    
    # 2. Uniqueness
    unique_smiles = set(valid_smiles)
    unique_count = len(unique_smiles)
    
    # 3. Novelty (External)
    novel_count = np.nan
    if training_smiles_set is not None and unique_count > 0:
        novel_smiles = [s for s in unique_smiles if s not in training_smiles_set]
        novel_count = len(novel_smiles)

    # 4. Internal Diversity (IntDiv)
    div_score = 0.0
    if unique_count > 1:
        try:
            mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
            fps = [mfpgen.GetFingerprint(Chem.MolFromSmiles(s)) for s in unique_smiles]
            
            # Subsample if set is extremely large to avoid O(N^2) memory/time blowup
            if len(fps) > 500:
                idx = np.random.choice(len(fps), 500, replace=False)
                fps = [fps[i] for i in idx]
                
            distances = []
            for i in range(len(fps)):
                sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
                distances.extend([1.0 - x for x in sims])
            div_score = np.mean(distances) if distances else 0.0
        except Exception as e:
            print(f"Error calculating Internal Diversity: {e}")

    return {
        "Validity": valid_count / total,
        "Uniqueness": unique_count / valid_count if valid_count > 0 else 0.0,
        "Novelty": novel_count / unique_count if unique_count > 0 and not np.isnan(novel_count) else 0.0,
        "IntDiv": div_score
    }


def calculate_properties_pipeline(input_csv_path, output_csv_path, output_summary_csv_path=None, local_ref_path="data/moses_test.csv"):
    """Main pipeline execution for step 2: molecular properties and summary metrics."""
    df = pd.read_csv(input_csv_path)
    print(f"Loaded {len(df)} rows from {input_csv_path} for property calculations.")
    
    # Pre-calculate lookup dictionary for all unique inhibitors present
    unique_inhibitors = df[["test_dataset", "inhibitor"]].drop_duplicates()
    inhibitor_lookups = {}
    for _, row in unique_inhibitors.iterrows():
        name = row["test_dataset"]
        smiles = row["inhibitor"]
        if isinstance(smiles, str) and smiles != "Unknown":
            ref_data = prepare_inhibitor_data(smiles)
            if ref_data:
                inhibitor_lookups[name] = ref_data
                
    catalog = get_safety_catalog()
    
    print("Calculating 2D descriptors and reference similarities...")
    tqdm.pandas()
    metrics = df.progress_apply(
        lambda r: get_molecular_metrics_for_row(r, catalog, inhibitor_lookups), 
        axis=1
    )
    
    # Filter out rows where molecular metrics failed (invalid SMILES)
    valid_mask = metrics.notna()
    df_valid = df[valid_mask].reset_index(drop=True)
    metrics_valid = pd.DataFrame(list(metrics[valid_mask]))
    
    # Combine original data with calculated metrics
    combined = pd.concat([df_valid, metrics_valid], axis=1)
    
    # Rank properties
    print("Calculating global consensus ranks...")
    combined = rank_all_properties(combined)
    
    # Save ranked properties
    os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
    combined.to_csv(output_csv_path, index=False)
    print(f"Saved calculated molecular properties to {output_csv_path}")
    
    # Calculate set-level summary macro-metrics if output path provided
    if output_summary_csv_path:
        print("Calculating set-level generative summary metrics...")
        external_ref = load_external_reference_set(local_ref_path)
        
        summary_results = []
        for (method, test_dataset), group in df.groupby(["method", "test_dataset"]):
            smiles_list = group["SMILES"].tolist()
            macro = calculate_macro_metrics(smiles_list, training_smiles_set=external_ref)
            macro["method"] = method
            macro["test_dataset"] = test_dataset
            summary_results.append(macro)
            
        summary_df = pd.DataFrame(summary_results)
        os.makedirs(os.path.dirname(output_summary_csv_path), exist_ok=True)
        summary_df.to_csv(output_summary_csv_path, index=False)
        print(f"Saved generative macro metrics summary to {output_summary_csv_path}")
        
    return combined
