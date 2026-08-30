import os
import sys
import argparse
import pandas as pd
import numpy as np
import torch
import shutil

# Paths setup
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
WORKSPACE_ROOT = os.path.dirname(os.path.dirname(BASE_DIR))
DATA_DIR = os.path.join(BASE_DIR, "data")
REST_METHODS_DIR = os.path.join(DATA_DIR, "test_datasets/rest_methods")
OUTPUTS_DIR = os.path.join(BASE_DIR, "outputs/biaee")
BIAEE_ROOT = os.path.join(WORKSPACE_ROOT, "BiAAE")

# Import setup
sys.path.append(BIAEE_ROOT)
sys.path.append(os.path.join(BIAEE_ROOT, "training_scripts"))

from models import BiAAE
from networks import FCDiscriminator, FinetunedEncoder, RNNEncoder

class CondGenExperiment(BiAAE):
    def __init__(self):
        super().__init__(dataset='lincs_rnn')
        rnn = RNNEncoder(out_dim=88)
        self.mine_enc = FinetunedEncoder(rnn, out_dim=self.z_dim)
        self.mine_fc = FCDiscriminator(in_dim=2 * self.z_dim)

DATASETS = [
    "lincs1000_1",
    "lincs1000_2",
    "cardiomyopathy",
    "parkinsons",
    "hcc",
    "glaucoma",
    "glaucoma_bulk"
]

def set_global_seed(seed):
    import random
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)

def main():
    parser = argparse.ArgumentParser(description="Orchestrator to run BiAEE on test datasets")
    parser.add_argument("--dataset", choices=DATASETS, help="Run on a specific dataset")
    parser.add_argument("--run-all", action="store_true", help="Run all datasets sequentially")
    
    args = parser.parse_args()
    if not (args.dataset or args.run_all):
        parser.print_help()
        sys.exit(1)
        
    datasets = DATASETS if args.run_all else [args.dataset]
    
    # 1. Load reference metadata and control NPZ mapping
    gene_meta = pd.read_csv(os.path.join(DATA_DIR, "LINCS/GSE70138/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz"), sep="	")
    npz_mapping = np.load(os.path.join(BIAEE_ROOT, "data/lincs/robust_normalized_controls.npz"), allow_pickle=True)
    master_symbols = npz_mapping['column_names'].flatten()[:978].astype(str)
    entrez_to_symbol = dict(zip(gene_meta['pr_gene_id'].astype(str), gene_meta['pr_gene_symbol']))
    
    # 2. Load model (using the training_scripts context CWD)
    original_cwd = os.getcwd()
    os.chdir(os.path.join(BIAEE_ROOT, "training_scripts"))
    
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = CondGenExperiment()
    checkpoint_path = "../saved_models/lincs_rnn_biaae.ckpt"
    checkpoint = torch.load(checkpoint_path, map_location=device)
    model.load_state_dict(checkpoint)
    model.to(device).eval()
    print("BiAEE model loaded successfully!")
    
    # Restore CWD
    os.chdir(original_cwd)
    
    # 3. Create outputs directory
    os.makedirs(OUTPUTS_DIR, exist_ok=True)
    
    # 4. Run inference for selected datasets
    for d in datasets:
        src_path = os.path.join(REST_METHODS_DIR, f"{d}.csv")
        if not os.path.exists(src_path):
            print(f"Error: Preprocessed file {src_path} does not exist. Run preprocessing first.")
            continue
            
        print(f"\nProcessing {d} using BiAEE model...")
        df_sig = pd.read_csv(src_path)
        
        # Align columns
        if "sig_id" in df_sig.columns:
            expression_data = df_sig.drop(columns=['sig_id'])
        else:
            expression_data = df_sig
            
        expression_data.columns = [entrez_to_symbol.get(col, col) for col in expression_data.columns]
        aligned_data = expression_data.reindex(columns=master_symbols).fillna(0.0)
        
        sample_idx = 0
        gene_tensor = torch.tensor(aligned_data.iloc[sample_idx].values).float().view(1, 978).to(device)
        
        # Add dose column (using 10.0 standard value)
        dose_val = torch.tensor([[10.0]]).to(device)
        model_input = torch.cat([gene_tensor, dose_val], dim=-1)
        
        # Repeat for 128 candidates
        batch_input = model_input.repeat(128, 1)
        
        set_global_seed(42)
        with torch.no_grad():
            generated_smiles = model.sample(batch_input)
            
        df_results = pd.DataFrame(generated_smiles, columns=['SMILES'])
        dest_path = os.path.join(OUTPUTS_DIR, f"{d}.csv")
        df_results.to_csv(dest_path, index=False)
        print(f"Saved {len(df_results)} SMILES to {dest_path}")

if __name__ == "__main__":
    main()
