# paccmann rl has to be cloned from https://github.com/PaccMann/paccmann_rl.git and the environment has to be setup
import os
import sys
import argparse
import pickle
import json
import pandas as pd
import numpy as np
import torch

# Paths setup
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
WORKSPACE_ROOT = os.path.dirname(os.path.dirname(BASE_DIR))
DATA_DIR = os.path.join(BASE_DIR, "data")
PACCMANN_TEST_DIR = os.path.join(DATA_DIR, "test_datasets/paccmann_rl")
OUTPUTS_DIR = os.path.join(BASE_DIR, "outputs/paccmann_rl")
PACCMANN_ROOT = os.path.join(WORKSPACE_ROOT, "paccmann_rl")

# Import setup
sys.path.append(os.path.join(PACCMANN_ROOT, "code"))

DATASETS = [
    "lincs1000_1",
    "lincs1000_2",
    "cardiomyopathy",
    "parkinsons",
    "hcc",
    "glaucoma",
    "glaucoma_bulk"
]

def main():
    parser = argparse.ArgumentParser(description="Orchestrator to run PaccMann RL on test datasets")
    parser.add_argument("--dataset", choices=DATASETS, help="Run on a specific dataset")
    parser.add_argument("--run-all", action="store_true", help="Run all datasets sequentially")
    
    args = parser.parse_args()
    if not (args.dataset or args.run_all):
        parser.print_help()
        sys.exit(1)
        
    datasets = DATASETS if args.run_all else [args.dataset]
    
    # Change working directory to paccmann_rl directory context for relative loading inside packages
    original_cwd = os.getcwd()
    os.chdir(PACCMANN_ROOT)
    
    from paccmann_omics.encoders import ENCODER_FACTORY
    from paccmann_chemistry.models import StackGRUEncoder, StackGRUDecoder, TeacherVAE
    
    # Setup CPU device to bypass library device bugs (as configured in run_paccmann_hcc.py)
    device = torch.device("cpu")
    
    # Model parameters and configuration
    PVAE_PATH   = "models/pvae/"
    BIASED_PATH = "biased_models/breast_paccmann_rl/"
    SVAE_PATH   = "models/svae/"  
    SMILES_LANG = "data/smiles_language_chembl_gdsc_ccle.pkl"
    GENE_LIST   = "data/2128_genes.pkl"
    N_MOLECULES = 100
    TEMPERATURE = 2
    
    # Load model prerequisites
    with open(SMILES_LANG, "rb") as f:
        smiles_language = pickle.load(f)
    with open(GENE_LIST, "rb") as f:
        gene_list = pickle.load(f)
        
    # Load PVAE model
    with open(f"{PVAE_PATH}/model_params.json") as f:
        pvae_params = json.load(f)
    pvae = ENCODER_FACTORY["dense"](pvae_params)
    full_state_dict = torch.load(
        f"{PVAE_PATH}/weights/best_rec_vae.pt",
        map_location=device
    )
    encoder_state_dict = {}
    for key, value in full_state_dict.items():
        if key.startswith("encoder."):
            new_key = key[8:]  
            encoder_state_dict[new_key] = value
    pvae.load_state_dict(encoder_state_dict)
    pvae.eval().to(device)
    
    # Load biased generator
    with open(f"{BIASED_PATH}/model_params.json") as f:
        gen_params = json.load(f)
    with open(f"{SVAE_PATH}/model_params.json") as f:
        svae_params = json.load(f)
        
    svae_encoder = StackGRUEncoder(svae_params)
    svae_decoder = StackGRUDecoder(svae_params)
    generator = TeacherVAE(svae_encoder, svae_decoder)
    generator.load_state_dict(
        torch.load(f"{BIASED_PATH}/weights/gen_30.pt", map_location=device)
    )
    generator.eval().to(device)
    
    # Load normalizers
    with open("pvae_gene_mean.pkl", "rb") as f:
        gene_mean = pickle.load(f)
    with open("pvae_gene_std.pkl", "rb") as f:
        gene_std = pickle.load(f)
        
    # Restore CWD
    os.chdir(original_cwd)
    
    # Create outputs directory
    os.makedirs(OUTPUTS_DIR, exist_ok=True)
    
    for d in datasets:
        src_path = os.path.join(PACCMANN_TEST_DIR, f"{d}.csv")
        if not os.path.exists(src_path):
            print(f"Error: Preprocessed file {src_path} does not exist. Run preprocessing first.")
            continue
            
        print(f"\nProcessing {d} using PaccMann RL...")
        df_sig = pd.read_csv(src_path)
        
        # Align, scale, and normalize (loaded directly from preprocessed paccmann_rl datasets)
        df_aligned = df_sig.reindex(columns=gene_list, fill_value=0.0)
        df_scaled = (df_aligned - gene_mean) / (gene_std + 1e-8)
        df_scaled = df_scaled.fillna(0.0)
        df_scaled = df_scaled.apply(pd.to_numeric, errors='coerce').fillna(0.0)
        df_scaled = df_scaled.clip(-10.0, 10.0)
        
        gex_tensor = torch.tensor(df_scaled.values, dtype=torch.float32).to(device)
        gex_tensor = torch.nan_to_num(gex_tensor, nan=0.0)
        
        # Encode
        with torch.no_grad():
            outputs = pvae(gex_tensor)
            if isinstance(outputs, tuple) and len(outputs) >= 2:
                mu, logvar = outputs[0], outputs[1]
            else:
                mu, logvar = outputs, torch.zeros_like(outputs)
            logvar = torch.clamp(logvar, min=-10.0, max=10.0)
            mu = torch.clamp(mu, min=-10.0, max=10.0)
            std = torch.exp(0.5 * logvar)
            
        # Sample & Generate SMILES
        results = []
        with torch.no_grad():
            for i in range(N_MOLECULES):
                eps = torch.randn_like(std) * TEMPERATURE
                latent = mu + eps * (std + 1e-4)
                latent = torch.nan_to_num(latent, nan=0.0)
                latent = torch.clamp(latent, min=-5.0, max=5.0)
                
                try:
                    gen_output = generator.generate(
                        latent,
                        prime_input=torch.tensor([smiles_language.start_index], device=device),
                        end_token=smiles_language.stop_index,
                        generate_len=gen_params.get("generate_len", 100),
                        temperature=2
                    )
                    token_ids = gen_output[0] if isinstance(gen_output, tuple) else gen_output
                    if isinstance(token_ids, map):
                        token_ids = list(token_ids)
                        
                    token_list = []
                    def flatten_tokens(item):
                        if torch.is_tensor(item):
                            token_list.extend([int(x) for x in item.reshape(-1).tolist()])
                        elif isinstance(item, (list, tuple)):
                            for val in item:
                                flatten_tokens(val)
                        elif isinstance(item, np.ndarray):
                            token_list.extend([int(x) for x in item.flatten().tolist()])
                        else:
                            token_list.append(int(item))
                            
                    flatten_tokens(token_ids)
                    smi = smiles_language.token_indexes_to_smiles(token_list)
                    results.append(smi)
                except RuntimeError:
                    continue
                    
        # Format results (saving all raw generated SMILES)
        df_results = pd.DataFrame(results, columns=['SMILES'])
        dest_path = os.path.join(OUTPUTS_DIR, f"{d}.csv")
        df_results.to_csv(dest_path, index=False)
        print(f"Saved {len(df_results)} SMILES to {dest_path}")

if __name__ == "__main__":
    main()
