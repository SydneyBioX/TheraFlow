import os
import sys
import argparse
import pandas as pd
import torch

# Paths setup
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
WORKSPACE_ROOT = os.path.dirname(os.path.dirname(BASE_DIR))
DATA_DIR = os.path.join(BASE_DIR, "data")
REST_METHODS_DIR = os.path.join(DATA_DIR, "test_datasets/rest_methods")
OUTPUTS_DIR = os.path.join(BASE_DIR, "outputs/bicev")
BICEV_ROOT = os.path.join(WORKSPACE_ROOT, "BiCEV")

# Import setup
sys.path.append(BICEV_ROOT)

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
    parser = argparse.ArgumentParser(description="Orchestrator to run BiCEV on test datasets")
    parser.add_argument("--dataset", choices=DATASETS, help="Run on a specific dataset")
    parser.add_argument("--run-all", action="store_true", help="Run all datasets sequentially")
    
    args = parser.parse_args()
    if not (args.dataset or args.run_all):
        parser.print_help()
        sys.exit(1)
        
    datasets = DATASETS if args.run_all else [args.dataset]
    
    # Change working directory to BiCEV directory context for weight loading
    original_cwd = os.getcwd()
    os.chdir(BICEV_ROOT)
    
    from model.cea_model import BiCEV
    
    print("Loading BiCEV model...")
    model = BiCEV(
        cla_enc_weight='weight/cla_encoder_weight.ckpt', 
        cla_dec_weight='weight/cla_decoder_weight.ckpt'
    )
    model.eval()
    
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device)
    print("BiCEV model loaded successfully!")
    
    # Restore CWD
    os.chdir(original_cwd)
    
    # Create outputs directory
    os.makedirs(OUTPUTS_DIR, exist_ok=True)
    
    for d in datasets:
        src_path = os.path.join(REST_METHODS_DIR, f"{d}.csv")
        if not os.path.exists(src_path):
            print(f"Error: Preprocessed file {src_path} does not exist. Run preprocessing first.")
            continue
            
        print(f"\nProcessing {d} using BiCEV model...")
        df_sig = pd.read_csv(src_path)
        
        # Grab values (excluding first ID column)
        expression_values = df_sig.iloc[:, 1:].values[0]
        input_tensor = torch.tensor(expression_values, dtype=torch.float32).to(device)
        
        print("Generating molecules...")
        with torch.no_grad():
            generated_output = model.predict_from_single_expression(input_tensor, amount=50)
            
        df_results = pd.DataFrame(generated_output[0], columns=['SMILES'])
        dest_path = os.path.join(OUTPUTS_DIR, f"{d}.csv")
        df_results.to_csv(dest_path, index=False)
        print(f"Saved {len(df_results)} SMILES to {dest_path}")

if __name__ == "__main__":
    main()
