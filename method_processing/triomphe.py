import os
import sys
import argparse
import shutil
import numpy as np
import pandas as pd
import torch
from rdkit import Chem
import selfies as sf

# Paths setup
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
WORKSPACE_ROOT = os.path.dirname(os.path.dirname(BASE_DIR))
DATA_DIR = os.path.join(BASE_DIR, "data")
TRIOMPHE_TEST_DIR = os.path.join(DATA_DIR, "test_datasets/triomphe")
OUTPUTS_DIR = os.path.join(BASE_DIR, "outputs/triomphe")
TRIOMPHE_ROOT = os.path.join(WORKSPACE_ROOT, "Triomphe")
CODE_DIR = os.path.join(TRIOMPHE_ROOT, "software/code/structure_generator")

# Import setup
sys.path.append(CODE_DIR)
sys.path.append(os.path.join(TRIOMPHE_ROOT, "software/code/VAE"))

DATASETS = [
    "lincs1000_1",
    "lincs1000_2",
    "cardiomyopathy",
    "parkinsons",
    "hcc",
    "glaucoma",
    "glaucoma_bulk"
]

def get_names(d):
    if d == "lincs1000_1":
        return "Lincs1000_1_cp.csv", "Lincs1000_1_pert.csv", "lincs1000_1.txt", "lincs1000_1_mols.csv"
    elif d == "lincs1000_2":
        return "Lincs1000_2_cp.csv", "Lincs1000_2_pert.csv", "lincs1000_2.txt", "lincs1000_2_mols.csv"
    elif d == "glaucoma_bulk":
        return "glaucoma_bulk_cp.csv", "glaucoma_bulk_pert.csv", "glaucoma_bulk.txt", "glaucoma_bulk_mols.csv"
    else:
        return f"{d}_cp.csv", f"{d}_pert.csv", f"{d}.txt", f"{d}_mols.csv"

def extract_cp_profile(d, sig_info, gene_info, pert_info, cp_example, gctx_path, landmark_genes):
    from cmapPy.pandasGEXpress.parse import parse
    
    print(f"[{d}] Extracting control profiles (X_cp) from raw GSE92742 files...")
    if d == "lincs1000_1":
        cell_line, pert = "MCF7", "AKT-inhibitor-1-2"
        relevant_sigs = sig_info[(sig_info['cell_id'] == cell_line) & (sig_info['pert_iname'] == pert)]
    elif d == "lincs1000_2":
        cell_line, pert = "A375", "venlafaxine"
        relevant_sigs = sig_info[(sig_info['cell_id'] == cell_line) & (sig_info['pert_iname'] == pert)]
    else:
        if d in ["glaucoma", "glaucoma_bulk"]:
            keywords = ["timolol", "brimonidine", "latanoprost"]
        elif d == "hcc":
            keywords = ["sorafenib"]
        elif d == "parkinsons":
            keywords = ["levodopa", "bromocriptine", "pramipexole", "ropinirole", "selegiline", "amantadine"]
        else:  # cardiomyopathy
            keywords = ["carvedilol", "metoprolol", "lisinopril", "enalapril", "valsartan", "losartan", "spironolactone"]
        relevant_sigs = sig_info[sig_info['pert_iname'].str.lower().isin(keywords)]
        
    sig_ids = relevant_sigs['sig_id'].tolist()
    if not sig_ids:
        raise ValueError(f"No signatures found for dataset {d}")
        
    gctoo = parse(gctx_path, rid=landmark_genes, cid=sig_ids)
    expr_df = gctoo.data_df.T 
    expr_df.columns = [f"hsa:{col}" for col in expr_df.columns]
    
    merged_df = relevant_sigs.merge(expr_df, left_on='sig_id', right_index=True)
    merged_df = merged_df.merge(pert_info[['pert_id', 'canonical_smiles', 'inchi_key']], on='pert_id', how='left')
    
    merged_df = merged_df.rename(columns={
        'canonical_smiles': 'Unnamed: 0',
        'inchi_key': 'InChIKey',
        'pert_id': 'drug_id',
        'cell_id': 'cell_line',
        'pert_idose': 'dose',
        'pert_itime': 'time'
    })
    
    final_df = merged_df.reindex(columns=cp_example.columns.tolist())
    return final_df

def main():
    parser = argparse.ArgumentParser(description="Orchestrator to run Triomphe on test datasets")
    parser.add_argument("--dataset", choices=DATASETS, help="Run on a specific dataset")
    parser.add_argument("--run-all", action="store_true", help="Run all datasets sequentially")
    
    args = parser.parse_args()
    if not (args.dataset or args.run_all):
        parser.print_help()
        sys.exit(1)
        
    datasets = DATASETS if args.run_all else [args.dataset]
    
    # Create outputs directory
    os.makedirs(OUTPUTS_DIR, exist_ok=True)
    
    # Change working directory to Triomphe structure_generator directory context
    original_cwd = os.getcwd()
    os.chdir(CODE_DIR)
    
    from SELFIES_VAE_monodirectional import VAE, Encoder, Decoder, Decoder_input
    from SELFIES_onehot_generator import Onehot_Generator
    
    # Onehot to SELFIES function
    def onehot_decode_to_selfies(encoded_data, label):
        label[0] = ""
        restoration_selfies = label[encoded_data[:, 0]]
        for i in list(range(1, encoded_data.shape[1])):
            restoration_selfies += label[encoded_data[:, i]]
        return restoration_selfies
        
    # Convert latent vector to character string decoder
    def VAE_decoder(latent_vec, model):
        latent_vec = latent_vec.detach().float()
        decoder_input = torch.reshape(latent_vec, (-1, 1, 256))
        for sequence_idx in list(range(99)):
            pred_char = model.decoder(decoder_input)[:, -1, :]
            pred_char = torch.reshape(pred_char, (pred_char.size(0), 1, pred_char.size(1)))
            pred_tag = torch.argmax(pred_char, 2)
            new_onehot = torch.zeros(pred_char.size(0), 1, pred_char.size(2))
            for i in list(range(pred_char.size(0))):
                new_onehot[i, 0, pred_tag[i].item()] = 1
            new_input = model.dim_converter(new_onehot)
            decoder_input = torch.cat([decoder_input, new_input], dim=1)
        decoder_output = model.decoder(decoder_input)
        X_pred = torch.argmax(decoder_output, 2).detach().numpy()
        return X_pred
        
    print("Loading LINCS GSE92742 raw metadata...")
    lincs_gse92742_dir = os.path.join(DATA_DIR, "LINCS/GSE92742")
    sig_info_path = os.path.join(lincs_gse92742_dir, "GSE92742_Broad_LINCS_sig_info.txt.gz")
    gene_info_path = os.path.join(lincs_gse92742_dir, "GSE92742_Broad_LINCS_gene_info.txt.gz")
    pert_info_path = os.path.join(lincs_gse92742_dir, "GSE92742_Broad_LINCS_pert_info.txt.gz")
    cp_example_path = os.path.join(TRIOMPHE_ROOT, "software/data/L1000_ideal_cp_example.csv")
    gctx_path = os.path.join(lincs_gse92742_dir, "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx")
    
    sig_info = pd.read_csv(sig_info_path, sep="\t", low_memory=False)
    gene_info = pd.read_csv(gene_info_path, sep="\t")
    pert_info = pd.read_csv(pert_info_path, sep="\t")
    cp_example = pd.read_csv(cp_example_path)
    
    landmark_genes = gene_info[gene_info['pr_is_lm'].astype(int) == 1]['pr_gene_id'].astype(str).tolist()
    
    print("Loading SELFIES VAE model...")
    vae = VAE(Encoder(), Decoder(), Decoder_input())
    vae.load_state_dict(torch.load("../../model/SELFIES_VAE_monodirectional.pth", map_location=torch.device('cpu')))
    vae.eval()
    
    label = list(np.array(pd.read_csv("../../model/VAE_SELFIES_label.csv", index_col=0)["Symbol"]))
    generated_num = 30
    
    for d in datasets:
        src_path = os.path.join(TRIOMPHE_TEST_DIR, f"{d}_pert.csv")
        if not os.path.exists(src_path):
            print(f"Error: Preprocessed file {src_path} does not exist. Run preprocessing first.")
            continue
            
        cp_file, pert_file, txt_file, mols_file = get_names(d)
        
        print(f"\n[{d}] Preparing Triomphe input dataset...")
        # Copy to Triomphe/software/data/{pert_file}
        tgt_in_path = os.path.join(TRIOMPHE_ROOT, "software/data", pert_file)
        try:
            shutil.copy2(src_path, tgt_in_path)
        except shutil.SameFileError:
            pass
        print(f"[{d}] Triomphe input prepared at {tgt_in_path}")
        
        # -------------------------------------------------------------
        # Source Molecule Selection
        # -------------------------------------------------------------
        print(f"[{d}] Performing source molecule selection...")
        X_cp = extract_cp_profile(d, sig_info, gene_info, pert_info, cp_example, gctx_path, landmark_genes)
        
        # Save X_cp to the standard cp_file path to keep environment status intact
        cp_path = os.path.join(TRIOMPHE_ROOT, "software/data", cp_file)
        X_cp.to_csv(cp_path, index=False)
        
        # Re-load X_cp with index_col=0 to match target vec selection
        X_cp = pd.read_csv(cp_path, index_col=0)
        X_lof = pd.read_csv(tgt_in_path, index_col=0)
        
        target_vec = np.array(X_lof.iloc[:, 8:])
        cp_vec = X_cp.iloc[:, 8:]
        
        # Estimate correlation coefficient
        corr = np.corrcoef(cp_vec, target_vec)
        ko_cp_corrr_matric = corr[cp_vec.shape[0]:, :cp_vec.shape[0]]
        
        # Get top 10 correlations
        min_idx_list = np.argsort(ko_cp_corrr_matric.flatten())
        max_idx_list = min_idx_list[::-1]
        coordinate = np.unravel_index(max_idx_list[:10], ko_cp_corrr_matric.shape)
        
        input_smiles = []
        for i in list(range(len(coordinate[0]))):
            smi_str = cp_vec.index[coordinate[1][i]]
            mol = Chem.MolFromSmiles(smi_str)
            if mol is not None:
                input_smiles.append(Chem.MolToSmiles(mol))
            else:
                input_smiles.append("No SMILES")
                
        top_correlations = [ko_cp_corrr_matric[coordinate[0][i], coordinate[1][i]] for i in list(range(len(input_smiles)))]
        
        txt_out_path = os.path.join(TRIOMPHE_ROOT, "software/result", txt_file)
        os.makedirs(os.path.dirname(txt_out_path), exist_ok=True)
        with open(txt_out_path, "w") as writer:
            writer.write("Correlation\tSMILES\n")
            for i in list(range(len(input_smiles))):
                writer.write(str(top_correlations[i]) + "\t" + input_smiles[i] + "\n")
        print(f"[{d}] Saved correlation to {txt_out_path}")
        
        # -------------------------------------------------------------
        # Novel Chemical Structure Generation
        # -------------------------------------------------------------
        print(f"[{d}] Generating novel molecules using VAE...")
        parent_smiles_list = pd.read_table(txt_out_path)
        parent_smiles_list = parent_smiles_list[parent_smiles_list.SMILES != "No SMILES"].reset_index(drop=True)
        
        # Remove too long SELFIES
        sel_smiles_tag = np.array([len([char for char in sf.split_selfies(sf.encoder(smiles))]) for smiles in parent_smiles_list["SMILES"]]) <= 100
        parent_smiles_list = parent_smiles_list[sel_smiles_tag].reset_index(drop=True)
        
        # Select unique parents to avoid redundant generation
        unique_parents = parent_smiles_list["SMILES"].unique()
        
        pred_smiles = []
        for parent_idx in list(range(min(5, len(unique_parents)))):
            parent_smiles = unique_parents[parent_idx]
            generator = Onehot_Generator([sf.encoder(parent_smiles)], symbol_label=label, max_length=100)
            X_onehot = generator.onehot(mode="TensorFlow")
            latent, logvar = vae.encoder(torch.from_numpy(X_onehot).float())
            latent = latent.detach().numpy()
            logvar = logvar.detach().numpy()
            var = np.exp(logvar)
            
            np.random.seed(0)
            random_latent_vectors = torch.from_numpy(np.random.normal(latent, var * 3, (generated_num, 256)))
            pred_encoded = VAE_decoder(random_latent_vectors, vae)
            selfies_matrix = onehot_decode_to_selfies(pred_encoded, np.array(label).astype(object))
            
            for pred_selfies in selfies_matrix:
                pred_mol = Chem.MolFromSmiles(sf.decoder(pred_selfies))
                if pred_mol is not None:
                    pred_smiles.append(Chem.MolToSmiles(pred_mol))
                else:
                    pred_smiles.append("None")
                    
        generated_all = pd.DataFrame(np.zeros((len(pred_smiles), 2)), columns=["Source_mol", "Generated_mol"])
        source_mols = []
        for parent_idx in list(range(min(5, len(unique_parents)))):
            parent_smiles = unique_parents[parent_idx]
            source_mols.extend([parent_smiles] * generated_num)
            
        generated_all["Source_mol"] = source_mols
        generated_all["Generated_mol"] = pred_smiles
        
        mols_out_path = os.path.join(TRIOMPHE_ROOT, "software/result", mols_file)
        generated_all.to_csv(mols_out_path, index=False)
        print(f"[{d}] Saved generated molecules to {mols_out_path}")
        
        # Save output to TheraFlow outputs folder
        df_clean = generated_all[['Generated_mol']].rename(columns={'Generated_mol': 'SMILES'})
        df_clean = df_clean[~df_clean['SMILES'].isin(["None", "No SMILES", ""])]
        dest_path = os.path.join(OUTPUTS_DIR, f"{d}.csv")
        df_clean.to_csv(dest_path, index=False)
        print(f"[{d}] Saved clean outputs to: {dest_path}")
        
    os.chdir(original_cwd)

if __name__ == "__main__":
    main()
