#!/usr/bin/env python3
"""
Unified Preprocessing Script for omics2mol Datasets
===================================================
This script compiles preprocessing logic from:
  - SmilesGEN: SmilesGEN/datasets/data_processing_smilesgen.ipynb
  - GExMolGen: GExMolGen/processing_test/processing.ipynb
  - TransGEM: TransGEM/processing_datasets.ipynb
  - paccmann_rl: paccmann_rl/processing_data.ipynb

Usage:
  python run_data_preprocessing.py --run-all
  python run_data_preprocessing.py --smilesgen --transgem
"""

import os
import sys
import pickle
import argparse
import pandas as pd
import numpy as np

# Ensure parent directory is accessible for modules
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(BASE_DIR)
if not os.path.exists(os.path.join(PARENT_DIR, "SmilesGEN")) and os.path.exists(os.path.join(os.path.dirname(PARENT_DIR), "SmilesGEN")):
    PARENT_DIR = os.path.dirname(PARENT_DIR)

def preprocess_smilesgen():
    print("==================================================")
    print("Starting SmilesGEN Preprocessing...")
    print("==================================================")
    try:
        from cmapPy.pandasGEXpress.parse import parse
    except ImportError:
        print("Error: cmapPy is required for SmilesGEN preprocessing. Please install it.")
        return
    
    raw_dir = os.path.join(PARENT_DIR, "SmilesGEN/datasets/LINCS/raw")
    sig_info_path = os.path.join(raw_dir, "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz")
    gene_info_path = os.path.join(raw_dir, "GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz")
    gctx_path = os.path.join(raw_dir, "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx")
    output_dir = os.path.join(PARENT_DIR, "data/test_datasets/smilesgen")
    os.makedirs(output_dir, exist_ok=True)
    
    # Check files
    for path in [sig_info_path, gene_info_path, gctx_path]:
        if not os.path.exists(path):
            print(f"Error: Required file not found: {path}")
            return

    sig_info = pd.read_csv(sig_info_path, sep="\t", low_memory=False)
    gene_info = pd.read_csv(gene_info_path, sep="\t")
    landmark_genes = gene_info[gene_info['pr_is_lm'].astype(int) == 1]['pr_gene_id'].astype(str).tolist()
    
    def get_and_save(cell, pert, name, mean):
        sig_ids = sig_info[(sig_info['cell_id'] == cell) & (sig_info['pert_iname'] == pert)]['sig_id'].tolist()
        if not sig_ids:
            print(f"Warning: No signatures found for {cell} + {pert}")
            return
        gctoo = parse(gctx_path, rid=landmark_genes, cid=sig_ids)
        df = gctoo.data_df.T
        if mean:
            df = df.mean(axis=0).to_frame().T
            df.insert(0, 'sig_id', f"{pert}_mean")
        out_file = os.path.join(output_dir, name)
        df.to_csv(out_file, index=False)
        print(f"  Saved: {out_file}")
        
        # Also copy to SmilesGEN folder for standard run
        tgt_dir = os.path.join(PARENT_DIR, "SmilesGEN/datasets/test_protein")
        os.makedirs(tgt_dir, exist_ok=True)
        df.to_csv(os.path.join(tgt_dir, f"{cell}_{pert}.csv" if not mean else f"{cell}_{pert}_mean.csv"), index=False)
        
    get_and_save("MCF7", "AKT-inhibitor-1-2", "Lincs1000_1.csv", mean=False)
    get_and_save("MCF7", "AKT-inhibitor-1-2", "Lincs1000_1_m.csv", mean=True)
    get_and_save("A375", "venlafaxine", "Lincs1000_2.csv", mean=False)
    get_and_save("A375", "venlafaxine", "Lincs1000_2_m.csv", mean=True)
    print("SmilesGEN Preprocessing Completed Successfully!\n")

def preprocess_gexmolgen():
    print("==================================================")
    print("Starting GExMolGen Preprocessing...")
    print("==================================================")
    try:
        from cmapPy.pandasGEXpress.parse import parse
    except ImportError:
        print("Error: cmapPy is required for GExMolGen preprocessing.")
        return
        
    gex_dir = os.path.join(PARENT_DIR, "GExMolGen/processing_test")
    gctx_path = os.path.join(gex_dir, "GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx")
    inst_info_path = os.path.join(gex_dir, "GSE92742_Broad_LINCS_inst_info.txt.gz")
    gene_info_path = os.path.join(gex_dir, "GSE92742_Broad_LINCS_gene_info.txt.gz")
    ref_path = os.path.join(PARENT_DIR, "GExMolGen/server_test_ctl_AKT2.csv")
    output_dir = os.path.join(PARENT_DIR, "data/test_datasets/gexmolgen")
    os.makedirs(output_dir, exist_ok=True)
    
    for path in [gctx_path, inst_info_path, gene_info_path, ref_path]:
        if not os.path.exists(path):
            print(f"Error: Required file not found: {path}")
            return
            
    ref_columns = pd.read_csv(ref_path, index_col=0).columns.tolist()
    inst_info = pd.read_csv(inst_info_path, sep="\t", low_memory=False)
    gene_info = pd.read_csv(gene_info_path, sep="\t", low_memory=False)
    
    def process_pair(cell_line, perturbation, out_suffix):
        test_meta = inst_info[(inst_info['cell_id'] == cell_line) & (inst_info['pert_iname'] == perturbation)]
        if test_meta.empty:
            print(f"Warning: No samples found for {perturbation} in {cell_line}")
            return
        relevant_plates = test_meta['rna_plate'].unique()
        ctl_meta = inst_info[(inst_info['cell_id'] == cell_line) & 
                             (inst_info['pert_type'].isin(['ctl_vehicle', 'ctl_vector'])) & 
                             (inst_info['rna_plate'].isin(relevant_plates))]
        
        all_ids = list(test_meta['inst_id'].unique()) + list(ctl_meta['inst_id'].unique())
        raw_data = parse(gctx_path, cid=all_ids).data_df
        gene_map = dict(zip(gene_info['pr_gene_id'].astype(str), gene_info['pr_gene_symbol']))
        raw_data.index = raw_data.index.astype(str).map(gene_map)
        raw_data = raw_data[~raw_data.index.duplicated(keep='first')].T
        
        pert_rows = []
        ctl_rows = []
        for plate in relevant_plates:
            p_ids = test_meta[test_meta['rna_plate'] == plate]['inst_id'].tolist()
            c_ids = ctl_meta[ctl_meta['rna_plate'] == plate]['inst_id'].tolist()
            if not p_ids or not c_ids:
                continue
            ctl_mean = raw_data.loc[c_ids].mean(axis=0)
            for p_id in p_ids:
                pert_rows.append(raw_data.loc[p_id])
                ctl_rows.append(ctl_mean)
                
        df_pert = pd.DataFrame(pert_rows).reindex(columns=ref_columns, fill_value=0)
        df_ctl = pd.DataFrame(ctl_rows).reindex(columns=ref_columns, fill_value=0)
        
        pert_out = os.path.join(output_dir, f"Lincs1000_{out_suffix}_pert.csv")
        ctl_out = os.path.join(output_dir, f"Lincs1000_{out_suffix}_ctl.csv")
        df_pert.to_csv(pert_out)
        df_ctl.to_csv(ctl_out)
        print(f"  Saved matched pair: {pert_out} & {ctl_out}")
        
    process_pair("MCF7", "AKT-inhibitor-1-2", "1")
    process_pair("A375", "venlafaxine", "2")
    print("GExMolGen Preprocessing Completed Successfully!\n")

def preprocess_transgem():
    print("==================================================")
    print("Starting TransGEM Preprocessing...")
    print("==================================================")
    test_pc_path = os.path.join(PARENT_DIR, "TransGEM/data/test_PC.csv")
    smilesgen_dir = os.path.join(PARENT_DIR, "data/test_datasets/smilesgen")
    output_dir = os.path.join(PARENT_DIR, "data/test_datasets/transgem")
    os.makedirs(output_dir, exist_ok=True)
    
    if not os.path.exists(test_pc_path):
        print(f"Error: TransGEM test_PC.csv not found at {test_pc_path}")
        return
        
    test_pc = pd.read_csv(test_pc_path)
    
    def process_file(src_name, dest_name):
        src_path = os.path.join(smilesgen_dir, src_name)
        if not os.path.exists(src_path):
            print(f"Warning: SmilesGEN file {src_path} not found. Running SmilesGEN preprocessing first is recommended.")
            return
        df = pd.read_csv(src_path).transpose().reset_index()
        df.columns = ['gene', 'gene_e']
        df = df[1:]
        df['gene'] = df['gene'].astype('int64')
        df['gene_e'] = df['gene_e'].astype('float64')
        df = test_pc[['gene']].merge(df, on='gene', how='inner')
        
        dest_path = os.path.join(output_dir, dest_name)
        df.to_csv(dest_path, index=False)
        print(f"  Saved TransGEM: {dest_path}")
        
        # Also copy to TransGEM dataset folder for standard run
        tgt_path = os.path.join(PARENT_DIR, "TransGEM/data", f"test_{dest_name}")
        df.to_csv(tgt_path, index=False)
        
    process_file("Lincs1000_1_m.csv", "Lincs1000_1.csv")
    process_file("Lincs1000_2_m.csv", "Lincs1000_2.csv")
    process_file("glaucoma.csv", "glaucoma.csv")
    print("TransGEM Preprocessing Completed Successfully!\n")

def preprocess_paccmann():
    print("==================================================")
    print("Starting paccmann_rl Preprocessing...")
    print("==================================================")
    try:
        from cmapPy.pandasGEXpress.parse import parse
    except ImportError:
        print("Error: cmapPy is required for paccmann_rl preprocessing.")
        return
        
    paccmann_dir = os.path.join(PARENT_DIR, "paccmann_rl")
    gene_list_path = os.path.join(paccmann_dir, "data/2128_genes.pkl")
    
    raw_dir = os.path.join(PARENT_DIR, "SmilesGEN/datasets/LINCS/raw")
    sig_info_path = os.path.join(raw_dir, "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz")
    gene_info_path = os.path.join(raw_dir, "GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz")
    gctx_path = os.path.join(raw_dir, "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx")
    output_dir = os.path.join(PARENT_DIR, "data/test_datasets/paccmann_rl")
    os.makedirs(output_dir, exist_ok=True)
    
    if not os.path.exists(gene_list_path):
        print(f"Error: paccmann_rl genes pkl file not found: {gene_list_path}")
        return
        
    with open(gene_list_path, "rb") as f:
        paccmann_gene_symbols = pickle.load(f)
        
    sig_info = pd.read_csv(sig_info_path, sep="\t", low_memory=False)
    gene_info = pd.read_csv(gene_info_path, sep="\t")
    
    symbol_to_entrez = gene_info.set_index('pr_gene_symbol')['pr_gene_id'].astype(str).to_dict()
    paccmann_entrez = [symbol_to_entrez[sym] for sym in paccmann_gene_symbols if sym in symbol_to_entrez]
    entrez_to_symbol = gene_info.set_index('pr_gene_id').astype(str)['pr_gene_symbol'].to_dict()
    
    def process_paccmann(cell_line, perturbation, dest_name):
        cell_line_sigs = sig_info[(sig_info['cell_id'] == cell_line) & (sig_info['pert_iname'] == perturbation)]['sig_id'].tolist()
        if not cell_line_sigs:
            print(f"Warning: No signatures found for {cell_line} + {perturbation}")
            return
        gctoo = parse(gctx_path, rid=paccmann_entrez, cid=cell_line_sigs)
        df = gctoo.data_df.T
        df.columns = df.columns.map(lambda x: entrez_to_symbol.get(int(x), x))
        df = df.reindex(columns=paccmann_gene_symbols, fill_value=0.0)
        df = df.mean(axis=0).to_frame().T
        df.insert(0, 'sig_id', f'{cell_line}_{perturbation}_mean')
        
        dest_path = os.path.join(output_dir, dest_name)
        df.to_csv(dest_path, index=False)
        print(f"  Saved paccmann_rl: {dest_path}")
        
    process_paccmann("MCF7", "AKT-inhibitor-1-2", "lincs1000_1.csv")
    process_paccmann("A375", "venlafaxine", "lincs1000_2.csv")
    print("paccmann_rl Preprocessing Completed Successfully!\n")

def main():
    parser = argparse.ArgumentParser(description="Unified Preprocessing Script for Molecule Generation Datasets")
    parser.add_argument("--smilesgen", action="store_true", help="Preprocess SmilesGEN")
    parser.add_argument("--gexmolgen", action="store_true", help="Preprocess GExMolGen")
    parser.add_argument("--transgem", action="store_true", help="Preprocess TransGEM")
    parser.add_argument("--paccmann", action="store_true", help="Preprocess paccmann_rl")
    parser.add_argument("--run-all", action="store_true", help="Preprocess all datasets")
    
    args = parser.parse_args()
    
    if not (args.smilesgen or args.gexmolgen or args.transgem or args.paccmann or args.run-all):
        parser.print_help()
        sys.exit(1)
        
    if args.run_all or args.smilesgen:
        preprocess_smilesgen()
    if args.run_all or args.gexmolgen:
        preprocess_gexmolgen()
    if args.run_all or args.transgem:
        preprocess_transgem()
    if args.run_all or args.paccmann:
        preprocess_paccmann()

if __name__ == "__main__":
    main()
