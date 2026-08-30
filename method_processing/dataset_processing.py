#!/usr/bin/env python3
"""
Unified Download and Preprocessing Script
=========================================
Handles downloading of LINCS and GEO datasets, parses and processes LINCS
data in Python, and invokes the R script for edgeR/limma disease data processing.
"""

import os
import sys
import gzip
import shutil
import pickle
import argparse
import subprocess
import urllib.request
import pandas as pd
import numpy as np

# Base paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(BASE_DIR)
if not os.path.exists(os.path.join(PARENT_DIR, "SmilesGEN")) and os.path.exists(os.path.join(os.path.dirname(PARENT_DIR), "SmilesGEN")):
    MODEL_ROOT = os.path.dirname(PARENT_DIR)
else:
    MODEL_ROOT = PARENT_DIR

DATA_DIR = os.path.join(BASE_DIR, "data")

GEO_URLS = {
    "GSE70138": {
        "dest_dir": os.path.join(DATA_DIR, "LINCS/GSE70138"),
        "files": {
            "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz": 
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz",
            "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz":
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz",
            "GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz":
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz"
        }
    },
    "GSE92742": {
        "dest_dir": os.path.join(DATA_DIR, "LINCS/GSE92742"),
        "files": {
            "GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx.gz":
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx.gz",
            "GSE92742_Broad_LINCS_inst_info.txt.gz":
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_inst_info.txt.gz",
            "GSE92742_Broad_LINCS_gene_info.txt.gz":
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_gene_info.txt.gz"
        }
    },
    "GSE143233": {
        "dest_dir": os.path.join(DATA_DIR, "bulk_transcriptomics_data/GSE143233"),
        "files": {
            "GSE143233_patients_normalHCC_and_sorafenib_resistantHCC_genes_fpkm_tracking.txt.gz":
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE143nnn/GSE143233/suppl/GSE143233_patients_normalHCC_and_sorafenib_resistantHCC_genes_fpkm_tracking.txt.gz"
        }
    },
    "GSE110716": {
        "dest_dir": os.path.join(DATA_DIR, "bulk_transcriptomics_data/GSE110716"),
        "files": {
            "GSE110716_RAW.tar":
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE110nnn/GSE110716/suppl/GSE110716_RAW.tar",
            "GSE110716_series_matrix.txt.gz":
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE110nnn/GSE110716/matrix/GSE110716_series_matrix.txt.gz"
        }
    },
    "GSE55296": {
        "dest_dir": os.path.join(DATA_DIR, "bulk_transcriptomics_data/GSE55296"),
        "files": {
            "GSE55296_count_data.txt.gz":
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE55nnn/GSE55296/suppl/GSE55296_count_data.txt.gz",
            "GSE55296_processed_data_readme.txt":
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE55nnn/GSE55296/suppl/GSE55296_processed_data_readme.txt"
        }
    },
    "GSE236302": {
        "dest_dir": os.path.join(DATA_DIR, "bulk_transcriptomics_data/GSE236302"),
        "files": {
            "GSE236302_raw_counts.csv.gz":
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE236nnn/GSE236302/suppl/GSE236302_raw_counts.csv.gz",
            "GSE236302_sample_metadata.csv.gz":
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE236nnn/GSE236302/suppl/GSE236302_sample_metadata.csv.gz"
        }
    }
}

def show_progress(count, block_size, total_size):
    if total_size > 0:
        percent = min(int(count * block_size * 100 / total_size), 100)
        downloaded = count * block_size / (1024 * 1024)
        total_mb = total_size / (1024 * 1024)
        sys.stdout.write(f"\r  Downloading... {percent}% ({downloaded:.1f}/{total_mb:.1f} MB)")
        sys.stdout.flush()

def download_dataset(accession):
    print(f"\nDownloading dataset: {accession}")
    meta = GEO_URLS[accession]
    dest_dir = meta["dest_dir"]
    os.makedirs(dest_dir, exist_ok=True)
    
    for filename, url in meta["files"].items():
        dest_path = os.path.join(dest_dir, filename)
        if os.path.exists(dest_path) or (filename.endswith(".gz") and os.path.exists(dest_path[:-3])):
            print(f"  [Skipped] {filename} already exists.")
            continue
            
        print(f"  Fetching: {url}")
        try:
            urllib.request.urlretrieve(url, dest_path, show_progress)
            print("\n  Download complete.")
        except Exception as e:
            print(f"\n  Error downloading {filename}: {e}")
            if os.path.exists(dest_path):
                try:
                    os.remove(dest_path)
                except:
                    pass
            continue
            
        if filename.endswith(".gz"):
            try:
                # Only unzip GCTX files, read others directly using gzip reader or R
                if "Level5_COMPZ" in filename or "Level3_INF" in filename:
                    unzipped = dest_path[:-3]
                    print(f"  Unzipping: {filename} -> {os.path.basename(unzipped)}")
                    with gzip.open(dest_path, 'rb') as f_in:
                        with open(unzipped, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    os.remove(dest_path)
            except Exception as e:
                print(f"  Error unzipping {filename}: {e}")

def extract_gse110716_metadata():
    metadata_path = os.path.join(DATA_DIR, "bulk_transcriptomics_data", "GSE110716_clinical_metadata.csv")
    if os.path.exists(metadata_path):
        return
    print("\nExtracting GSE110716 clinical metadata using R...")
    r_cmd = (
        'library(GEOquery); '
        'gse <- getGEO("GSE110716", GSEMatrix=TRUE, destdir="' + os.path.join(DATA_DIR, "bulk_transcriptomics_data") + '"); '
        'write.csv(pData(gse[[1]]), "' + metadata_path + '")'
    )
    try:
        subprocess.run(["Rscript", "-e", r_cmd], check=True)
    except Exception as e:
        print(f"Error extracting metadata: {e}")

def preprocess_rest():
    print("==================================================")
    print("Starting Rest Methods Preprocessing...")
    print("==================================================")
    try:
        from cmapPy.pandasGEXpress.parse import parse
    except ImportError:
        print("Error: cmapPy is required for rest methods preprocessing.")
        return
        
    raw_dir = os.path.join(DATA_DIR, "LINCS/GSE70138")
    sig_info_path = os.path.join(raw_dir, "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz")
    gene_info_path = os.path.join(raw_dir, "GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz")
    gctx_path = os.path.join(raw_dir, "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx")
    output_dir = os.path.join(DATA_DIR, "test_datasets/rest_methods")
    os.makedirs(output_dir, exist_ok=True)
    
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
        
        # # Copy to model folder for standard run compatibility
        # tgt_dir = os.path.join(MODEL_ROOT, "SmilesGEN/datasets/test_protein")
        # os.makedirs(tgt_dir, exist_ok=True)
        # # Save under exact model hardcoded paths
        # df.to_csv(os.path.join(tgt_dir, f"{cell}_{pert}.csv" if not mean else f"{cell}_{pert}_mean.csv"), index=False)
        # if pert == "AKT-inhibitor-1-2":
        #     df.to_csv(os.path.join(tgt_dir, f"{cell}_AKT1.csv" if not mean else f"{cell}_AKT1_mean.csv"), index=False)
            
    get_and_save("MCF7", "AKT-inhibitor-1-2", "lincs1000_1.csv", mean=True)
    get_and_save("A375", "venlafaxine", "lincs1000_2.csv", mean=True)
    print("Rest Methods Preprocessing Completed Successfully!\n")

def preprocess_gexmolgen():
    print("==================================================")
    print("Starting GExMolGen Preprocessing...")
    print("==================================================")
    try:
        from cmapPy.pandasGEXpress.parse import parse
    except ImportError:
        print("Error: cmapPy is required for GExMolGen preprocessing.")
        return
        
    gex_dir = os.path.join(DATA_DIR, "LINCS/GSE92742")
    gctx_path = os.path.join(gex_dir, "GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx")
    inst_info_path = os.path.join(gex_dir, "GSE92742_Broad_LINCS_inst_info.txt.gz")
    gene_info_path = os.path.join(gex_dir, "GSE92742_Broad_LINCS_gene_info.txt.gz")
    ref_path = os.path.join(MODEL_ROOT, "GExMolGen/server_test_ctl_AKT2.csv")
    output_dir = os.path.join(DATA_DIR, "test_datasets/gexmolgen")
    os.makedirs(output_dir, exist_ok=True)
    
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
        
        # Save under standard suffix formats (both suffix 2 and 3 for venlafaxine to match both scripts/notebooks)
        for sfx in ([out_suffix] if out_suffix == "1" else ["2", "3"]):
            pert_out = os.path.join(output_dir, f"lincs1000_{sfx}_pert.csv")
            ctl_out = os.path.join(output_dir, f"lincs1000_{sfx}_ctl.csv")
            df_pert.to_csv(pert_out)
            df_ctl.to_csv(ctl_out)
            print(f"  Saved matched pair: {pert_out} & {ctl_out}")
            
    process_pair("MCF7", "AKT-inhibitor-1-2", "1")
    process_pair("A375", "venlafaxine", "2")
    print("GExMolGen Preprocessing Completed Successfully!\n")

def preprocess_paccmann():
    print("==================================================")
    print("Starting paccmann_rl Preprocessing...")
    print("==================================================")
    try:
        from cmapPy.pandasGEXpress.parse import parse
    except ImportError:
        print("Error: cmapPy is required for paccmann_rl preprocessing.")
        return
        
    paccmann_dir = os.path.join(MODEL_ROOT, "paccmann_rl")
    gene_list_path = os.path.join(DATA_DIR, "2128_genes.pkl")
    raw_dir = os.path.join(DATA_DIR, "LINCS/GSE70138")
    sig_info_path = os.path.join(raw_dir, "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz")
    gene_info_path = os.path.join(raw_dir, "GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz")
    gctx_path = os.path.join(raw_dir, "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx")
    output_dir = os.path.join(DATA_DIR, "test_datasets/paccmann_rl")
    os.makedirs(output_dir, exist_ok=True)
    
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

def preprocess_transgem():
    print("==================================================")
    print("Starting TransGEM Preprocessing...")
    print("==================================================")
    test_pc_path = os.path.join(MODEL_ROOT, "TransGEM/data/test_PC.csv")
    rest_dir = os.path.join(DATA_DIR, "test_datasets/rest_methods")
    output_dir = os.path.join(DATA_DIR, "test_datasets/transgem")
    os.makedirs(output_dir, exist_ok=True)
    
    if not os.path.exists(test_pc_path):
        print(f"Error: TransGEM test_PC.csv not found at {test_pc_path}")
        return
        
    test_pc = pd.read_csv(test_pc_path)
    
    def process_file(src_name, dest_name):
        src_path = os.path.join(rest_dir, src_name)
        if not os.path.exists(src_path):
            print(f"Warning: file {src_path} not found. Skipping.")
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
        tgt_path = os.path.join(MODEL_ROOT, "TransGEM/data", f"test_{dest_name}")
        df.to_csv(tgt_path, index=False)
        
    process_file("lincs1000_1.csv", "lincs1000_1.csv")
    process_file("lincs1000_2.csv", "lincs1000_2.csv")
    process_file("cardiomyopathy.csv", "cardiomyopathy.csv")
    process_file("glaucoma.csv", "glaucoma.csv")
    process_file("glaucoma_bulk.csv", "glaucoma_bulk.csv")
    process_file("hcc.csv", "hcc.csv")
    process_file("parkinsons.csv", "parkinsons.csv")
    print("TransGEM Preprocessing Completed Successfully!\n")

def run_disease_processing_r():
    print("==================================================")
    print("Running Disease Preprocessing (R Script)...")
    print("==================================================")
    r_script = os.path.join(BASE_DIR, "dataset_processing.R")
    if not os.path.exists(r_script):
        print(f"Error: R script not found at {r_script}")
        return
    subprocess.run(["Rscript", r_script], check=True)

def main():
    parser = argparse.ArgumentParser(description="Unified Download and Preprocessing Script for omics2mol Datasets")
    parser.add_argument("--download", action="store_true", help="Download LINCS and GEO Disease datasets")
    parser.add_argument("--rest-methods", action="store_true", help="Preprocess rest methods LINCS")
    parser.add_argument("--gexmolgen", action="store_true", help="Preprocess GExMolGen LINCS")
    parser.add_argument("--paccmann", action="store_true", help="Preprocess paccmann_rl LINCS")
    parser.add_argument("--diseases", action="store_true", help="Preprocess Disease datasets via R script")
    parser.add_argument("--transgem", action="store_true", help="Preprocess TransGEM for LINCS and Disease")
    parser.add_argument("--run-all", action="store_true", help="Run all download and preprocessing steps")
    
    args = parser.parse_args()
    
    if not (args.download or args.rest_methods or args.gexmolgen or args.paccmann or args.diseases or args.transgem or args.run_all):
        parser.print_help()
        sys.exit(1)
        
    if args.run_all or args.download:
        for acc in GEO_URLS.keys():
            download_dataset(acc)
        extract_gse110716_metadata()
        
    if args.run_all or args.rest_methods:
        preprocess_rest()
        
    if args.run_all or args.gexmolgen:
        preprocess_gexmolgen()
        
    if args.run_all or args.paccmann:
        preprocess_paccmann()
        
    if args.run_all or args.diseases:
        run_disease_processing_r()
        
    if args.run_all or args.transgem:
        preprocess_transgem()

if __name__ == "__main__":
    main()
