#!/usr/bin/env python3
"""
Raw Data Access Script for LINCS and GEO Datasets
=================================================
This script facilitates the programmatic download and extraction of the raw 
LINCS Phase 1 and Phase 2 datasets from the NCBI GEO FTP server.

Datasets fetched:
  1. GSE70138 (Phase 2 LINCS Level 5 COMPZ, ~5.6 GB download):
     - GCTX: GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz
     - Sig Info: GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz
     - Gene Info: GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz

  2. GSE92742 (Phase 1 LINCS Level 3 INF, ~53 GB download):
     - GCTX: GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx.gz
     - Inst Info: GSE92742_Broad_LINCS_inst_info.txt.gz
     - Gene Info: GSE92742_Broad_LINCS_gene_info.txt.gz

Usage:
  python raw_data_access.py --dry-run
  python raw_data_access.py --gse70138
"""

import os
import sys
import gzip
import shutil
import argparse
import urllib.request

# Base paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(BASE_DIR)
if not os.path.exists(os.path.join(PARENT_DIR, "SmilesGEN")) and os.path.exists(os.path.join(os.path.dirname(PARENT_DIR), "SmilesGEN")):
    MODEL_ROOT = os.path.dirname(PARENT_DIR)
else:
    MODEL_ROOT = PARENT_DIR

# URLs on NCBI GEO supplementary FTP
GEO_URLS = {
    "GSE70138": {
        "dest_dir": os.path.join(MODEL_ROOT, "SmilesGEN/datasets/LINCS/raw"),
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
        "dest_dir": os.path.join(MODEL_ROOT, "GExMolGen/processing_test"),
        "files": {
            "GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx.gz":
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx.gz",
            "GSE92742_Broad_LINCS_inst_info.txt.gz":
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_inst_info.txt.gz",
            "GSE92742_Broad_LINCS_gene_info.txt.gz":
                "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_gene_info.txt.gz"
        }
    }
}

def show_progress(count, block_size, total_size):
    """Callback for urllib retrieve to print progress bar."""
    if total_size > 0:
        percent = int(count * block_size * 100 / total_size)
        percent = min(percent, 100)
        downloaded = count * block_size / (1024 * 1024)
        total_mb = total_size / (1024 * 1024)
        sys.stdout.write(f"\r  Downloading... {percent}% ({downloaded:.1f}/{total_mb:.1f} MB)")
        sys.stdout.flush()

def download_and_extract(accession, dry_run=False):
    print(f"\n==================================================")
    print(f"Dataset Access: {accession}")
    print(f"==================================================")
    
    meta = GEO_URLS[accession]
    dest_dir = meta["dest_dir"]
    
    if not dry_run:
        os.makedirs(dest_dir, exist_ok=True)
        
    for filename, url in meta["files"].items():
        dest_gz_path = os.path.join(dest_dir, filename)
        unzipped_name = filename[:-3]  # Strip ".gz"
        dest_unzipped_path = os.path.join(dest_dir, unzipped_name)
        
        # Check if already present unzipped
        if os.path.exists(dest_unzipped_path):
            print(f"  [Skipped] {unzipped_name} already exists at {dest_unzipped_path}")
            continue
            
        print(f"  Target File: {unzipped_name}")
        print(f"  Source URL:  {url}")
        
        if dry_run:
            print(f"  [Dry Run] Would download to: {dest_gz_path}")
            continue
            
        print(f"  Downloading to temporary archive: {dest_gz_path}...")
        try:
            urllib.request.urlretrieve(url, dest_gz_path, show_progress)
            print("\n  Download complete!")
        except Exception as e:
            print(f"\n  Error downloading {url}: {e}")
            continue
            
        print(f"  Extracting {filename} -> {unzipped_name}...")
        try:
            with gzip.open(dest_gz_path, 'rb') as f_in:
                with open(dest_unzipped_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            print(f"  Extraction successful! File ready at: {dest_unzipped_path}")
            # Clean up the .gz archive to save disk space
            os.remove(dest_gz_path)
            print(f"  Cleaned up temporary archive {filename}")
        except Exception as e:
            print(f"  Error extracting archive: {e}")
            if os.path.exists(dest_gz_path):
                os.remove(dest_gz_path)

def main():
    parser = argparse.ArgumentParser(description="Programmatic Download & Setup of Raw LINCS/GEO Datasets")
    parser.add_argument("--gse70138", action="store_true", help="Download GSE70138 (Phase 2 Level 5 COMPZ, ~5.6 GB)")
    parser.add_argument("--gse92742", action="store_true", help="Download GSE92742 (Phase 1 Level 3 INF, ~53 GB)")
    parser.add_argument("--download-all", action="store_true", help="Download both datasets (~59 GB total)")
    parser.add_argument("--dry-run", action="store_true", help="Print download URLs and locations without fetching")
    
    args = parser.parse_args()
    
    if not (args.gse70138 or args.gse92742 or args.download_all or args.dry-run):
        parser.print_help()
        sys.exit(0)
        
    dry_run = args.dry-run
    
    if args.download_all or args.gse70138:
        download_and_extract("GSE70138", dry_run)
    if args.download_all or args.gse92742:
        download_and_extract("GSE92742", dry_run)
        
    print("\nProcessing complete!")

if __name__ == "__main__":
    main()
