import os
import subprocess


def run_disease_reversal_pipeline(
    targets_csv, 
    potency_csv, 
    output_reversal_csv, 
    ppi_rds="data/ppi.rds", 
    mapping_file="data/chembl_uniprot_mapping.txt", 
    lincs_meta="Omics2Drug/data/lincs_gene_meta.csv", 
    test_datasets_dir="Omics2Drug/data/test_datasets", 
    workers=4
):
    """Launches the Rscript disease_reversal.R to compute DTSEA-based disease reversal scores."""
    # Locate the R script adjacent to this Python script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    r_script_path = os.path.join(script_dir, "disease_reversal.R")
    
    if not os.path.exists(r_script_path):
        raise FileNotFoundError(f"R script not found at {r_script_path}")
        
    cmd = [
        "Rscript", r_script_path,
        "--targets_csv", str(targets_csv),
        "--potency_csv", str(potency_csv),
        "--output_reversal_csv", str(output_reversal_csv),
        "--ppi_rds", str(ppi_rds),
        "--mapping_file", str(mapping_file),
        "--lincs_meta", str(lincs_meta),
        "--test_datasets_dir", str(test_datasets_dir),
        "--workers", str(workers)
    ]
    
    print(f"Launching R script: {' '.join(cmd)}")
    
    # Run command and capture output
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print("R script execution failed!")
        print("STDOUT:")
        print(result.stdout)
        print("STDERR:")
        print(result.stderr)
        raise RuntimeError(f"R script failed with exit code {result.returncode}")
        
    print("R script executed successfully.")
    print("R Output:")
    print(result.stdout)
    
    return output_reversal_csv
