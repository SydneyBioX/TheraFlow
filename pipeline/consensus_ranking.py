import os
import subprocess


def run_consensus_ranking_pipeline(input_csv, output_csv):
    """Launches the Rscript consensus_ranking.R to perform Pareto/RRA ranking and select top molecules."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    r_script_path = os.path.join(script_dir, "consensus_ranking.R")
    
    if not os.path.exists(r_script_path):
        raise FileNotFoundError(f"R script not found at {r_script_path}")
        
    cmd = [
        "Rscript", r_script_path,
        "--input_csv", str(input_csv),
        "--output_csv", str(output_csv)
    ]
    
    print(f"Launching consensus ranking R script: {' '.join(cmd)}")
    
    # Run command and capture output
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print("Consensus ranking R script execution failed!")
        print("STDOUT:")
        print(result.stdout)
        print("STDERR:")
        print(result.stderr)
        raise RuntimeError(f"R script failed with exit code {result.returncode}")
        
    print("Consensus ranking R script executed successfully.")
    print("R Output:")
    print(result.stdout)
    
    return output_csv
