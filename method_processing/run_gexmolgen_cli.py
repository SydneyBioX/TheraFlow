import os
import sys
import argparse
import shutil
import pandas as pd

# Ensure GExMolGen is in the path
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(BASE_DIR)
if not os.path.exists(os.path.join(PARENT_DIR, "GExMolGen")) and os.path.exists(os.path.join(os.path.dirname(PARENT_DIR), "GExMolGen")):
    PARENT_DIR = os.path.dirname(PARENT_DIR)
sys.path.append(os.path.join(PARENT_DIR, "GExMolGen"))

from functions import Standard

class GradioFileWrapper:
    def __init__(self, name):
        self.name = name

def main():
    parser = argparse.ArgumentParser(description="Run GExMolGen offline via CLI")
    parser.add_argument("--control", required=True, help="Path to control gene expression CSV")
    parser.add_argument("--pert", required=True, help="Path to perturbation/disease gene expression CSV")
    parser.add_argument("--output", required=True, help="Output path for generated smiles CSV")
    parser.add_argument("--beam_width", type=int, default=100, help="Beam width / number of molecules to generate")
    parser.add_argument("--seed", type=int, default=0, help="Random seed")
    
    args = parser.parse_args()
    
    # Check if input files exist
    if not os.path.exists(args.control):
        print(f"Error: Control file {args.control} not found.")
        sys.exit(1)
    if not os.path.exists(args.pert):
        print(f"Error: Perturbation file {args.pert} not found.")
        sys.exit(1)
        
    print(f"Running GExMolGen standard generation...")
    print(f"  Control: {args.control}")
    print(f"  Pert: {args.pert}")
    print(f"  Beam Width: {args.beam_width}")
    print(f"  Seed: {args.seed}")
    
    control_wrapper = GradioFileWrapper(args.control)
    pert_wrapper = GradioFileWrapper(args.pert)
    
    # Run the generator
    output_temp_csv, _ = Standard(
        control=control_wrapper,
        pert=pert_wrapper,
        beam_width=args.beam_width,
        seed=args.seed
    )
    
    # Copy the temp output to the desired location
    if os.path.exists(output_temp_csv):
        os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
        shutil.copy(output_temp_csv, args.output)
        print(f"Success! Generated molecules saved to {args.output}")
        # Clean up temp file
        os.remove(output_temp_csv)
    else:
        print(f"Error: Expected output file {output_temp_csv} was not created.")
        sys.exit(1)

if __name__ == "__main__":
    main()
