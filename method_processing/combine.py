#!/usr/bin/env python3
import os
import glob
import pandas as pd

def combine_outputs():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    outputs_dir = os.path.join(base_dir, "outputs")
    data_dir = os.path.join(base_dir, "data")
    out_csv = os.path.join(data_dir, "processed_data.csv")
    
    os.makedirs(data_dir, exist_ok=True)
    
    records = []
    
    # Iterate through all method directories
    methods = sorted(os.listdir(outputs_dir))
    for method in methods:
        method_path = os.path.join(outputs_dir, method)
        if not os.path.isdir(method_path):
            continue
        
        csv_files = sorted(glob.glob(os.path.join(method_path, "*.csv")))
        for csv_file in csv_files:
            test_dataset = os.path.splitext(os.path.basename(csv_file))[0]
            
            try:
                df = pd.read_csv(csv_file)
                # Find SMILES column (case-insensitive)
                smiles_col = None
                for col in df.columns:
                    if col.lower() == "smiles":
                        smiles_col = col
                        break
                
                if smiles_col is None:
                    print(f"Warning: No SMILES column in {csv_file}")
                    continue
                
                smiles_series = df[smiles_col].dropna().astype(str).str.strip()
                smiles_series = smiles_series[smiles_series != ""]
                
                for s in smiles_series:
                    records.append({
                        "SMILES": s,
                        "method": method,
                        "test_dataset": test_dataset
                    })
            except Exception as e:
                print(f"Error reading {csv_file}: {e}")
                
    combined_df = pd.DataFrame(records, columns=["SMILES", "method", "test_dataset"])
    combined_df.to_csv(out_csv, index=False)
    
    print(f"Successfully created: {out_csv}")
    print(f"Total rows: {len(combined_df)}")
    print("\nBreakdown by method:")
    print(combined_df["method"].value_counts())
    print("\nBreakdown by test_dataset:")
    print(combined_df["test_dataset"].value_counts())
    print("\nSample records:")
    print(combined_df.head(10))

if __name__ == "__main__":
    combine_outputs()
