# TheraFlow: Metrics Calculation Pipeline

A modular, robust, and parallelized Python/R pipeline for calculating molecular metrics, safety profiles, ChEMBL target similarities, potency scores, and running DTSEA-based disease reversal.

## Project Layout

When cloned or set up as a standalone repository, TheraFlow files should be structured as follows:

```
your-github-repo/
├── pyproject.toml              # Modern Python packaging configuration
├── README.md                   # Package documentation
├── __init__.py                 # Package initializer
├── utils.py                    # Output merging and RDKit 3D SDF conformation generation
├── molecular_properties.py     # Descriptors (QED, SA, LogP, TPSA, CSP3), PAINS, similarities
├── posebuster_analysis.py      # PoseBusters validation checks
├── potency_chembl.py           # Parallel ChEMBL API target and activity fetches
├── disease_reversal.py         # Python wrapper invoking the R DTSEA script
├── disease_reversal.R          # R script for biomaRt mapping and DTSEA calculations
├── consensus_ranking.py        # Python wrapper invoking the R consensus ranking script
├── consensus_ranking.R         # R script for Pareto front and RRA consensus ranking
├── main_pipeline.py            # Global pipeline orchestrator (CLI entrypoint)
└── data/                       # Directory containing reference and mapping datasets
    ├── ppi.rds                 # PPI network RDS file (originally from STRING DB, score > 700)
    ├── chembl_uniprot_mapping.txt # Target descriptions to UniProt IDs, can use BioMart in R to get the mappings
    ├── lincs_gene_meta.csv     # LINCS probe IDs to gene symbols, can be downloaded from LINCS1000
    └── test_datasets/          # Folder containing signature CSV files (e.g. glaucoma.csv)
```

---

## Installation

### Dependencies

#### Python Packages
* `pandas`
* `numpy`
* `rdkit`
* `joblib`
* `tqdm`
* `posebusters`
* `chembl-webresource-client`

#### R Packages
* `DTSEA`
* `biomaRt`
* `org.Hs.eg.db` (fallback)
* `limma`
* `data.table`
* `dplyr`
* `tidyr`
* `igraph`
* `furrr`

### Package Installation

Install the package in editable mode from the root directory of your cloned repository:

```bash
pip install -e .
```

Once installed, you can launch the pipeline directly from your shell command line using `theraflow`.

### Required Reference Datasets (User Supplied)

Because reference files are large, they are **not** included in the repository. You must create a `data/` directory in the repository root and place the following files there before running the pipeline:

1. **`data/ppi.rds`**: A Protein-Protein Interaction (PPI) network saved as an R `igraph` RDS object. You can build this from STRING DB (e.g. human interactions with experimental score > 700).
2. **`data/chembl_uniprot_mapping.txt`**: Tab-separated target mapping file.
3. **`data/lincs_gene_meta.csv`**: LINCS probe mapping metadata.
4. **`data/test_datasets/`**: A folder containing the disease signatures CSV files matching the datasets you want to evaluate (e.g., `glaucoma.csv` containing LINCS gene z-score profiles, or control/perturbed files for Excape datasets).

---

## Pipeline Execution Steps

The pipeline runs the following sequential steps, which can be run all at once or individually:

1. **Step 1: Merge Outputs**: Combines generative method folder outputs (SmilesGen, GexMolGen, TransGEM, Bicev, Triomphe, etc.), standardizes columns to `SMILES`, and attaches inhibitor references.
2. **Step 2: Calculate Properties**: Computes 2D chemical descriptors (QED, SA Score, MW, LogP, TPSA, CSP3 fraction, rotatable bonds, etc.), PAINS/BRENK/NIH safety filters, reference similarities, group rankings, and set-level summaries (Validity, Uniqueness, Novelty, Internal Diversity).
3. **Step 3: Filter & Generate Conformations**: Filters compounds (Lipinski Rules of 5 + Complexity checks) and generates 3D conformations using RDKit's ETKDGv3.
4. **Step 4: PoseBusters**: Runs PoseBusters validation checks (using config `"mol"`) on the generated SDF conformation file.
5. **Step 5: ChEMBL Lookups**: Queries the ChEMBL API client in parallel to fetch target associations (60% similarity threshold) and consensus potency scores (80% similarity threshold).
6. **Step 6: Disease Reversal**: Standardizes targets, maps SwissProt IDs to HGNC symbols, runs parallel DTSEA enrichment calculations in R, and performs consensus potency scaling.
7. **Step 7: Merge Master Dataframe**: Combines the output results from all steps into a single aligned final master dataframe.
8. **Step 8: Consensus Ranking**: Applies a 15-dimensional Pareto front optimization (using rPref) and Robust Rank Aggregation (RRA) ranking to select the top 10 consensus molecules per method and dataset.

---

## CLI Usage

Run the entire pipeline end-to-end:

```bash
theraflow \
  --method_outputs_dir results/method_outputs \
  --output_dir results/pipeline_runs
```

### Argument Flags

* `--method_outputs_dir`: Directory containing raw method subfolders (default: `../../results/method_outputs`).
* `--output_dir`: Output directory for generated pipeline runs (default: `../../results/pipeline_runs`).
* `--steps`: Command-separated steps to run (default: `1,2,3,4,5,6,7,8`).
* `--ppi_rds`: Path to PPI rds network file (default: `data/ppi.rds`). Note: The default dataset was sourced from STRING DB and filtered to high-confidence experimental interactions (confidence score > 700).
* `--mapping_file`: Path to ChEMBL-UniProt mapping text file (default: `data/chembl_uniprot_mapping.txt`).
* `--lincs_meta`: Path to LINCS gene meta CSV (default: `Omics2Drug/data/lincs_gene_meta.csv`).
* `--test_datasets_dir`: Directory containing disease signature signature files (default: `Omics2Drug/data/test_datasets`).
* `--n_jobs_3d`: Number of parallel jobs for RDKit 3D embedding (default: `4`).
* `--workers_r`: Number of parallel workers for R DTSEA scoring (default: `4`).
* `--workers_chembl`: Number of threads for ChEMBL API client requests (default: `5`).

---

## Python API Usage

You can also import individual steps of the pipeline in your Python code:

```python
from pathlib import Path
import os

# 1. Merge raw outputs
from utils import merge_outputs
merge_outputs(
    method_outputs_dir="results/method_outputs",
    output_csv_path="results/pipeline_runs/step1_merged.csv"
)

# 2. Run molecular properties calculations
from molecular_properties import calculate_properties_pipeline
calculate_properties_pipeline(
    input_csv_path="results/pipeline_runs/step1_merged.csv",
    output_csv_path="results/pipeline_runs/step2_properties.csv",
    output_summary_csv_path="results/pipeline_runs/step2_summary.csv"
)
```
