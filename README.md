# TheraFlow: A Multi-Stage Framework for evaluating generative drug models from transcriptomics and prioritizing candidates


TheraFlow is a comprehensive, modular, and parallelized replication package for:
1. **Molecule Preprocessing and Generation (`method_processing/`)**: Orchestrates and executes 10 generative methods on 7 transcriptomic datasets.
2. **Metrics Calculation and Benchmarking (`pipeline/`)**: Evaluates generated molecules across properties, safety filters, conformation validation (PoseBusters), activity/potency queries (ChEMBL), disease reversal scoring (DTSEA), and multi-dimensional consensus ranking.

---

## 📂 Project Repository Layout

This repository is organized as follows:
```
TheraFlow/
├── pyproject.toml              # Modern Python packaging configuration
├── README.md                   # Repository root documentation (this file)
├── .gitignore                  # Standard Git ignore configurations
│
├── method_processing/          # Molecule Generation & Preprocessing Sub-package
│   ├── README.md               # Preprocessing and model generation documentation
│   ├── raw_data_access.py      # Programmatic download of LINCS raw GEO datasets
│   ├── dataset_processing.py   # Preprocesses transcriptomic signatures
│   ├── run_all_generation.py   # Master orchestrator script to run all 10 methods
│
└── pipeline/                   # Downstream Benchmarking & Metrics Package (Python/R)
    ├── README.md               # Downstream evaluation documentation
    ├── utils.py                # 3D SDF conformation generators
    ├── molecular_properties.py # Descriptors and safety filters
    ├── posebuster_analysis.py  # PoseBusters validation checks
    ├── potency_chembl.py       # ChEMBL activity fetchers
    ├── disease_reversal.py     # R DTSEA wrapper
    ├── disease_reversal.R      # R disease reversal calculations script
    ├── consensus_ranking.py    # Python wrapper for ranking script
    ├── consensus_ranking.R     # R consensus ranking calculations script
    └── main_pipeline.py        # Master pipeline orchestrator (CLI entrypoint)
```

---

## 🚀 Getting Started

### 1. Installation
Install the TheraFlow package in editable mode from the repository root:
```bash
pip install -e .
```
Once installed, the CLI command `theraflow` is accessible directly from your shell.

### 2. Run Preprocessing and Generation
Prepare target profiles and execute generative models:
```bash
cd method_processing
# Download raw datasets (fetches ~5.6 GB Phase 2 data and ~53 GB Phase 1 data)
python raw_data_access.py --gse70138
# Preprocess signatures
python dataset_processing.py --run-all
# Generate compounds
python run_all_generation.py --run-all
```
For more information, see the [method_processing/README.md](file:///users/stgrad/rojashrj/Projects/omics2mol/TheraFlow/method_processing/README.md).

### 3. Run Downstream Evaluation and Benchmarking
Analyze generated molecule libraries, validate conformations, query potency, perform disease reversal scoring, and rank results:
```bash
theraflow --method_outputs_dir <path_to_outputs> --output_dir <path_to_results>
```
For more information, see the [pipeline/README.md](file:///users/stgrad/rojashrj/Projects/omics2mol/TheraFlow/pipeline/README.md).

> [!IMPORTANT]
> **Note on parent directory dependency**: Most of the method processing is configured using relative paths (e.g. `../SmilesGEN`, `../Triomphe`, etc.) pointing to the parent folder. This is because the respective methods were cloned directly from their official GitHub repositories into the main workspace.
