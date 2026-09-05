# TheraFlow: A Multi-Stage Framework for Evaluating Generative Drug Models from Transcriptomics and Prioritizing Candidates

**TheraFlow** is an end-to-end, modular, and parallelized benchmarking and prioritization framework designed to evaluate deep generative models for transcriptomics-driven drug discovery. It bridges the gap between de novo molecular generation from disease differential gene expression signatures and downstream clinical / biophysical prioritization.

---

## 🌟 Key Highlights

- 🧪 **10 Deep Generative Methods**: Standardized execution and reimplementation pipelines for SmilesGEN, BiAEE, BiCEV, PaccMann RL, GxRNN, TransGEM, Triomphe, FAME, GEx2SGen, and GexMolGen.
- 🧬 **7 Benchmark Transcriptomic Datasets**: Automated NCBI GEO download, differential expression computation (limma/edgeR), and landmark gene alignment across LINCS L1000 Phase 1/2 and 4 disease states (Cardiomyopathy, Parkinson's Disease, Hepatocellular Carcinoma, and Glaucoma).
- 🔬 **7-Stage Downstream Evaluation & Prioritization**:
  1. **Generative Performance**: Macro-metrics (Validity, Uniqueness, MOSES Novelty, and Internal Diversity).
  2. **2D Physicochemical & Safety**: Lipinski properties, QED, SA Score, and structural alerts (PAINS, BRENK, NIH).
  3. **3D Conformation & Physical Validity**: ETKDGv3 embedding, MMFF94 force field optimization, and PoseBusters validity checks.
  4. **Bioactivity & Target Profiling**: ChEMBL analog searches ($\ge 60\%$ and $\ge 80\%$), target prediction, IC50/EC50 potency retrieval, and Mechanism of Action (MoA).
  5. **DTSEA Disease Reversal**: Drug-Target Set Enrichment Analysis (DTSEA) in R scaled by compound potency and mode of action.
  6. **11D Pareto Optimization & Consensus Ranking**: Multi-objective Pareto frontier identification and Robust Rank Aggregation (RRA).
  7. **Structural Docking & QM/MM Simulation**: Multi-GPU AlphaFold 3 batch complex prediction, pocket interface cluster extraction, and ORCA/xTB quantum chemistry setup.
- ⚡ **One-Click Master Orchestrator**: `run_pipeline.py` seamlessly executes all evaluation stages, passes data between layers, and merges calculated annotations directly back into the candidate dataset with automatic backup.

---

## 📂 Repository Organization

```
TheraFlow/
├── pyproject.toml              # Modern Python packaging configuration (pip install -e .)
├── requirements.txt            # Unified Python dependency specification
├── README.md                   # Top-level framework documentation (this file)
├── run_pipeline.py             # Master end-to-end evaluation orchestrator
│
├── method_processing/          # Molecule Generation & Preprocessing Sub-package
│   ├── README.md               # Method execution & dataset preprocessing documentation
│   ├── dataset_processing.py   # Unified NCBI GEO download & DEG signature processing
│   ├── dataset_processing.R    # Differential expression analysis in R (limma/edgeR)
│   ├── data/                   # Reference databases and preprocessed input signatures
│   │   ├── test_datasets/      # Landmark-aligned profiles formatted for each model
│   │   └── LINCS/              # LINCS GSE70138 and GSE92742 signature caches
│   ├── outputs/                # Generated raw SMILES libraries for all methods & datasets
│   ├── smilesgen.py            # SmilesGEN orchestrator runner
│   ├── biaee.py                # BiAEE orchestrator runner
│   ├── bicev.py                # BiCEV orchestrator runner
│   ├── paccmann_rl.py          # PaccMann RL orchestrator runner
│   ├── gxrnn.py                # GxRNN orchestrator runner
│   ├── transgem.py             # TransGEM orchestrator runner
│   ├── triomphe.py             # Triomphe orchestrator runner
│   ├── fame.py                 # FAME orchestrator runner (Self-contained)
│   ├── gex2sgen.py             # GEx2SGen orchestrator runner (Self-contained)
│   └── combine.py              # Consolidates all generated outputs into processed_data.csv
│
├── pipeline/                   # Downstream Evaluation, Ranking & Docking Package
│   ├── README.md               # Detailed evaluation pipeline documentation
│   ├── 1.generative_perf.py    # Stage 1: Macro-level generative quality metrics
│   ├── 2.molecular_properties.py # Stage 2: 2D descriptors, safety alerts & similarities
│   ├── 3.posebuster_convert_sdf.py # Stage 3: 3D conformation embedding & PoseBusters
│   ├── 4.potency_chembl.py     # Stage 4: ChEMBL analog searches & potency retrieval
│   ├── 5.disease_reversal.py   # Stage 5: Python wrapper for DTSEA calculations
│   ├── disease_reversal.R      # Stage 5: DTSEA disease reversal calculations in R
│   ├── 6.merging.py            # Stage 6: Python wrapper for consensus ranking & AF3 prep
│   ├── consensus_ranking.R     # Stage 6: 11D Pareto optimization & RRA in R
│   ├── 7.docking_md_setup.py   # Stage 7: AlphaFold 3 batch docking & ORCA QM/MM setup
│   ├── take_interface.py       # Stage 7: Binding site cluster extractor & PDB2PQR parser
│   └── run_all.py              # Master pipeline wrapper (equivalent to run_pipeline.py)
│
└── outputs/                    # Output directory for all evaluation metrics and rankings
    ├── set_properties.csv
    ├── molecular_properties_metrics.csv
    ├── filtered_3d.sdf
    ├── molecular_properties_with_posebusters.csv
    ├── chembl_targets.csv
    ├── potency_scores.csv
    ├── reversal_res.csv
    ├── merged_master_dataframe.csv
    ├── consensus_top.csv
    ├── af3_batch_input.json
    └── docking_qm_metrics.csv
```

---

## 🛠️ Installation & Setup

### 1. Python Environment Setup
Install the Python dependencies in a clean virtual or Conda environment:
```bash
# Clone the repository
git clone [https://github.com/<your-org>/TheraFlow.git](https://github.com/SydneyBioX/TheraFlow.git)
cd TheraFlow

# Install Python dependencies
pip install -r requirements.txt

# (Optional) Install the package in editable mode
pip install -e .
```

### 2. R Dependencies (for DTSEA & Consensus Ranking)
Ensure R ($\ge 4.0$) is installed and install the required CRAN and Bioconductor packages:
```R
install.packages(c("dplyr", "readr", "data.table", "tidyr", "rPref", "RobustRankAggreg", "jsonlite", "igraph", "purrr", "stringr", "BiocManager"))
BiocManager::install("DTSEA")
```

### 3. External Structural Tools (Optional, for Stage 7)
- **AlphaFold 3**: For batch structural complex modeling.
- **PDB2PQR** (`pdb2pqr`): For PARSE protonation and partial charge assignment.
- **ORCA / xTB**: For quantum mechanics optimization and electronic property calculations.

---

## ⚡ Quick Start

### Step 1: Preprocessing & Candidate Generation
To download raw signatures, compute disease DEGs, align input profiles, and run candidate generation across models:
```bash
# 1. Download raw GEO datasets and preprocess signatures
python TheraFlow/method_processing/dataset_processing.py

# 2. Run generative models (e.g. SmilesGEN, BiAEE, etc.)
python TheraFlow/method_processing/smilesgen.py --run-all
python TheraFlow/method_processing/biaee.py --run-all

# 3. Consolidate generated SMILES into master dataset
python TheraFlow/method_processing/combine.py
```
*See [method_processing/README.md](file:///users/stgrad/rojashrj/Projects/omics2mol/TheraFlow/method_processing/README.md) for full details on method environments and runner scripts.*

### Step 2: Master Downstream Evaluation & Ranking
Run the complete multi-stage evaluation pipeline end-to-end with a single command:
```bash
python TheraFlow/run_pipeline.py   --input TheraFlow/method_processing/data/processed_data.csv   --output_dir TheraFlow/outputs   --output_merged TheraFlow/method_processing/data/processed_data.csv
```
This executes Stages 1 through 6, produces top consensus candidates (`consensus_top.csv`), exports AlphaFold 3 batch input JSON (`af3_batch_input.json`), and automatically merges all calculated properties and rank flags back into `processed_data.csv` (with automatic backup).

*See [pipeline/README.md](file:///users/stgrad/rojashrj/Projects/omics2mol/TheraFlow/pipeline/README.md) for individual stage commands, CLI flags, and data schemas.*

---

## 🔬 Benchmark Generative Models & Datasets

### Generative Models Supported:
| Model | Representation | Architecture / Methodology | Runner Script |
|---|---|---|---|
| **SmilesGEN** | SMILES | Reinforcement Learning & RNN | `smilesgen.py` |
| **BiAEE** | SMILES | Bidirectional Adversarial Autoencoder | `biaee.py` |
| **BiCEV** | SMILES | Conditional Variational Autoencoder | `bicev.py` |
| **PaccMann RL** | SMILES | Multimodal RL for Anticancer Compounds | `paccmann_rl.py` |
| **GxRNN** | SMILES | Expression-conditioned RNN | `gxrnn.py` |
| **TransGEM** | SMILES | Transformer-based Expression-to-Molecule | `transgem.py` |
| **Triomphe** | SELFIES | Transcriptome-conditioned Deep Generator | `triomphe.py` |
| **FAME** | SELFIES | Feature-Augmented Molecular Exploration (Reimplemented) | `fame.py` |
| **GEx2SGen** | SELFIES | Gene Expression to SELFIES Generator (Reimplemented) | `gex2sgen.py` |
| **GexMolGen** | SMILES / Web | Conditional Generative Framework | `gexmolgen` |

### Benchmark Transcriptomic Datasets:
1. `lincs1000_1`: LINCS L1000 Phase 1 landmark perturbation signatures.
2. `lincs1000_2`: LINCS L1000 Phase 2 landmark perturbation signatures.
3. `cardiomyopathy`: Dilated cardiomyopathy disease differential expression profile.
4. `parkinsons`: Parkinson's disease substantia nigra differential expression profile.
5. `hcc`: Hepatocellular carcinoma tumor vs. normal differential expression profile.
6. `glaucoma`: Primary open-angle glaucoma trabecular meshwork profile.
7. `glaucoma_bulk`: Glaucoma bulk RNA-seq differential signature.

---

## 📊 Summary of Evaluation Outputs

All pipeline outputs are structured under `TheraFlow/outputs/`:

| Output File | Description |
|---|---|
| `set_properties.csv` | Distributional quality metrics (Validity, Uniqueness, MOSES Novelty, Internal Diversity). |
| `molecular_properties_metrics.csv` | 2D Physicochemical descriptors, SA score, QED, PAINS/BRENK/NIH structural alerts. |
| `filtered_3d.sdf` | Energy-minimized 3D molecular conformations generated by ETKDGv3. |
| `molecular_properties_with_posebusters.csv` | InChIKey-deduplicated table with PoseBusters physical validity flags. |
| `chembl_targets.csv` | Matched ChEMBL analog compounds, target genes, UniProt IDs, and Mechanism of Action. |
| `potency_scores.csv` | Retrieved experimental bioactivity metrics (IC50, EC50, pChEMBL, nM potency). |
| `reversal_res.csv` | DTSEA disease DEG reversal scores (Enrichment Score, NES, p-value, FDR q-value). |
| `merged_master_dataframe.csv` | Unified master table combining all 2D, 3D, ChEMBL potency, and reversal metrics. |
| `consensus_top.csv` | Top prioritized candidates per method and dataset from 11D Pareto & RRA ranking. |
| `af3_batch_input.json` | Ready-to-run batch input specification for AlphaFold 3 complex structure prediction. |
| `docking_qm_metrics.csv` | AlphaFold 3 ligand pLDDT, complex ranking score, and ORCA/xTB QM electronic descriptors. |
