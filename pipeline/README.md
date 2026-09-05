# TheraFlow Downstream Evaluation & Prioritization Pipeline

The `TheraFlow/pipeline/` package is an end-to-end, modular benchmarking, filtering, target profiling, disease reversal scoring, multi-objective consensus ranking, and structural simulation framework for generative drug discovery models evaluated on transcriptomic signatures.

---

## 📋 Table of Contents
1. [Pipeline Architecture](#-pipeline-architecture)
2. [Prerequisites & Environment Setup](#-prerequisites--environment-setup)
3. [Master Orchestrator (One-Click Execution)](#-master-orchestrator-one-click-execution)
4. [Step-by-Step Execution Guide](#-step-by-step-execution-guide)
   - [Step 1: Generative Performance Macro-Metrics](#step-1-generative-performance-macro-metrics-1generative_perfpy)
   - [Step 2: 2D Physicochemical Properties, Safety Filters & Similarities](#stage-2-2d-physicochemical-properties-safety-filters--similarities-2molecular_propertiespy)
   - [Step 3: InChIKey Deduplication, 3D Embedding & PoseBusters](#stage-3-inchikey-deduplication-3d-embedding--posebusters-3posebuster_convert_sdfpy)
   - [Step 4: ChEMBL Target Profiling, Potency & Mechanism of Action](#stage-4-chembl-target-profiling-potency--mechanism-of-action-4potency_chemblpy)
   - [Step 5: DTSEA Disease DEG Reversal Scoring](#stage-5-dtsea-disease-deg-reversal-scoring-5disease_reversalpy--disease_reversalr)
   - [Step 6: Master Merging, Multi-Objective Pareto Optimization & RRA Consensus Ranking](#stage-6-master-merging-multi-objective-pareto-optimization--rra-consensus-ranking-6mergingpy--consensus_rankingr)
   - [Step 7: AlphaFold 3 Docking, Pocket Interface Extraction & QM/MM Setup](#stage-7-alphafold-3-docking-pocket-interface-extraction--qmmm-setup-7docking_md_setuppy--take_interfacepy)
5. [Outputs & Data Dictionary](#-outputs--data-dictionary)

---

## 🏗️ Pipeline Architecture

#### Figure insert here

---

## ⚙️ Prerequisites & Environment Setup

### 1. Python Environment
Install core Python dependencies via pip:
```bash
pip install -r TheraFlow/requirements.txt
```

Key Python libraries:
- `rdkit` (2D descriptors, fingerprints, ETKDGv3 conformation embedding, MMFF94 minimization)
- `posebusters` (Physical and chemical validity validation)
- `chembl-webresource-client` (Programmatic ChEMBL bioactivity and target API access)
- `pandas`, `numpy`, `scipy`, `joblib`, `tqdm` (Multiprocessing & data manipulation)

### 2. R Environment (for Steps 5 & 6)
Steps 5 (DTSEA) and 6 (Consensus Ranking) call optimized R scripts. Ensure R ($\ge 4.0$) is installed with the following packages:

```R
install.packages(c("dplyr", "readr", "data.table", "tidyr", "rPref", "RobustRankAggreg", "jsonlite", "igraph", "purrr", "stringr", "BiocManager"))
BiocManager::install("DTSEA")
```

### 3. Structural Tools (for Step 7, Optional)
- **AlphaFold 3**: Clone repository and configure inference virtual environment.
- **PDB2PQR** (`pdb2pqr`): Used for PARSE force-field protonation and partial charge calculation at pH 7.4.
- **ORCA** / **xTB**: Quantum chemistry calculation packages.

---

## ⚡ Master Orchestrator (One-Click Execution)

You can execute the entire pipeline (Step 1 through 6, and optionally Step 7) with a single command. The master orchestrator runs all steps sequentially, passes intermediate files between steps, generates consensus candidate rankings, and merges all calculated metrics back into your input candidate dataset with an automatic timestamped backup:

```bash
# Run full end-to-end evaluation & merge annotations back into processed_data.csv
python TheraFlow/run_pipeline.py   --input TheraFlow/method_processing/data/processed_data.csv   --output_dir TheraFlow/outputs   --output_merged TheraFlow/method_processing/data/processed_data.csv
```

### Master CLI Options (`run_pipeline.py` / `pipeline/run_all.py`):
| Option | Type | Default | Description |
|---|---|---|---|
| `--input` | `str` | `TheraFlow/method_processing/data/processed_data.csv` | Input candidate SMILES dataset. |
| `--output_dir` | `str` | `TheraFlow/outputs` | Directory to store all intermediate and final outputs. |
| `--output_merged` | `str` | `TheraFlow/method_processing/data/processed_data.csv` | Path to output dataset merged with all computed annotations. |
| `--backup` | `bool` | `True` | Automatically create a `.backup.csv` file before overwriting target dataset. |
| `--steps` | `str` | `1,2,3,4,5,6` | Comma-separated list of step numbers to run (e.g., `--steps 1,2,3`). |
| `--n_jobs_3d` | `int` | `4` | Parallel worker processes for 3D ETKDGv3 embedding in step 3. |
| `--chembl_workers` | `int` | `5` | Concurrent threads for ChEMBL queries in step 4. |
| `--r_workers` | `int` | `2` | Parallel processes for DTSEA calculations in step 5. |
| `--top_n_consensus` | `int` | `10` | Number of top consensus candidates selected per method/dataset in step 6. |
| `--run_af3` | `flag` | `False` | Launch AlphaFold 3 batch docking inference in step 7. |
| `--extract_qm` | `flag` | `False` | Extract AF3 docking and QM metrics in step 7. |

---

## 🚀 Step-by-Step Execution Guide

You can also run each pipeline step independently as a modular command-line tool.

---

### Step 1: Generative Performance Macro-Metrics (`1.generative_perf.py`)
Computes distribution-level quality and diversity metrics grouped by generative `method` and `test_dataset`.

- **Metrics Computed**:
  - **Validity**: Proportion of SMILES strings parseable by RDKit into valid chemical graphs.
  - **Uniqueness**: Proportion of unique canonical SMILES among valid molecules.
  - **Novelty (MOSES)**: Fraction of valid molecules not present in the MOSES reference chemical library (automatically downloaded and cached).
  - **Internal Diversity (IntDiv)**: $1 - rac{1}{|S|^2} \sum_{m_1, m_2 \in S} T(m_1, m_2)$ using Morgan fingerprints (radius 2, 2048 bits).

```bash
python TheraFlow/pipeline/1.generative_perf.py   --input TheraFlow/method_processing/data/processed_data.csv   --output TheraFlow/outputs/set_properties.csv   --ref_cache TheraFlow/outputs/ref_cache
```

---

### Step 2: 2D Physicochemical Properties, Safety Filters & Similarities (`2.molecular_properties.py`)
Evaluates single-molecule 2D descriptors, drug-likeness rules, structural alerts, and reference similarity.

- **Metrics Computed**:
  - **Physicochemical Properties**: Molecular Weight (MW), $\log P$ (MolLogP), Topological Polar Surface Area (TPSA), Quantitative Estimate of Drug-likeness (QED), Synthetic Accessibility Score (SA Score), Fraction $C_{sp3}$ ($F_{Csp3}$), H-Bond Acceptors (HBA), H-Bond Donors (HBD), Rotatable Bonds, Ring Count, Aromatic Rings.
  - **Safety & Structural Alerts**: PAINS (Pan-Assay Interference Compounds), BRENK (unwanted reactive groups), and NIH filter alerts.
  - **Reference Similarities**: Max and mean Tanimoto similarity against known benchmark reference inhibitors/perturbations for each test dataset.

```bash
python TheraFlow/pipeline/2.molecular_properties.py   --input TheraFlow/method_processing/data/processed_data.csv   --output TheraFlow/outputs/molecular_properties_metrics.csv   --output_filtered TheraFlow/outputs/molecular_properties_filtered.csv
```

---

### Step 3: InChIKey Deduplication, 3D Embedding & PoseBusters (`3.posebuster_convert_sdf.py`)
Performs canonical InChIKey deduplication, Lipinski & complexity filtering, generates 3D geometries, and validates physical plausibility.

- **Methodology**:
  - Deduplicates molecules via canonical InChIKey.
  - Filters molecules meeting basic drug-likeness criteria ($	ext{MW} \le 800$, $	ext{Rotatable Bonds} \le 15$).
  - Generates 3D coordinates using RDKit **ETKDGv3** followed by **MMFF94** force field energy minimization.
  - Saves optimized 3D structures to SDF format.
  - Executes **PoseBusters** validation suite to verify chemical consistency, bond lengths, angles, steric clashes, and planar aromatic rings.

```bash
python TheraFlow/pipeline/3.posebuster_convert_sdf.py   --input TheraFlow/outputs/molecular_properties_metrics.csv   --output_sdf TheraFlow/outputs/filtered_3d.sdf   --output_csv TheraFlow/outputs/molecular_properties_with_posebusters.csv   --n_jobs 4
```

---

### Step 4: ChEMBL Target Profiling, Potency & Mechanism of Action (`4.potency_chembl.py`)
Queries the ChEMBL database to retrieve structural analogs, predicted protein targets, bioactivities, and Mechanism of Action (MoA).

- **Methodology**:
  - Performs similarity searches against ChEMBL at two similarity cutoffs: $\ge 60\%$ (moderate similarity) and $\ge 80\%$ (high similarity).
  - Retrieves associated human target proteins (UniProt IDs, Gene Symbols) and experimental bioactivity measurements ($	ext{IC}_{50}$, $	ext{EC}_{50}$, $K_i$, $K_d$, $	ext{pChEMBL}$).
  - Extracts annotated Mechanism of Action (MoA) and pharmacological action types (inhibitor, agonist, antagonist).

```bash
python TheraFlow/pipeline/4.potency_chembl.py   --input TheraFlow/outputs/molecular_properties_filtered.csv   --output_targets TheraFlow/outputs/chembl_targets.csv   --output_potency TheraFlow/outputs/potency_scores.csv   --max_workers 5
```

---

### Step 5: DTSEA Disease DEG Reversal Scoring (`5.disease_reversal.py` / `disease_reversal.R`)
Executes Drug-Target Set Enrichment Analysis (DTSEA) in R to quantify each candidate compound's ability to reverse disease-specific differential expression signatures.

- **Methodology**:
  - Uses the human Protein-Protein Interaction (PPI) network (`full_human_ppi.rds`).
  - Weights drug target perturbations based on ChEMBL bioactivity ($-\log_{10}(	ext{Potency})$) and drug action mode (+1 for agonist/activator, -1 for inhibitor/antagonist).
  - Computes DTSEA Enrichment Score (ES), Normalized Enrichment Score (NES), p-value, and FDR q-value against disease DEG profiles.

```bash
python TheraFlow/pipeline/5.disease_reversal.py   --targets_csv TheraFlow/outputs/chembl_targets.csv   --potency_csv TheraFlow/outputs/potency_scores.csv   --output TheraFlow/outputs/reversal_res.csv   --ppi_rds TheraFlow/method_processing/data/full_human_ppi.rds   --lincs_meta TheraFlow/method_processing/data/LINCS/lincs_gene_meta.csv   --test_datasets_dir TheraFlow/method_processing/data/test_datasets/rest_methods   --workers 2
```

---

### Step 6: Master Merging, Multi-Objective Pareto Optimization & RRA Consensus Ranking (`6.merging.py` / `consensus_ranking.R`)
Consolidates all computed metric layers into a master table, performs Pareto frontier optimization across 11 property objectives, and applies Robust Rank Aggregation (RRA).

- **Methodology**:
  - Merges 2D properties, PoseBusters validity, ChEMBL targets/potency, and DTSEA reversal scores into `merged_master_dataframe.csv`.
  - Performs 11-dimensional Pareto optimization via `rPref` (maximizing QED, Potency, Reversal NES, PoseBusters pass rate, and minimizing SA Score, Molecular Weight, and Alert penalties).
  - Aggregates rank lists across criteria using `RobustRankAggreg` (RRA) to produce consensus prioritization scores.
  - Selects top $N$ candidates per method and dataset (`consensus_top.csv`).
  - Automatically exports top candidates and target sequences into AlphaFold 3 batch input JSON format (`af3_batch_input.json`).

```bash
python TheraFlow/pipeline/6.merging.py   --properties_csv TheraFlow/outputs/molecular_properties_metrics.csv   --posebusters_csv TheraFlow/outputs/molecular_properties_with_posebusters.csv   --targets_csv TheraFlow/outputs/chembl_targets.csv   --potency_csv TheraFlow/outputs/potency_scores.csv   --reversal_csv TheraFlow/outputs/reversal_res.csv   --master_csv TheraFlow/outputs/merged_master_dataframe.csv   --output_consensus TheraFlow/outputs/consensus_top.csv   --output_af3_json TheraFlow/outputs/af3_batch_input.json   --top_n 10
```

---

### Step 7: AlphaFold 3 Docking, Pocket Interface Extraction & QM/MM Setup (`7.docking_md_setup.py` / `take_interface.py`)
Performs 3D structural validation of top consensus candidate-target complexes via AlphaFold 3 batch docking, automated active site pocket interface extraction, and QM/MM input preparation.

#### 1. Batch AlphaFold 3 Inference (Multi-GPU Parallelization)
```bash
python TheraFlow/pipeline/7.docking_md_setup.py   --run_af3   --af3_json TheraFlow/outputs/af3_batch_input.json   --af3_dir alphafold3   --af3_output TheraFlow/outputs/af3_predictions   --gpus 0 1 2 3
```

#### 2. Pocket Interface Extraction & ORCA / xTB QM Setup
Converts AF3 CIF predictions to PDB, runs `pdb2pqr` (PARSE forcefield, pH 7.4) to determine protonation states and partial charges, extracts residue clusters within 4.0 Å of the ligand using `take_interface.py`, and generates ORCA (DFT B3LYP/def2-SVP) and GFN2-xTB inputs:
```bash
python TheraFlow/pipeline/7.docking_md_setup.py   --setup_qm   --cif_dir TheraFlow/outputs/af3_predictions   --qm_output TheraFlow/outputs/docking_md_runs   --cutoff 4.0   --nproc 16
```

#### 3. Extract Docking Confidence & QM Descriptors
Parses AF3 prediction summaries and QM calculation logs into a consolidated metrics table:
```bash
python TheraFlow/pipeline/7.docking_md_setup.py   --extract_metrics   --af3_output TheraFlow/outputs/af3_predictions   --output_metrics TheraFlow/outputs/docking_qm_metrics.csv
```

---

## 📊 Outputs & Data Dictionary

| Output File | Primary Source Step | Key Columns & Annotations |
|---|---|---|
| `set_properties.csv` | Step 1 | `method`, `test_dataset`, `validity`, `uniqueness`, `novelty_moses`, `internal_diversity` |
| `molecular_properties_metrics.csv` | Step 2 | `SMILES`, `MW`, `LogP`, `TPSA`, `QED`, `SA_Score`, `F_Csp3`, `HBA`, `HBD`, `pains_alert`, `brenk_alert`, `nih_alert`, `sim_mean_ref`, `sim_max_ref` |
| `molecular_properties_filtered.csv` | Step 2 | Subset passing Lipinski Rule of 5 and zero severe structural alerts. |
| `filtered_3d.sdf` | Step 3 | 3D coordinate conformers optimized with MMFF94 force field in SDF format. |
| `molecular_properties_with_posebusters.csv` | Step 3 | `inchikey`, `posebusters_passed`, `pb_passes_valence_check`, `pb_passes_geometry_check`, `pb_passes_clash_check` |
| `chembl_targets.csv` | Step 4 | `SMILES`, `chembl_id`, `target_chembl_id`, `target_gene_symbol`, `uniprot_id`, `moa_action_type`, `similarity` |
| `potency_scores.csv` | Step 4 | `SMILES`, `best_potency_nm`, `pChEMBL_mean`, `pChEMBL_max`, `has_nanomolar_activity`, `num_targets` |
| `reversal_res.csv` | Step 5 | `SMILES`, `test_dataset`, `dtsea_es`, `dtsea_nes`, `dtsea_pvalue`, `dtsea_padj`, `reversal_score` |
| `merged_master_dataframe.csv` | Step 6 | Complete unified dataset joining all 2D properties, 3D flags, ChEMBL potency, and DTSEA reversal metrics. |
| `consensus_top.csv` | Step 6 | Top ranked candidate molecules per method and dataset prioritized by Pareto frontier depth and RRA score. |
| `af3_batch_input.json` | Step 6 | Formatted JSON specifications for batch AlphaFold 3 complex structure prediction. |
| `docking_qm_metrics.csv` | Step 7 | `ranking_score`, `ligand_plddt`, `ligand_rmsd`, `homo_lumo_gap_ev`, `dipole_moment_debye`, `scf_energy_hartree` |
