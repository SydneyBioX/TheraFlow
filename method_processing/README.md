# Method Processing: Standalone Molecule Generation Orchestration
This folder contains the orchestration scripts for generating hit-like molecules using **10 different cross-modal methods** across **7 distinct gene expression datasets**.

---

## 📂 Folder Structure & Prerequisites
To run the models standalone, clone their repositories and structure the parent directory of this folder (`method_processing/`) as follows:

```
omics2mol/                  <-- Main Project Directory
│
├── method_processing/      <-- THIS FOLDER (Orchestration & Templates)
│   ├── README.md           <-- Documentation
│   ├── raw_data_access.py  <-- Programmatic download of GEO LINCS datasets
│   ├── dataset_processing.py   <-- Unified Preprocessing
│   ├── run_all_generation.py   <-- Master Runner
│   └── run_gexmolgen_cli.py    <-- GExMolGen Offline CLI
│
├── data/
│   └── test_datasets/      <-- Folder containing dataset .csv files
│       ├── smilesgen/      <-- VAE / FAME / GEx2sGen / TransGEM inputs
│       ├── paccmann_rl/    <-- Paccmann inputs
│       └── gexmolgen/      <-- GExMolGen control/disease pairs
│
└── [Model Repositories]    <-- Cloned model folders (see cloning instructions below)
    ├── SmilesGEN/
    ├── TransGEM/
    ├── Triomphe/
    └── ...
```

---

## 📥 Cloning Instructions
Run the following commands in the parent folder to clone all required repositories:

```bash
# 1. SmilesGEN
git clone https://github.com/hliulab/SmilesGEN.git

# 2. TransGEM
git clone https://github.com/hzauzqy/TransGEM.git

# 3. paccmann_rl
git clone https://github.com/PaccMann/paccmann_rl.git

# 4. BiAAE
git clone https://github.com/insilicomedicine/BiAAE.git

# 5. BiCEV
git clone https://github.com/ChemEXL/BiCEV.git

# 6. GExMolGen
git clone https://github.com/Bunnybeibei/GExMolGen.git
```

*Note: Models such as Triomphe, GxRNN, and the FAME/GEx2sGen reimplementations should be extracted or cloned into the main directory as well.*

---


## 📥 Accessing Raw LINCS/GEO Datasets
Because the raw supplementary GCTX files are large (~5.6 GB for Phase 2 Level 5 COMPZ and ~53 GB for Phase 1 Level 3 INF), they are not versioned in Git. You can programmatically download and unpack them to the correct model folders using:

```bash
# Print download URLs and targets without downloading
python raw_data_access.py --dry-run

# Download & extract GSE70138 dataset (used by SmilesGEN, paccmann_rl, etc. ~5.6 GB)
python raw_data_access.py --gse70138
```

## 🛠️ Data Preprocessing
Each model expects its datasets in a specific path inside `data/test_datasets/`.
You can preprocess and structure all raw dataset signatures into their required formats with the unified script **`dataset_processing.py`**:

To preprocess all datasets at once:
```bash
python dataset_processing.py --run-all
```

Or preprocess for specific models:
```bash
python dataset_processing.py --smilesgen --transgem
```

---

## 🚀 Running the Master Orchestrator

The main runner script is **`run_all_generation.py`**. It automatically maps the correct Python interpreters (with pre-configured conda environments), prepares dataset inputs, dynamically instantiates templates, executes generation inside the corresponding model directories, and logs outcomes.

### Execution Command:
To run a **specific method** on a **specific dataset**:
```bash
python run_all_generation.py --method <method_name> --dataset <dataset_name>
```

To run **all combinations** sequentially:
```bash
python run_all_generation.py --run-all
```

### Options:
*   **`--method`**: `smilesgen_vae_gen`, `gexmolgen`, `triomphe`, `transgem`, `paccmann_rl`, `gxrnn`, `biaee`, `bicev`, `fame`, `gex2sgen`
*   **`--dataset`**: `Lincs1000_1`, `Lincs1000_2`, `cardiomyopathy`, `glaucoma`, `glaucoma_bulk`, `hcc`, `parkinsons`

---


---

> [!IMPORTANT]
> **Workspace Structure and Cloned Repositories**:
> The orchestration runner scripts are configured using relative parent paths (e.g. `../SmilesGEN/`, `../TransGEM/`, etc.). This requires that the official GitHub repositories for the respective methods were cloned directly into the parent directory.
