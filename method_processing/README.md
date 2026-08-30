# 🧬 TheraFlow Method Processing

This directory serves as the centralized orchestration and execution hub for preparing test datasets and running all 9 deep generative models for drug candidate generation from transcriptomic signatures.

---

## 📊 1. Dataset Preprocessing

The preprocessing pipeline downloads raw datasets, computes disease differential gene expression signatures, and aligns Entrez IDs or Gene Symbols to the landmark shapes expected by each target model.

To run the unified download and preprocessing script:
```bash
python TheraFlow/method_processing/dataset_processing.py
```

### ⚙️ Pipeline Steps:
* 📥 **NCBI GEO Retrieval**: Downloads raw LINCS Level 5 signatures and metadata (GSE70138 & GSE92742).
* 🧪 **Differential Expression**: Invokes internal R modules (`edgeR`/`limma`) on GSE raw count data to derive z-scored log-fold changes for disease states (Cardiomyopathy, Parkinson's, HCC, Glaucoma).
* 🔌 **Control Extraction**: Filters and extracts disease-specific drug perturbation control signatures.
* 💾 **Unified Alignment**: Outputs aligned landmark gene datasets directly to `data/test_datasets/`.

---

## 🚀 2. Running the 9 Generative Methods

The corresponding GitHub repository or software of each method was accessed into the parent directory of `TheraFlow` and the required package was installed in a clean Conda environment.

| Method | Conda Environment | Runner Script |
|---|---|---|
| **SmilesGEN** | `smilesgen` | `smilesgen.py` |
| **BiAEE** | `biaee` | `biaee.py` |
| **BiCEV** | `bicev` | `bicev.py` |
| **PaccMann RL** | `paccmann_rl` | `paccmann_rl.py` |
| **GxRNN** | `gxrnn` | `gxrnn.py` |
| **TransGEM** | `transgem` | `transgem.py` |
| **Triomphe** | `triomphe` | `triomphe.py` |
| **FAME** | `fame` | `fame.py` |
| **GEx2SGen** | `gex2sgen` | `gex2sgen.py` |
| **GexMolGen** | `gexmolgen` | `PARENTDIR/gexmolgen/server.py` |

---

### 1. SmilesGEN (Molecular Optimization)
* **Execution**:
  ```bash
  git clone https://github.com/hliulab/SmilesGEN.git
  python method_processing/smilesgen.py --run-all
  ```

### 2. BiAEE
* **Execution**:
  ```bash
  git clone https://github.com/insilicomedicine/BiAAE.git
  conda activate biaee
  python method_processing/biaee.py --run-all
  ```

### 3. BiCEV
* **Execution**:
  ```bash
  git clone https://github.com/ChemEXL/BiCEV.git
  conda activate bicev
  python method_processing/bicev.py --run-all
  ```

### 4. PaccMann RL
* **Execution**:
  ```bash
  git clone https://github.com/PaccMann/paccmann_rl.git
  conda activate paccmann_rl
  python method_processing/paccmann_rl.py --run-all
  ```

### 5. GxRNN
* **Source**: The source code was accessed from [Nagoya CS GxRNN Page](https://yamanishi.cs.i.nagoya-u.ac.jp/gxrnn/)
* **Execution**:
  ```bash
  conda activate gxrnn
  python method_processing/gxrnn.py --run-all
  ```

### 6. TransGEM
* **Execution**:
  ```bash
  git clone https://github.com/hzauzqy/TransGEM.git
  conda activate transgem
  python method_processing/transgem.py --run-all
  ```

### 7. Triomphe
* **Source**: The source code was accessed from [Kyutech Triomphe Page](http://labo.bio.kyutech.ac.jp/~yamani/triomphe/)
* **Execution**:
  ```bash
  conda activate triomphe
  python method_processing/triomphe.py --run-all
  ```

### 8. FAME
* **Source**: This code is reimplemented (process details given in `/reimplementation/fame` directory)
* **Execution**:
  ```bash
  conda activate fame
  python method_processing/fame.py --run-all
  ```

### 9. GEx2SGen
* **Source**: This code is reimplemented (process details given in `/reimplementation/gex2sgen` directory)
* **Execution**:
  ```bash
  conda activate gex2sgen
  python method_processing/gex2sgen.py --run-all

### 10. GexMolGen
* **Execution**:
  ```bash
  git clone https://github.com/Bunnybeibei/GexMolGen.git
  conda activate gexmolgen
  python server.py
  ```
  Upload the required data from the data/test_datasets/gexmolgen/ into the Standard tab of the webserver

---

## 🛠️ 3. CLI Arguments (Standardized for All Runners)

Each python runner script supports the following command-line interface:
* `--run-all`: Run the model sequentially on all 7 benchmark datasets.
* `--dataset <name>`: Run only on a specific benchmark dataset. Available choices:
  `lincs1000_1`, `lincs1000_2`, `cardiomyopathy`, `parkinsons`, `hcc`, `glaucoma`, `glaucoma_bulk`.

---

## 📁 4. Directory Structure (once the processing and methods are run)

```
TheraFlow/method_processing/
├── data/                         # Preprocessed datasets and reference metadata
│   ├── test_datasets/            # Subfolders containing aligned signatures for each model
│   │   ├── gexmolgen/            # GexMolGen input profiles (for FAME & GEx2SGen)
│   │   ├── transgem/             # TransGEM input profiles
│   │   ├── triomphe/             # Triomphe input profiles
│   │   ├── paccmann_rl/          # PaccMann RL input profiles
│   │   └── rest_methods/         # Inputs for SmilesGEN, BiAEE, BiCEV, GxRNN, FAME and GEx2SGen
│   └── LINCS/                    # Local LINCS reference databases
│       ├── GSE70138/             # GSE70138 compz profiles
│       └── GSE92742/             # GSE92742 compz profiles
├── outputs/                      # Output generated SMILES files for each method & dataset
│   ├── smilesgen/
│   ├── biaee/
│   ├── bicev/
│   ├── paccmann_rl/
│   ├── gxrnn/
│   ├── transgem/
│   ├── triomphe/
│   ├── fame/
│   └── gex2sgen/
│   └── gexmolgen/
├── dataset_processing.py         # Unified download & preprocessing script
├── smilesgen.py                  # SmilesGEN orchestrator runner
├── biaee.py                      # BiAEE orchestrator runner
├── bicev.py                      # BiCEV orchestrator runner
├── paccmann_rl.py                # PaccMann RL orchestrator runner
├── gxrnn.py                      # GxRNN orchestrator runner
├── transgem.py                   # TransGEM orchestrator runner
├── triomphe.py                   # Triomphe orchestrator runner
├── fame.py                       # FAME orchestrator runner (Self-contained)
└── gex2sgen.py                   # GEx2SGen orchestrator runner (Self-contained)
```
