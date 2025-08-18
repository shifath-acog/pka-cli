## 🚀 Getting Started (GRT-ANI-Autograph & GRT-ANI Pipelines)

### 🔧 Launching the Docker Container

To run any of the pipelines, launch a container using the following Docker image:

```bash
docker pull dockerhub.aganitha.ai:4443/chem/pka_pipeline:ovh5
```

---

## 🧠 Pipeline Variants

This repository supports two primary variants of the ANI-based pipeline:

- **GRT-ANI-Autograph (ANI + Autograph)**
- **GRT-ANI (ANI + Hierarchical Clustering)**

Each variant uses its own Jupyter notebook and specific setup.

---

### 📦 Requirements

#### For **GRT-ANI-Autograph**

Before running the pipeline, ensure the following:

- Inside the container, activate the `ani_env` Conda environment.
- Install the required version of `pandas`:

```bash
conda activate ani_env
pip install pandas==1.3.5
```

> ⚠️ This version is required by `Autograph.py`.

---

### 📁 Directory Structure

#### ✅ For **GRT-ANI-Autograph**

Ensure the following files and directories are present in the `/work` directory **inside the container**:

```
/work
│
├── run_all.ipynb                             # Main notebook (edit parameters here)
├── pka_pipeline_TD.ipynb                     # Torsional Diffusion-based conformer generation
├── full_pipeline_pka_gnnis_1.ipynb           # GNNIS-based conformer generation
├── full_pipeline_pka_ani_autograph.ipynb     # Final stage: RDKit-based + ANI optimization + Autograph clustering
├── full_pipeline_pka_DFT.ipynb               # DFT Calculation for pKa estimation
├── Autograph/                                # Clone from: https://github.com/TanemuraKiyoto/AutoGraph
└── csv/                                      # Input directory with CSV file containing SMILES
```

---

#### ✅ For **GRT-ANI**

Ensure the following files and directories are present in the `/work` directory **inside the container**:

```
/work
│
├── run_all.ipynb                             # Main notebook (edit parameters here)
├── pka_pipeline_TD.ipynb                     # Torsional Diffusion-based conformer generation
├── full_pipeline_pka_gnnis_1.ipynb           # GNNIS-based conformer generation
├── full_pipeline_pka_ani.ipynb               # Final stage: RDKit-based + ANI optimization + hierarchical clustering
├── full_pipeline_pka_DFT.ipynb               # DFT Calculation for pKa estimation
└── csv/                                      # Input directory with CSV file containing SMILES
```
