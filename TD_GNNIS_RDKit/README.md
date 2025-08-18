## 🚀 Quick Start Guide (GRT pipeline)
### GRT = TD + GNNIS + RDKit
### 🔧 Launching the Docker Container

To run the pipeline, pull and start the Docker container using the image below:

```bash
docker pull dockerhub.aganitha.ai:4443/chem/pka_pipeline:ovh5
```

## 📁 Required Files in /work Directory
Ensure the following files and directories are available inside the container under ```/work```:

```bash
/work
│
├── run_all.ipynb                          # Main notebook to be executed (edit parameters here)
├── pka_pipeline_TD.ipynb                  # Conformer generation using Torsional Diffusion
├── full_pipeline_pka_gnnis_1.ipynb        # Conformer generation using GNNIS
├── full_pipeline_pka_rdkit_TD_gnnis_v2.ipynb  # Final step: RDKit-based conformers + pKa estimation
└── csv/                                   #  input CSV with SMILES
```

## ⚙️ Parameter Configuration
Update the following parameters in run_all.ipynb before execution:

```bash
inp_smiles = row['smiles']
num_conf_rdkit: int = 500           # Number of conformers to generate using RDKit
num_conf_TD: int = 500              # Number of conformers to generate using Torsional Diffusion
num_conf_gnnis: int = 500           # Number of conformers to generate using GNNIS
E_avg_proton: float = -275.755
pKa_EXP: float = 24.88
Solvent: str = "DMSO"
dielectric_value = 46.826
```

## ▶️ Running the Pipeline
After updating the parameters, run the run_all.ipynb notebook to begin the full pKa estimation workflow.
