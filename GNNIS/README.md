
## 🚀 Quick Start Guide (GNNIS Pipeline)

### 🔧 Launching the Docker Container

Launch a container using the following Docker image:

```bash
docker pull dockerhub.aganitha.ai:4443/chem/pka_pipeline:ovh5
```

---

### 📁 Inside the Container

Ensure the following files and directories are present in the `/work` directory:

```
/work
│
├── run.ipynb               # Main notebook to execute the pipeline
├── full_pipeline_pka_gnnis-1.ipynb   # Notebook for GNNIS-based conformer generation
├── full_pipeline_pka_gnnis-2_v2.ipynb # Full pipeline for pKa estimation
└── csv                    #  input CSV
```

---

### ▶️ Running the Pipeline

Once inside the container, open and execute the `run.ipynb` notebook to initiate the pipeline.
