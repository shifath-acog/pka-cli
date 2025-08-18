## 🚀 Quick Start Guide (RDKit Pipeline)

### 🔧 Launching the Docker Container

Launch a container using the following Docker image:

```bash
docker pull dockerhub.aganitha.ai:4443/chem/pka_pipeline:ovh5
```

---

### 📁 Inside the Container
Ensure the following files and directories are present in the /work directory:
```bash
/work
│
├── run_scr.ipynb               # The main notebook to be executed
├── full_pipeline_pka_v0.3.ipynb # The second notebook (modify name as needed)
└── csv                       # input CSV file with SMILES
```
---

### ▶️ Running the Pipeline

Once inside the container, open and execute the `run_scr.ipynb` notebook to initiate the pipeline.

