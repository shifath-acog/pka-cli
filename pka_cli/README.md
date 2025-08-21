# pKa Estimation CLI

This command-line interface provides a streamlined way to estimate the pKa of small molecules using various computational chemistry workflows. The tool is designed to be run inside a Docker container to manage its complex dependencies across multiple conda environments.

## Setup

### Step 1: Clone the Repository

```bash
git clone git@github.com:aganitha/pKa-estimation.git
```

### Step 2: Navigate to the Cli Directory

All subsequent commands should be run from the pka-cli project directory.

```bash
cd ~pKa-estimation/pka-cli
```

### Step 3: Run the Docker Container

The container uses the pre-built `pka_pipeline:ovh4` image from the registry. This command will start the container, mount your current pka-cli directory to `/work`, and allocate GPU resources.


```bash
# Make sure you are in the pka-cli directory from the previous step
docker run -d -v "$PWD":/work --gpus all --name pka-cli \
  --label user="$USER" --label description="pka cli" \
  dockerhub.aganitha.ai:4443/chem/pka_pipeline:ovh4 tail -f /dev/null
```

*The `tail -f /dev/null` command keeps the container running in the background.*


## How to Run the CLI

After starting the container, follow these steps to run the pKa estimation workflows.

### Step 1: Enter the Container

Open an interactive shell inside the running `pka-cli` container:

```bash
docker exec -it pka-cli bash
```

### Step 2: Navigate to the Work Directory

All commands should be run from the `/work` directory, which is the root of the mounted project.

```bash
cd /work
```

### Step 3: Execute a Workflow

Run the desired command.

**General Structure:**
```bash
PYTHONPATH=/work python pka_cli/main.py <METHOD> [OPTIONS...]
```

---

## Available Methods & Examples

The CLI supports three different methods for pKa estimation.

### 1. MolConSUL (RDKit/MMFF)

This is the standard workflow using RDKit for conformer generation and the MMFF94 force field for optimization.

**Example 1a: Single SMILES (Default Parameters)**
```bash
PYTHONPATH=/work python pka_cli/main.py molconsul --smiles "c1ccccc1O" --output_dir "cli-test-molconsul-default"
```

**Example 1b: Single SMILES (Custom Parameters)**
```bash
PYTHONPATH=/work python pka_cli/main.py molconsul --smiles "c1ccccc1O" --output_dir "cli-test-molconsul-custom" --num_confs 200 --dielectric 78.3 --pka_exp 9.95
```

**Example 1c: Batch from CSV**
```bash
PYTHONPATH=/work python pka_cli/main.py molconsul --csv "batch_input.csv" --output_dir "cli-test-molconsul-batch"
```

### 2. GNNIS (Graph Neural Network)

This workflow uses the GNNImplicitSolvent model for fast conformer generation. The CLI automatically handles switching to the required `GNNImplicitSolvent` conda environment.

**Example 2a: Single SMILES (Default Parameters)**
```bash
PYTHONPATH=/work python pka_cli/main.py gnnis --smiles "c1ccccc1O" --output_dir "cli-test-gnnis-default"
```

**Example 2b: Single SMILES (Custom Parameters)**
```bash
PYTHONPATH=/work python pka_cli/main.py gnnis --smiles "c1ccccc1O" --output_dir "cli-test-gnnis-custom" --num_confs 25 --solvent_name "DMSO" --pka_exp 9.95
```

**Example 2c: Batch from CSV**
```bash
PYTHONPATH=/work python pka_cli/main.py gnnis --csv "batch_input.csv" --output_dir "cli-test-gnnis-batch"
```

### 3. GRT (Ensemble Method)

This workflow generates a comprehensive set of conformers by combining results from GNNIS, RDKit, and Torsional Diffusion (TD).

**Example 3a: Single SMILES (Default Parameters)**
```bash
PYTHONPATH=/work python pka_cli/main.py grt --smiles "c1ccccc1O" --output_dir "cli-test-grt-default"
```

**Example 3b: Single SMILES (Custom Parameters)**
```bash
PYTHONPATH=/work python pka_cli/main.py grt --smiles "c1ccccc1O" --output_dir "cli-test-grt-custom" --num_confs_td 10 --num_confs_gnnis 10 --num_confs_rdkit 10 --pka_exp 9.95
```

**Example 3c: Batch from CSV**
```bash
PYTHONPATH=/work python pka_cli/main.py grt --csv "batch_input.csv" --output_dir "cli-test-grt-batch"
```
