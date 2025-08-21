### 1. What steps is the notebook running in order and are we following the same for our CLI?

Yes, the CLI is meticulously designed to replicate the exact workflow of the `MolConSUL/full_pipeline_pka_v0.3.ipynb` notebook. The process is divided into three main modules:

**Module 1: Conformer Generation and Optimization**
*   **Input & Setup:** Define molecular parameters (SMILES, dielectric, conformer count) and perform initial cleanup of previous run data.
*   **Conformer Generation (RDKit-ETKDG):** Generate initial 3D conformers.
*   **MMFF Optimization:** Optimize conformer geometries using the MMFF94 force field.
*   **Data Saving:** Save optimized conformers to SDF and their energies to CSV.
*   **SMILES Conversion:** Convert 3D conformers back to SMILES strings.
*   **Feasibility Check (Tanimoto Similarity):** Filter conformers based on Tanimoto similarity to the input SMILES, separating feasible from infeasible geometries and updating energy data.
*   **Relative Energy Calculation:** Compute relative MMFF energies for feasible conformers.

**Module 2: Clustering**
*   **RMSD Matrix Calculation:** Use Open Babel's `obrms` to compute pairwise RMSD for feasible conformers.
*   **RMSD Matrix Processing:** Prepare the RMSD data for clustering.
*   **Hierarchical Clustering:** Group conformers using hierarchical clustering.
*   **Optimal Cluster Determination:** Identify the optimal number of clusters using silhouette scores.
*   **Cluster Representative Identification:** Select the minimum energy conformer from each cluster as its representative.
*   **Representative Saving:** Save cluster representatives to individual and concatenated SDF files.

**Module 3: DFT & pKa Calculation**
*   **Neutral DFT Calculation:** Run `dft_main.py` to perform DFT optimization and energy calculation on neutral cluster representatives, generating `.xyz` and `.charge` files.
*   **Crucial SDF Renaming:** **(This was the missing piece in the CLI initially!)** The notebook renames the cluster representative `.sdf` files to include the basis set and XC functional (e.g., `rep_of_cluster_X_aug-cc-pVDZ_M06-2X.sdf`). This is done by copying the original `.sdf` files to new names.
*   **Ionized Species Generation:** Run `ionization_tool_with_ranking/main.py` to generate deprotonated/protonated forms from the *renamed* `.sdf` files and their corresponding `.charge` files. This outputs new ionized `.sdf` files into a `treated` subdirectory.
*   **Move Ionized Species:** Transfer the newly generated ionized `.sdf` files from the `treated` subdirectory to the main `cluster_rep_conformers` directory.
*   **Charged DFT Calculation:** Run `dft_main.py` again on the ionized `.sdf` files to perform DFT optimization and energy calculation, generating ionized `.xyz` files.
*   **DFT Data Accumulation:** Collect all `.xyz` files (neutral, deprotonated, protonated) and extract their energies, basis sets, XC functionals, and ionization states into a pandas DataFrame.
*   **pKa Calculation:** Group data by conformer, apply Boltzmann weighting, and calculate acidic/basic pKa values.
*   **Result Saving:** Save detailed and summary pKa results to CSV files.
*   **Execution Time:** Report the total pipeline execution time.

The CLI now successfully executes all these steps, ensuring full parity with the notebook's scientific workflow.

### 2. In the CLI run command, why did you use the ionization path (`/work/ionization_tool_with_ranking`) in `PYTHONPATH`? Can't we do it like the notebook does?

You're asking about the `-e PYTHONPATH=/work:/work/ionization_tool_with_ranking` part of the `docker exec` command.

The notebook's initial cells often include `sys.path.append('/work')`. When a Python script (like `ionization_tool_with_ranking/main.py`) is executed directly, its own directory is automatically added to `sys.path` as the first entry. So, when `main.py` runs, `/work/ionization_tool_with_ranking` is implicitly in `sys.path`.

The issue arises because `ionization_tool_with_ranking/main.py` imports modules like `protodeproto.acidic`. `protodeproto` is a *sub-package* located within the `ionization_tool_with_ranking` directory (i.e., `/work/ionization_tool_with_ranking/protodeproto`).

If `PYTHONPATH` only contained `/work`, Python would search for `protodeproto` directly under `/work` (i.e., `/work/protodeproto`), which doesn't exist. It needs to find `protodeproto` as a sub-package of `ionization_tool_with_ranking`.

By explicitly setting `PYTHONPATH=/work:/work/ionization_tool_with_ranking`, we ensure that:
1.  `/work` is in the path, allowing imports like `atk_conformer_generation_pipeline.utils`.
2.  `/work/ionization_tool_with_ranking` is also in the path, allowing `main.py` to correctly locate and import `protodeproto` as a top-level package.

While Jupyter notebooks can sometimes handle these paths more flexibly due to their execution environment, for a direct `subprocess.run` call, being explicit with `PYTHONPATH` guarantees that all necessary modules are found.

### 3. For the SDF files, what was the notebook doing? Was it after generating those SDF files, was it resaving them with added suffix? And which code was generating these SDF, XYZ, and charge files? Help me understand all this questions.

This was indeed a crucial point of confusion and the source of many errors! You've hit on a critical detail of the workflow.

Here's the detailed breakdown of how SDF, XYZ, and charge files are generated and handled:

1.  **Initial Conformer SDF (`initial_generated_conformers.sdf`):
    *   **Generated by:** `generate_conformers()` (from `atk_conformer_generation_pipeline/utils.py`) and saved by `save_conformers_to_sdf()`.
    *   **Purpose:** Stores the raw 3D conformers from RDKit.

2.  **Optimized Conformer SDF (`optimized_generated_conformers.sdf`):
    *   **Generated by:** `mmff_optimize_conformers()` (from `atk_conformer_generation_pipeline/utils.py`) and saved by `save_conformers_to_sdf()`.
    *   **Purpose:** Stores the MMFF-optimized 3D conformers.

3.  **Feasible/Infeasible Geometries SDFs (`feasible_geometries.sdf`, `infeasible_geometries.sdf`):
    *   **Generated by:** `process_conformers()` (from `atk_conformer_generation_pipeline/utils.py`).
    *   **Purpose:** These are subsets of the optimized conformers, filtered based on Tanimoto similarity. Only feasible geometries are passed to the next steps.

4.  **Cluster Representative SDFs (e.g., `rep_of_cluster_1.sdf`):
    *   **Generated by:** `write_cluster_representatives()` (from `atk_conformer_generation_pipeline/utils.py`).
    *   **Purpose:** For each cluster, the minimum energy conformer is selected and saved as an individual `.sdf` file (e.g., `rep_of_cluster_1.sdf`, `rep_of_cluster_2.sdf`, etc.) within the `cluster_rep_conformers` directory. A concatenated `cluster_rep_conformers.sdf` is also created.
    *   **Crucial Detail:** At this stage, these individual `.sdf` files **do NOT** have the basis set or XC functional in their names. They are simply `rep_of_cluster_X.sdf`.

5.  **Neutral XYZ and Charge Files (e.g., `rep_of_cluster_1_aug-cc-pVDZ_M06-2X.xyz`, `rep_of_cluster_1_aug-cc-pVDZ_M06-2X.charge`):
    *   **Generated by:** The **first call** to `dft_main.py` (`!python /work/dft_main.py DFT ...`).
    *   **Purpose:** `dft_main.py` reads the `rep_of_cluster_X.sdf` files, performs DFT geometry optimization and energy calculation, and outputs the results as `.xyz` files (containing optimized geometry and energy) and `.charge` files (containing Mulliken charges). These files **DO** include the basis set and XC functional in their names.

6.  **SDF Renaming/Copying (The "Resaving with Suffix" Part!):
    *   **What the notebook was doing:** This was the **key missing step** in our initial CLI implementation that caused many errors. Immediately after the first `dft_main.py` call, the notebook has a block of Python code (which we added to `run_dft_calculation` in the CLI) that iterates through the `rep_of_cluster_X.sdf` files. For each, it **copies** (not resaves, but creates a new file) it to a new filename that **includes the basis set and XC functional** (e.g., `rep_of_cluster_1.sdf` is copied to `rep_of_cluster_1_aug-cc-pVDZ_M06-2X.sdf`).
    *   **Purpose:** This step is critical because the `ionization_tool_with_ranking/main.py` script expects its input `.sdf` files to have the basis set and XC functional in their names so it can correctly find their corresponding `.charge` files (which were generated by `dft_main.py` with those suffixes).

7.  **Ionized SDF Generation (e.g., `rep_of_cluster_1_deprotonated.sdf`):
    *   **Generated by:** The call to `ionization_tool_with_ranking/main.py`.
    *   **Purpose:** This script reads the **renamed** `.sdf` files (e.g., `rep_of_cluster_1_aug-cc-pVDZ_M06-2X.sdf`) and their corresponding `.charge` files. It then uses the `protodeproto` library to generate deprotonated and/or protonated forms, saving them as new `.sdf` files (e.g., `rep_of_cluster_1_deprotonated.sdf`) in a `treated` subdirectory.

8.  **Ionized XYZ File Generation (e.g., `rep_of_cluster_1_deprotonated.xyz`):
    *   **Generated by:** The **second call** to `dft_main.py` (`!python /work/dft_main.py DFT_DEPROTO ...`).
    *   **Purpose:** `dft_main.py` reads the newly generated ionized `.sdf` files (e.g., `rep_of_cluster_1_deprotonated.sdf`) and performs DFT optimization and energy calculation. It outputs `.xyz` files for these ionized species.
    *   **Important Note:** `dft_main.py` does **NOT** add basis/XC information to these ionized `.xyz` filenames. They remain simple (e.g., `rep_of_cluster_1_deprotonated.xyz`). This was another point of mismatch we resolved by adjusting the regex in `create_dataframe_from_xyz_files`.

This intricate dance of file naming and processing was indeed the root cause of many of the `FileNotFoundError` and `Filename pattern did not match` issues we faced.

### How to Run the CLI

The pKa estimation CLI is designed to be run within a Docker container named `updated-pka`. Ensure your current working directory (where this `info.md` file is located) is mounted to `/work` inside the container.

**If you are already `exec`'d into the Docker container:**

If you have already entered the Docker container (e.g., using `docker exec -it updated-pka bash`), you can run the CLI commands directly without the `docker exec` prefix and without specifying `/work` for paths, as your current working directory inside the container is `/work`.

**General Command Structure (inside container):**

```bash
PYTHONPATH=/work:/work/ionization_tool_with_ranking python pka_cli/main.py molconsul \
  --smiles "<SMILES_STRING>" \
  --output_dir "<OUTPUT_FOLDER_NAME>" \
  [--dielectric <VALUE>] \
  [--num_confs <VALUE>] \
  [--e_avg_proton <VALUE>] \
  [--pka_exp <VALUE>]
```

Or for batch processing from a CSV file (inside container):

```bash
PYTHONPATH=/work:/work/ionization_tool_with_ranking python pka_cli/main.py molconsul \
  --csv "<PATH_TO_CSV_FILE>" \
  --output_dir "<OUTPUT_FOLDER_NAME>" \
  [--dielectric <VALUE>] \
  [--num_confs <VALUE>] \
  [--e_avg_proton <VALUE>] \
  [--pka_exp <VALUE>]
```

**Example Commands (inside container):**

1.  **Run for a single SMILES string:**
    ```bash
    PYTHONPATH=/work:/work/ionization_tool_with_ranking python pka_cli/main.py molconsul \
      --smiles "NC(=O)CO" \
      --output_dir "my_molecule_output"
    ```

2.  **Run for a single SMILES string with custom parameters:**
    ```bash
    PYTHONPATH=/work:/work/ionization_tool_with_ranking python pka_cli/main.py molconsul \
      --smiles "CCO" \
      --output_dir "ethanol_pka" \
      --dielectric 78.355 \
      --num_confs 100 \
      --e_avg_proton -270.28 \
      --pka_exp 15.9
    ```

3.  **Run for a batch of SMILES strings from a CSV file:**
    (Assuming you have a `molecules.csv` file in your current working directory inside the container)
    ```bash
    PYTHONPATH=/work:/work/ionization_tool_with_ranking python pka_cli/main.py molconsul \
      --csv "molecules.csv" \
      --output_dir "batch_pka_results"
    ```

**Important Notes:**
*   The `--output_dir` path will be created relative to your current working directory inside the container (which is `/work`).
*   The `--csv` path should also be relative to your current working directory inside the container.
*   The `updated-pka` Docker container must be running and your host machine's current working directory must be mounted to `/work` within the container.
*   Optional parameters (`--dielectric`, `--num_confs`, `--e_avg_proton`, `--pka_exp`) can be adjusted as needed. If not provided, default values will be used.

**From outside the container:**

```bash
docker exec updated-pka bash -c 'PYTHONPATH=/work python3 /work/pka_cli/main.py molconsul --smiles "NC(=O)CO" --output_dir "cli_test_gemini_2" --dielectric 46.826 --num_confs 500 --e_avg_proton -277.60 --pKa_exp 23.00'
```
