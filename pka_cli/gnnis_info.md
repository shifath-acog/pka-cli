# GNNIS pKa Estimation Workflow

This document provides an overview of the GNNIS (Graph Neural Network Implicit Solvent) pKa estimation workflow as implemented in the `pka_cli` tool.

## Implementation Overview

The GNNIS workflow is a two-part process that leverages both a specialized GNN model for conformer generation and a common workflow for subsequent calculations. The implementation is designed to handle the specific conda environment requirements of the GNNIS library.

### Environment Handling

The GNNIS workflow requires the `GNNImplicitSolvent` conda environment. The `pka_cli` tool handles this automatically:

1.  When you execute the `gnnis` command, the script first checks if it is running within the `GNNImplicitSolvent` environment.
2.  If it is not, the script will automatically re-launch itself using `conda run -n GNNImplicitSolvent ...`.
3.  This ensures that the GNNIS-specific code is executed with the correct dependencies, while the rest of the CLI remains in the base environment.

### Workflow Steps

1.  **Module 1: Conformer Generation (GNNIS)**
    *   This step is handled by the `run_gnnis_conformer_generation` function in `pka_cli/utils/gnnis_helpers.py`.
    *   It takes a SMILES string and the number of conformers to generate as input.
    *   It uses the `minimize_mol` function from the `GNNImplicitSolvent` library to generate and optimize the conformers.
    *   The output of this step is an SDF file with the optimized conformers and a CSV file with their energies.

2.  **Module 2 & 3: Clustering, DFT, and pKa Calculation**
    *   Once the conformers are generated, the workflow proceeds to the common clustering, DFT, and pKa calculation steps.
    *   These steps are handled by the `run_common_pka_workflow` function in `pka_cli/utils/common_pka_workflow.py`.
    *   This part of the workflow is shared with the `molconsul` command.

## Example Commands

To run the GNNIS workflow, use the `gnnis` subcommand. Here are some examples:

**Run for a single SMILES string:**

```bash
docker exec updated-pka bash -c 'PYTHONPATH=/work python3 /work/pka_cli/main.py gnnis --smiles "NC(=O)CO" --output_dir "cli-gnnis-test"'
```

**Run with a different number of conformers:**

*Note: The DFT calculation in the GNNIS workflow is memory-intensive. If you encounter out-of-memory errors, try reducing the number of conformers.* 

```bash
docker exec updated-pka bash -c 'PYTHONPATH=/work python3 /work/pka_cli/main.py gnnis --smiles "NC(=O)CO" --output_dir "cli-gnnis-test" --num_confs 10'
```

**Run for a batch of SMILES from a CSV file:**

*The CSV file must contain 'smiles' and 'id' columns.*

```bash
docker exec updated-pka bash -c 'PYTHONPATH=/work python3 /work/pka_cli/main.py gnnis --csv "path/to/your/file.csv" --output_dir "cli-gnnis-batch-test"'
```
