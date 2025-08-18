import os
import subprocess
import pandas as pd
import numpy as np
import shutil
import time
import sys
from typing import List, Tuple, Dict, Optional

from rdkit import Chem
from rdkit.Chem import AllChem

# Assuming GNNImplicitSolvent is installed and available in the environment
from GNNImplicitSolvent import minimize_mol

from atk_conformer_generation_pipeline.utils import *
from atk_conformer_generation_pipeline.variables import *
from pka_cli.utils.common_pka_workflow import run_common_pka_workflow

# =====================================================================
# MODULE 1: GNNIS Conformer Generation
# =====================================================================

def run_gnnis_conformer_generation(smiles: str, num_confs: int, solvent: str, output_dir: str):
    """
    Generates and optimizes conformers using the GNNImplicitSolvent method.
    This function encapsulates the logic from the `full_pipeline_pka_gnnis-1.ipynb` notebook.
    """
    print("Step 1: Generating and Optimizing Conformers with GNNIS...")
    
    # Go into the output directory to run this part
    original_cwd = os.getcwd()
    os.chdir(output_dir)

    try:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, useExpTorsionAnglePrefs=False)

        # The core GNNIS minimization
        minimized_mol, energies = minimize_mol(mol, solvent)

        # Save energies to CSV
        energies_kcal = [e / 4.184 for e in energies]
        energy_df = pd.DataFrame({
            'conformer_id': range(len(energies_kcal)),
            'energy_in_kcalpermol': energies_kcal
        })
        energy_df.to_csv(opt_conf_energy_csv, index=False)
        print(f"Saved GNNIS energies to {opt_conf_energy_csv}")

        # Save conformers to SDF
        save_conformers_to_sdf(minimized_mol, opt_conf_sdf)
        print(f"Saved GNNIS conformers to {opt_conf_sdf}")
        
        os.chdir(original_cwd)
        return minimized_mol, energy_df

    except Exception as e:
        print(f"[ERROR] An error occurred during GNNIS conformer generation: {e}")
        os.chdir(original_cwd)
        return None, None


# =====================================================================
# MAIN WORKFLOW FUNCTION
# =====================================================================

def run_gnnis_workflow(smiles, output_dir, args, base_python_executable):
    """
    Main workflow for the GNNIS pKa estimation pipeline.
    """
    print(f"Running GNNIS workflow for {smiles}")
    print(f"Output will be saved in: {output_dir}")
    
    os.makedirs(output_dir, exist_ok=True)
    
    opt_mol, energy_df = run_gnnis_conformer_generation(
        smiles=smiles,
        num_confs=args.num_confs,
        solvent=args.solvent_name,
        output_dir=output_dir
    )

    if opt_mol is None or energy_df is None:
        print("GNNIS conformer generation failed. Aborting workflow.")
        return

    original_cwd = os.getcwd()
    os.chdir(output_dir)

    # Clean up previous RMSD file to prevent issues
    if os.path.exists(pairwise_RMSDs_dat):
        os.remove(pairwise_RMSDs_dat)

    try:
        print("\nProcessing optimized conformers (Tanimoto, etc.)...")
        convert_conformers_to_smiles(opt_conf_sdf, opt_conf_SMILES_file)
        infeasible_geom_DF, energy_df = process_conformers(
            opt_conf_SMILES_file,
            opt_conf_sdf,
            feasible_geometries_sdf,
            infeasible_geometries_sdf,
            similarity_output_csv,
            infeasible_geometries_csv,
            smiles,
            opt_mol.GetNumConformers(),
            energy_df
        )
        rel_energy_DF = calculate_relative_energies(energy_df, feasible_geometries_csv)

        print("\nStep 2 & 3: Clustering, DFT & pKa Calculation...")
        run_common_pka_workflow(
            opt_mol=opt_mol,
            energy_df=energy_df,
            smiles=smiles,
            output_dir=output_dir,
            dielectric_value=args.dielectric,
            e_avg_proton=args.e_avg_proton,
            pka_exp=args.pka_exp,
            opt_conf_sdf=opt_conf_sdf,
            base_python_executable=base_python_executable
        )

    finally:
        os.chdir(original_cwd)
        print(f"--- GNNIS Workflow for {smiles} complete ---")


