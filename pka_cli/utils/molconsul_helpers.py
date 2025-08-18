import os
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from collections import Counter
from scipy.spatial.distance import squareform, is_valid_dm
import glob
import re
import shutil
import time
import sys
from typing import List, Tuple, Dict, Optional


from atk_conformer_generation_pipeline.utils import *
from atk_conformer_generation_pipeline.variables import *

# Import the common workflow for Module 2 & 3
from pka_cli.utils.common_pka_workflow import run_common_pka_workflow

# =====================================================================
# MAIN WORKFLOW FUNCTION
# =====================================================================

def run_molconsul_workflow(smiles, output_dir, args):
    print(f"Running workflow for {smiles}")
    print(f"Output will be saved in: {output_dir}")
    
    os.makedirs(output_dir, exist_ok=True)
    original_cwd = os.getcwd()
    os.chdir(output_dir)

    # Capture the start time
    start_time: float = time.time()

    # Set recursion limit
    sys.setrecursionlimit(10000)

    # Define files and directories to remove (from notebook)
    file_and_dir_to_remove_corrected: List[str]=[
        init_conf_sdf,
        opt_conf_sdf,
        opt_conf_energy_csv,
        opt_conf_SMILES_file,
        similarity_output_csv,
        feasible_geometries_csv,
        infeasible_geometries_csv,
        feasible_geometries_sdf,
        infeasible_geometries_sdf,
        pairwise_RMSDs_dat,
        pairwise_RMSDs_csv,
        cluster_reps_csv,
        cluster_reps_sdf,
        cluster_reps_dir,
        clusters_RMSD_stats_csv,
        clusters_energy_stats_csv,
        opt_cluster_reps_csv
    ]
    remove_paths(file_and_dir_to_remove_corrected)

    # --- MODULE 1 ---
    print("Step 1: Generating and Optimizing Conformers...")
    initial_mol = generate_conformers(smiles, args.num_confs)
    save_conformers_to_sdf(initial_mol, init_conf_sdf) # Added this line as per notebook

    opt_mol, conformer_energies = mmff_optimize_conformers(initial_mol)
    save_conformers_to_sdf(opt_mol, opt_conf_sdf)
    energy_df = pd.DataFrame(list(conformer_energies.items()), columns=['conformer_id', 'energy_in_kcalpermol'])
    energy_df.to_csv(opt_conf_energy_csv, index=False)
    print(f"Successfully generated and optimized {opt_mol.GetNumConformers()} conformers.")

    # Added missing Module 1 steps from notebook
    print("\nProcessing optimized conformers...")
    convert_conformers_to_smiles(opt_conf_sdf, opt_conf_SMILES_file)
    infeasible_geom_DF, energy_df = process_conformers(
        opt_conf_SMILES_file,
        opt_conf_sdf,
        feasible_geometries_sdf,
        infeasible_geometries_sdf,
        similarity_output_csv,
        infeasible_geometries_csv,
        smiles, # inp_smiles from args
        opt_mol.GetNumConformers(), # num_opt_conf
        energy_df
    )
    rel_energy_DF = calculate_relative_energies(energy_df, feasible_geometries_csv) # This output CSV name is confusing, notebook uses it for feasible_geom_energies.csv


    # --- MODULE 2 & 3: Use the common workflow ---
    print("\nStep 2 & 3: Clustering, DFT & pKa Calculation...")
    run_common_pka_workflow(
        opt_mol,
        energy_df,
        smiles,
        output_dir,
        args.dielectric,
        args.e_avg_proton,
        args.pka_exp,
        opt_conf_sdf
    )

    os.chdir(original_cwd)
    
    # Capture the end time
    end_time: float = time.time()
    # Calculate the execution time in seconds
    execution_time_seconds: float = end_time - start_time
    # Convert the execution time to minutes
    execution_time_minutes: int = int(execution_time_seconds // 60)
    print(f"Total execution time: {execution_time_minutes} minutes.")
    print(f"--- Workflow for {smiles} complete ---")

