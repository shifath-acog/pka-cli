import os
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from sklearn.metrics import silhouette_score
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

# =====================================================================
# MODULE 2: Clustering
# =====================================================================

def perform_clustering_and_select_reps(rmsd_matrix_file: str, energy_df: pd.DataFrame, opt_mol: Chem.Mol, num_clusters_to_pick: int = 5):
    rmsd_df = pd.read_csv(rmsd_matrix_file, header=None, index_col=0)
    rmsd_matrix = np.round(rmsd_df.to_numpy(dtype=float), 2)
    rmsd_matrix = (rmsd_matrix + rmsd_matrix.T) / 2
    condensed_matrix = squareform(rmsd_matrix)
    linkage_matrix = sch.linkage(condensed_matrix, method='ward')
    
    sil_scores = []
    upper_range = min(101, len(energy_df))
    if upper_range <= 2:
        optimal_clusters = 2
    else:
        range_n_clusters = list(range(2, upper_range))
        for n_clusters in range_n_clusters:
            labels = sch.fcluster(linkage_matrix, n_clusters, criterion='maxclust')
            if len(set(labels)) > 1:
                sil_scores.append(silhouette_score(rmsd_matrix, labels, metric='precomputed'))
            else:
                sil_scores.append(-1)
        if not sil_scores:
             optimal_clusters = 2
        else:
             optimal_clusters = range_n_clusters[np.argmax(sil_scores)]

    num_clusters = 5 if optimal_clusters > 5 else optimal_clusters
    print(f"Optimal number of clusters found: {optimal_clusters}. Selecting {num_clusters}.")

    cluster_labels = sch.fcluster(linkage_matrix, num_clusters, criterion='maxclust')
    clusters = {i: [] for i in range(1, num_clusters + 1)}
    for index, label in enumerate(cluster_labels):
        clusters[label].append(index)

    cluster_reps_list = []
    for clust_label, elements in clusters.items():
        if elements:
            clust_df = energy_df.loc[elements]
            min_energy_index = clust_df['energy_in_kcalpermol'].idxmin()
            min_energy_df = clust_df.loc[[min_energy_index]].copy()
            min_energy_df['cluster_id'] = clust_label
            cluster_reps_list.append(min_energy_df)
    
    cluster_reps_df = pd.concat(cluster_reps_list, ignore_index=True)
    return cluster_reps_df.sort_values(by='conformer_id', ascending=True)

# =====================================================================
# MODULE 3: DFT & pKa Calculation
# =====================================================================

def run_dft_calculation(cluster_reps_dir, dielectric_value, python_executable):
    print("Running DFT on neutral conformers...")
    dft_command_neutral = [python_executable, "/work/dft_main.py", "DFT", cluster_reps_dir, str(dielectric_value)]
    try:
        subprocess.run(dft_command_neutral, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] DFT calculation failed for neutral species.")
        print(f"COMMAND: {' '.join(e.cmd)}")
        print(f"STDOUT:\n{e.stdout}")
        print(f"STDERR:\n{e.stderr}")
        raise

    # Renaming SDF files as per notebook workflow
    basis_sets = ['aug-cc-pVDZ'] #aug-cc-pVDZ
    xc_functionals = ['M06-2X']

    sdf_files = sorted(glob.glob(os.path.join(cluster_reps_dir, "*.sdf")))

    for sdf_path in sdf_files:
        sdf_filename = os.path.basename(sdf_path)
        sdf_stem = os.path.splitext(sdf_filename)[0]

        for basis in basis_sets:
            for xc in xc_functionals:
                basis_tag = basis.replace("(", "").replace(")", "").replace("*", "").replace("+", "").replace("/", "-").replace(" ", "")
                xc_tag = xc.replace("(", "").replace(")", "").replace("*", "").replace("+", "").replace("/", "-").replace(" ", "")
                
                new_filename = f"{sdf_stem}_{basis_tag}_{xc_tag}.sdf"
                new_filepath = os.path.join(cluster_reps_dir, new_filename)
                shutil.copy(sdf_path, new_filepath)
                print(f"Renamed {sdf_filename} to {new_filename}")

    print("Generating ionized species...")
    ionization_command = [python_executable, os.path.join("/work", "ionization_tool_with_ranking", "main.py"), os.path.abspath(cluster_reps_dir)]
    try:
        subprocess.run(ionization_command, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Ionization failed.")
        print(f"COMMAND: {' '.join(e.cmd)}")
        print(f"STDOUT:\n{e.stdout}")
        print(f"STDERR:\n{e.stderr}")
        raise
    
    treated_dir = os.path.join(cluster_reps_dir, "treated")
    if os.path.exists(treated_dir):
        for file_name in os.listdir(treated_dir):
            shutil.move(os.path.join(treated_dir, file_name), cluster_reps_dir)
        shutil.rmtree(treated_dir)

    print("Running DFT on charged conformers...")
    dft_command_charged = [python_executable, "/work/dft_main.py", "DFT_DEPROTO", cluster_reps_dir, str(dielectric_value)]
    try:
        subprocess.run(dft_command_charged, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] DFT calculation failed for charged species.")
        print(f"COMMAND: {' '.join(e.cmd)}")
        print(f"STDOUT:\n{e.stdout}")
        print(f"STDERR:\n{e.stderr}")
        raise

def create_dataframe_from_xyz_files(folder_path):
    """
    Creates a pandas DataFrame from .xyz files in a given folder,
    extracting Basis, XC, and energy.

    Args:
        folder_path (str): The path to the folder containing the .xyz files.

    Returns:
        pandas.DataFrame: A DataFrame with 'file name', 'Basis', 'XC', and 'energy' columns.
    """
    data = []

    # Updated Regex to match the filename pattern and capture Basis and XC
    # This regex is more robust to handle basis sets and XC functionals
    # that may contain hyphens, numbers, or specific characters.
    # It assumes the structure is:
    # rep_of_cluster_N_BASIS_XC_optional_treatment.xyz
    filename_pattern = re.compile(
        r"rep_of_cluster_(\d+)"  # Cluster ID (Group 1)
        r"(?:_((?!deprotonated|protonated)[a-zA-Z0-9\\-]+)_([a-zA-Z0-9\\-]+))?"  # Optional Basis (Group 2) and XC (Group 3) - negative lookahead
        r"(?:_(deprotonated|protonated))?"  # Optional Ionization State (Group 4)
        r"\.xyz"
    )
    for filename in os.listdir(folder_path):
        if filename.endswith(".xyz"):
            match = filename_pattern.match(filename)
            if match:
                cluster_id = match.group(1)
                basis = match.group(2)
                xc = match.group(3)
                ionization_state = match.group(4) # This will be 'deprotonated', 'protonated', or None

                file_path = os.path.join(folder_path, filename)
                try:
                    with open(file_path, 'r') as f:
                        f.readline()  # Skip the first line (number of atoms)
                        energy_line = f.readline().strip()

                        energy = None
                        try:
                            energy = float(energy_line)
                        except ValueError:
                            energy_match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", energy_line)
                            if energy_match:
                                energy = float(energy_match.group(0))
                            else:
                                print(f"Could not extract energy from line: '{energy_line}' in file: {filename}")

                        data.append({
                            'file_name': filename,
                            'cluster_id': cluster_id,
                            'Basis': basis,
                            'XC': xc,
                            'ionization_state': ionization_state,
                            'energy_kcal_mol': energy
                        })
                except Exception as e:
                    print(f"Error processing file {filename}: {e}")
            else:
                print(f"Filename pattern did not match for: {filename}")

    df = pd.DataFrame(data)
    return df

def get_prefix(filename):
    # This regex extracts the base part of the filename for grouping
    # It should match 'rep_of_cluster_X_aug-cc-pVDZ_M06-2X' for neutral
    # and 'rep_of_cluster_X' for ionized
    match_neutral = re.match(r'(rep_of_cluster_\d+_aug-cc-pVDZ_M06-2X)', filename)
    if match_neutral:
        return match_neutral.group(1)
    
    match_ionized = re.match(r'(rep_of_cluster_\d+)', filename)
    if match_ionized:
        return match_ionized.group(1) + "_aug-cc-pVDZ_M06-2X" # Reconstruct the full prefix for grouping
    
    return filename # Fallback, though ideally all files should match

def calculate_pka(energy_df, e_avg_proton, pka_exp):
    if energy_df.empty:
        print("Warning: DFT energy DataFrame is empty. Cannot calculate pKa.")
        return np.nan, np.nan, pd.DataFrame()

    # Group by cluster_id to process each conformer's neutral and ionized forms
    output_rows = []
    for cluster_id, group in energy_df.groupby('cluster_id'):
        row_dict = {'Cluster_ID': cluster_id}
        
        # Initialize placeholders
        row_dict['Filename_original'] = None
        row_dict['Energy_original'] = None
        row_dict['Filename_deprotonated'] = None
        row_dict['Energy_deprotonated'] = None
        row_dict['Filename_protonated'] = None
        row_dict['Energy_protonated'] = None
        row_dict['Basis'] = None
        row_dict['XC'] = None

        for _, entry in group.iterrows():
            if entry['ionization_state'] is None: # This is the neutral species
                row_dict['Filename_original'] = entry['file_name']
                row_dict['Energy_original'] = entry['energy_kcal_mol']
                row_dict['Basis'] = entry['Basis']
                row_dict['XC'] = entry['XC']
            elif entry['ionization_state'] == 'deprotonated':
                row_dict['Filename_deprotonated'] = entry['file_name']
                row_dict['Energy_deprotonated'] = entry['energy_kcal_mol']
            elif entry['ionization_state'] == 'protonated':
                row_dict['Filename_protonated'] = entry['file_name']
                row_dict['Energy_protonated'] = entry['energy_kcal_mol']
        output_rows.append(row_dict)

    result_df = pd.DataFrame(output_rows)

    # Filter for only the columns that actually exist in result_df to avoid errors
    # This part remains the same as before, but now uses 'Cluster_ID' instead of 'Prefix'
    desired_columns = ['Cluster_ID', 'Basis', 'XC',
                       'Filename_original', 'Energy_original',
                       'Filename_deprotonated', 'Energy_deprotonated',
                       'Filename_protonated', 'Energy_protonated']

    existing_columns = [col for col in desired_columns if col in result_df.columns]
    result_df = result_df[existing_columns]

    # Boltzmann weighting
    T = 298  # Temperature in Kelvin
    k_B = 0.0019872041  # Boltzmann constant in kcal/mol·K

    if 'Energy_original' not in result_df.columns:
        print("Error: Original (neutral) conformer energies not found for pKa calculation.")
        return np.nan, np.nan, pd.DataFrame()

    E_shifted = result_df['Energy_original'] - result_df['Energy_original'].min()
    weights = np.exp(-E_shifted / (k_B * T))
    result_df['Weights'] = weights / np.sum(weights)

    pka_cal_acidic, pka_cal_basic = np.nan, np.nan

    if 'Energy_deprotonated' in result_df.columns:
        result_df.loc[:, 'pka (acidic)'] = 0.7322 * (
            e_avg_proton - (result_df['Energy_original'] - result_df['Energy_deprotonated'])
        )
        pka_cal_acidic = np.sum(result_df['pka (acidic)'] * result_df['Weights'])

    if 'Energy_protonated' in result_df.columns:
        result_df.loc[:, 'pka (basic)'] = 0.7322 * (
            e_avg_proton + (result_df['Energy_original'] - result_df['Energy_protonated'])
        )
        pka_cal_basic = np.sum(result_df['pka (basic)'] * result_df['Weights'])

    return pka_cal_acidic, pka_cal_basic, result_df

# =====================================================================
# MAIN WORKFLOW FUNCTION (Common Part)
# =====================================================================

def run_common_pka_workflow(opt_mol, energy_df, smiles, output_dir, dielectric_value, e_avg_proton, pka_exp, opt_conf_sdf, base_python_executable=None):
    # This function encapsulates Module 2 and Module 3
    # It expects opt_mol and energy_df from Module 1

    python_executable = base_python_executable or sys.executable

    # --- MODULE 2 ---
    print("\nStep 2: Clustering Conformers...")
    os.makedirs(cluster_reps_dir, exist_ok=True)

    calculate_rmsd(feasible_geometries_sdf, pairwise_RMSDs_dat)
    cluster_reps_df = perform_clustering_and_select_reps(pairwise_RMSDs_dat, energy_df, opt_mol)
    
    if cluster_reps_df.empty or len(cluster_reps_df) == 0:
        print("Clustering resulted in no representative conformers. Aborting pKa calculation.")
        return

    cluster_reps_df.to_csv(cluster_reps_csv, index=False)
    print(f"Found {len(cluster_reps_df)} representative conformers.")

    # This call now uses the passed-in argument, which is much safer
    write_cluster_representatives(opt_conf_sdf, cluster_reps_dir, cluster_reps_sdf, cluster_reps_df, cluster_reps_df, "rep_of_cluster_", ".sdf")

    # --- MODULE 3 ---
    print("\nStep 3: DFT & pKa Calculation...")
    run_dft_calculation(cluster_reps_dir, dielectric_value, python_executable)
    
    # Use the corrected function to create the dataframe
    dft_energy_df = create_dataframe_from_xyz_files(cluster_reps_dir)
    
    # Check if the dataframe is empty after the fix
    if dft_energy_df.empty:
        print("\nERROR: DFT energy dataframe is still empty after regex fix. Check DFT output files in the cluster_rep_conformers directory.")
        return

    pka_acidic, pka_basic, final_df = calculate_pka(dft_energy_df, e_avg_proton, pka_exp)
    
    final_df.to_csv("conformer_pka_results.csv", index=False)
    
    summary = {
        "Molecule": [smiles],
        "Solvent (dielectric)": [dielectric_value],
        "Num_Cluster_Reps": [len(cluster_reps_df)],
        "pKa_Exp": [pka_exp],
        "pKa_Calc_Acidic": [round(pka_acidic, 2) if not np.isnan(pka_acidic) else 'N/A'],
        "pKa_Calc_Basic": [round(pka_basic, 2) if not np.isnan(pka_basic) else 'N/A']
    }
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv("final_pka_summary.csv", index=False)
    
    print("\n--- Final Results ---")
    print(summary_df.to_string())


#
