import os
import shutil
import time
from typing import *
import subprocess
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from IPython.display import Image, display
from rdkit.Chem import SDWriter
from rdkit import DataStructs
import pandas as pd

def time_rand_seed() -> int:
    """Generate a random seed based on the current time.

    Returns:
        int: A random seed derived from the decimal part of the current time.
    """
    # Get the current time
    t: float = time.time()
    
    # Isolate the decimal part of current time
    t -= int(t)
    
    # Multiply the decimal part of current time with 1.0e6 and take only the integer part of the product
    seed: int = int(t * 1.0e6)
    
    # Return seed as output
    return seed

def display_2d_structure(smiles: str, filename: str = 'molecule.png'):
    """
    Generate and display the 2D structure of a molecule given its SMILES string.

    Args:
        smiles (str): The SMILES string of the molecule.
        filename (str): The filename to save the image. Default is 'molecule.png'.
    """
    # Convert SMILES string to RDKit molecule object
    mol: Chem.Mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")

    # Generate 2D coordinates
    AllChem.Compute2DCoords(mol)

    # Draw the molecule and save to file
    img = Draw.MolToImage(mol, size=(300, 300))
    img.save(filename)
    display(Image(filename))

def remove_path(path: str) -> None:
    """Remove a file or directory if it exists.

    Args:
        path (str): The path to the file or directory to remove.
    """
    # Check if the path exists
    if os.path.exists(path):
        # Check if the path is a file
        if os.path.isfile(path):
            # Remove the file
            os.remove(path)
            print(f"File {path} has been removed.")
        # Check if the path is a directory
        elif os.path.isdir(path):
            # Remove the directory and its contents
            shutil.rmtree(path)
            print(f"Directory {path} has been removed.")
    else:
        # Path does not exist
        print(f"The path {path} does not exist.")

def remove_paths(paths: List[str]) -> None:
    """Remove a files or directories if it exists.

    Args:
        path List[str]: The paths to the files or directories to remove.
    """
    for path in paths:
        # Check if the path exists
        if os.path.exists(path):
            # Check if the path is a file
            if os.path.isfile(path):
                # Remove the file
                os.remove(path)
                print(f"File {path} has been removed.")
            # Check if the path is a directory
            elif os.path.isdir(path):
                # Remove the directory and its contents
                shutil.rmtree(path)
                print(f"Directory {path} has been removed.")
        else:
            # Path does not exist
            print(f"The path {path} does not exist.")

def generate_conformers(smiles: str, num_conf: int,num_threads: int=10) -> Chem.Mol:
    """Generate conformers using RDKit-ETKDG.

    Args:
        smiles (str): The SMILES string representing the molecule.
        num_conf (int): The number of conformers to generate.

    Returns:
        Chem.Mol: The molecule object with generated conformers.
    """
    # Convert the SMILES string to a molecule
    mol: Chem.Mol = Chem.MolFromSmiles(smiles)
    
    # Add hydrogens to the molecule, which is essential to get good structures
    mol = Chem.AddHs(mol)
    
    # Embed multiple conformers using ETKDG
    params: AllChem.EmbedParameters = AllChem.ETKDGv3()  # EmbedParameters object for the ETKDG method - version 3 (macrocycles)
    params.randomSeed = time_rand_seed()  # Provide random seed for the generator based on the current time
    params.numThreads = num_threads  # Use ten threads
    
    # params.pruneRmsThresh = 0.1  # Prune conformers that are too similar
    
    # Generate conformers
    AllChem.EmbedMultipleConfs(mol, numConfs=num_conf, params=params)

    # Return the molecule object with generated conformers
    return mol

def mmff_optimize_conformers(opt_mol: Chem.Mol) -> Tuple[Chem.Mol, Dict[int, float]]:
    """Optimize the conformers of a molecule using the MMFF94 force field.

    Args:
        opt_mol (Chem.Mol): The molecule with conformers to optimize.

    Returns:
        Tuple[Chem.Mol, Dict[int, float]]: The optimized molecule and a dictionary
                                           containing MMFF energies of the conformers.
    """
    # An empty dictionary to store MMFF energies of optimized conformers
    mmff_energies: Dict[int, float] = {}

    # Loop over each conformer in the molecule
    for conf_id in range(opt_mol.GetNumConformers()):
        # Calculate MMFF properties
        mmff_props = AllChem.MMFFGetMoleculeProperties(opt_mol, mmffVariant='MMFF94')
        
        if mmff_props is None:
            # Raise an error if MMFF properties cannot be calculated for a conformer
            raise ValueError(f"MMFF properties could not be calculated for conformer {conf_id}")
            # Remove the conformer from the molecule object if its MMFF properties cannot be obtained
            opt_mol.RemoveConformer(conf_id)
        else:
            # Optimize the conformer
            AllChem.MMFFOptimizeMolecule(opt_mol, mmffVariant='MMFF94', confId=conf_id)
            
            # Obtain MMFF parameters for the conformer
            mmff_forcefield = AllChem.MMFFGetMoleculeForceField(opt_mol, mmff_props, confId=conf_id)
            
            # Calculate the energy of conformer using MMFF parameters
            energy: float = mmff_forcefield.CalcEnergy()
            
            # Store the conformer energy along with its index
            mmff_energies[conf_id] = energy

    # Return the molecule object containing optimized geometries and MMFF energies of those geometries
    return opt_mol, mmff_energies


def save_conformers_to_sdf(mol: Chem.Mol, filename: str) -> None:
    """
    Save a molecule with multiple conformers to an SDF file, each starting with the conformer ID.

    Args:
    mol (Chem.Mol): The molecule containing multiple conformers.
    filename (str): The name of the output SDF file.
    """
    writer = SDWriter(filename)
    
    for conf_id in range(mol.GetNumConformers()):
        # Set the conformer ID as a property
        mol.SetProp("_Name", f"conformer_{conf_id}")
        
        # Add the molecule with the current conformer to the writer
        writer.write(mol, confId=conf_id)
    
    writer.close()

    
def calculate_tanimoto_similarity(smiles1: str, smiles2: str):
    """
    Calculate the Tanimoto similarity coefficient between two SMILES strings using RDKit.

    Parameters:
        smiles1 (str): The first SMILES string.
        smiles2 (str): The second SMILES string.

    Returns:
        float: The Tanimoto similarity coefficient.
    """

    # Convert SMILES strings to RDKit Mol objects
    mol_1 = Chem.MolFromSmiles(smiles1)
    mol_2 = Chem.MolFromSmiles(smiles2)
    
    if mol_1 is not None and mol_2 is not None:
        invgen = AllChem.GetMorganFeatureAtomInvGen()
        ffpgen = AllChem.GetMorganGenerator(radius=2, atomInvariantsGenerator=invgen)
        ffp_1 = ffpgen.GetSparseCountFingerprint(mol_1)
        ffp_2 = ffpgen.GetSparseCountFingerprint(mol_2)

        similarity_score = DataStructs.TanimotoSimilarity(ffp_1, ffp_2)
        return similarity_score
    else:
        return None


def convert_conformers_to_smiles(opt_conf_sdf: str, opt_conf_smiles_file: str) -> None:
    """
    Converts 3D geometries of conformers into SMILES and saves them to a file.

    :param opt_conf_sdf: Path to the SDF file containing the optimized conformers.
    :param opt_conf_smiles_file: Path where the SMILES output should be saved.
    """

    # Step 1: Convert SDF to SMILES using obabel
    subprocess.run(f"obabel {opt_conf_sdf} -isdf -osmi -O dummy_SMILES_1.smi", shell=True, check=True)

    # Step 2: Extract the first column from dummy_SMILES_1.smi and save it to dummy_SMILES_2.smi
    subprocess.run("awk '{print $1}' dummy_SMILES_1.smi > dummy_SMILES_2.smi", shell=True, check=True)

    # Step 3: Copy dummy_SMILES_2.smi to the desired output file
    subprocess.run(f"cp -p dummy_SMILES_2.smi {opt_conf_smiles_file}", shell=True, check=True)

    # Step 4: Remove temporary files
    subprocess.run("rm dummy_SMILES_1.smi dummy_SMILES_2.smi", shell=True, check=True)


def calculate_rmsd(input_sdf: str, output_file: str) -> None:
    """
    Run the Open Babel `obrms` command to calculate the RMSD between geometries
    and append the output to the specified file.

    :param input_sdf: Path to the SDF file containing the geometries.
    :param output_file: Path to the file where the RMSD results will be appended.
    """
    command = f"obrms -m -x {input_sdf} >> {output_file}"
    
    # Execute the command
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"RMSD calculation completed and appended to {output_file}.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while calculating RMSD: {e}")


def process_conformers(
    opt_conf_SMILES_file: str,
    opt_conf_sdf: str,
    feasible_geometries_sdf: str,
    infeasible_geometries_sdf: str,
    similarity_output_csv: str,
    infeasible_geometries_csv: str,
    inp_smiles: str,
    num_opt_conf: int,
    energy_DF: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Process optimized conformers to calculate Tanimoto similarity and separate feasible and infeasible geometries.

    :param opt_conf_SMILES_file: File containing SMILES codes for the 3D conformers
    :param opt_conf_sdf: SDF file of the optimized conformers
    :param feasible_geometries_sdf: Output SDF file for feasible geometries
    :param infeasible_geometries_sdf: Output SDF file for infeasible geometries
    :param similarity_output_csv: Output CSV file for Tanimoto similarity results
    :param infeasible_geometries_csv: Output CSV file for infeasible geometries
    :param inp_smiles: Original SMILES of the molecule
    :param num_opt_conf: Number of optimized conformers
    :param energy_DF: DataFrame containing energy values of the conformers
    :return: Tuple of DataFrames (infeasible_geom_DF, energy_DF)
    """

    conf_id: int = 0
    tanimoto_similarity_LD: List[Dict] = []
    infeasible_geom_LD: List[Dict] = []

    # Open the file containing SMILES codes corresponding to all the generated 3D conformers
    with open(opt_conf_SMILES_file, 'r') as SMILES_file:
        opt_conf_SMILES: List[str] = SMILES_file.readlines()
    
    # Read the input SDF file
    suppl = Chem.SDMolSupplier(opt_conf_sdf, removeHs=False)
    
    # Initialize SDWriters for feasible and infeasible geometries
    feasible_writer = Chem.SDWriter(feasible_geometries_sdf)
    infeasible_writer = Chem.SDWriter(infeasible_geometries_sdf)

    # Run a loop over all optimized conformers
    for conf_id in range(num_opt_conf):
        mol = suppl[conf_id]
        # Sanitize the molecule (fix common errors)
        Chem.SanitizeMol(mol)
        curr_SMILES = opt_conf_SMILES[conf_id].strip()
        # Calculate the Tanimoto similarity between SMILES of current conformer and original SMILES
        tanimoto_similarity = calculate_tanimoto_similarity(inp_smiles, curr_SMILES)
        tanimoto_similarity_LD.append({
            'conformer_id': conf_id,
            'SMILES_from_3D_conformer': curr_SMILES,
            'tanimoto_similarity_with_original_SMILES': tanimoto_similarity
        })
        
        # Separate the chemically feasible and infeasible geometries
        if tanimoto_similarity == 1.0:
            # Write the current conformer to feasible geometries SDF file
            mol.SetProp("ID", f"conformer_{conf_id}")
            feasible_writer.write(mol)
        else:
            infeasible_geom_LD.append({
                'conformer_id': conf_id,
                'original_SMILES_of_the_molecule': inp_smiles,
                'SMILES_from_3D_conformer': curr_SMILES,
                'tanimoto_similarity': tanimoto_similarity
            })
            # Write the current conformer to infeasible geometries SDF file
            mol.SetProp("ID", f"conformer_{conf_id}")
            infeasible_writer.write(mol)
            # Drop the row corresponding to infeasible geometry in the energy dataframe
            energy_DF.drop(index=conf_id, inplace=True)

    # Close the SDWriters
    feasible_writer.close()
    infeasible_writer.close()

    # Convert all lists of dictionaries (LD) to dataframes (DF) and dump the latter into CSV files
    tanimoto_similarity_DF = pd.DataFrame(tanimoto_similarity_LD)
    tanimoto_similarity_DF.to_csv(similarity_output_csv, index=False)
    infeasible_geom_DF = pd.DataFrame(infeasible_geom_LD)
    infeasible_geom_DF.to_csv(infeasible_geometries_csv, index=False)
    energy_DF.reset_index(drop=True, inplace=True)
    
    # Return DataFrames
    return infeasible_geom_DF, energy_DF


def calculate_relative_energies(
    energy_DF: pd.DataFrame,
    output_csv: str
) -> pd.DataFrame:
    """
    Calculate the relative energies of conformers and write the results to a CSV file.

    :param energy_DF: DataFrame containing absolute energies of conformers
    :param output_csv: Output CSV file path for the DataFrame with relative energies
    :return: DataFrame containing relative energies in kcal/mol and kJ/mol
    """
    # Take a copy of the dataframe containing absolute energies
    rel_energy_DF = energy_DF.copy()
    
    # Convert the energy column from 'str' type to 'float' type
    rel_energy_DF['energy_in_kcalpermol'] = rel_energy_DF['energy_in_kcalpermol'].astype(float)
    
    # Compute the minimum of the energies
    min_energy = rel_energy_DF['energy_in_kcalpermol'].min()
    
    # Calculate relative energy in kcal/mol
    rel_energy_DF['rel_energy_in_kcalpermol'] = rel_energy_DF['energy_in_kcalpermol'] - min_energy
    
    # Convert relative energy to kJ/mol
    rel_energy_DF['rel_energy_in_kJpermol'] = rel_energy_DF['rel_energy_in_kcalpermol'] * 4.184
    
    # Write the relative energies to a CSV file
    rel_energy_DF.to_csv(output_csv, index=False)
    
    return rel_energy_DF


def write_cluster_representatives(
    opt_conf_sdf: str,
    cluster_reps_dir: str,
    cluster_reps_sdf: str,
    sorted_cluster_reps_DF: pd.DataFrame,
    cluster_reps_DF: pd.DataFrame,
    cluster_rep_prefix: str,
    conf_extension: str
) -> None:
    """
    Write the coordinates of cluster representative conformers to SDF files.

    :param opt_conf_sdf: Path to the SDF file containing optimized conformers
    :param cluster_reps_dir: Directory where cluster representative SDF files will be saved
    :param cluster_reps_sdf: Path to the concatenated SDF file for all cluster representatives
    :param sorted_cluster_reps_DF: DataFrame containing information on sorted cluster representatives
    :param cluster_reps_DF: DataFrame mapping conformers to cluster IDs
    :param cluster_rep_prefix: Prefix for naming individual SDF files for each cluster representative
    :param conf_extension: File extension for the individual SDF files
    :return: None
    """
    
    # Load the molecules (conformers) from the optimized conformer SDF file
    suppl = Chem.SDMolSupplier(opt_conf_sdf, removeHs=False)

    # Check if the cluster_reps_dir exists, delete it
    if os.path.exists(cluster_reps_dir):
        shutil.rmtree(cluster_reps_dir)

    # Create the directory afresh
    os.makedirs(cluster_reps_dir)

    # Create SD writers for the output SDF files
    writer_all_reps = Chem.SDWriter(cluster_reps_sdf)

    # Loop over all optimized conformers to write the coordinates of cluster representative conformers to an SDF file
    for conf_id, mol in enumerate(suppl):
        if mol is None:
            continue

        # Sanitize the molecule (fix common errors)
        Chem.SanitizeMol(mol)

        # If the current conformer is a cluster representative, write its coordinates to the (i) concatenated SDF file and (ii) individual SDF file
        if conf_id in sorted_cluster_reps_DF['conformer_id'].values:
            writer_all_reps.write(mol)

            clust_id = cluster_reps_DF.loc[cluster_reps_DF['conformer_id'] == conf_id, 'cluster_id'].values[0]  # Extract the 'cluster_id' corresponding to the 'conformer_id'
            rep_of_cluster_sdf = f"{cluster_rep_prefix}{clust_id}{conf_extension}"
            f3_path = os.path.join(cluster_reps_dir, rep_of_cluster_sdf)

            # Write individual cluster representative conformer to a separate SDF file
            writer_single_rep = Chem.SDWriter(f3_path)
            writer_single_rep.write(mol)
            writer_single_rep.close()

    # Close the main SD writer
    writer_all_reps.close()

    print("Completed Writing the SDF files of cluster representative conformers")


def calculate_min_rmsd(
    ref_confo_path: str,
    cluster_reps_sdf: str,
    output_sdf: str = "cluster_rep_conformers_vs_ref_conformer.sdf",
    dat_file: str = "rmsd-cluster_rep_conformers_vs_ref_conformer-fm_flags.dat"
) -> Optional[float]:
    """
    Calculate the minimum RMSD between the reference conformer and cluster representatives.

    :param ref_confo_path: Path to the reference conformer file
    :param cluster_reps_sdf: Path to the SDF file containing cluster representative conformers
    :param output_sdf: Path to the output SDF file for RMSD calculations
    :param dat_file: Path to the file where RMSD results are stored
    :return: Minimum RMSD value between the reference conformer and cluster representatives, or None if an error occurs
    """
    try:
        # Run the RMSD calculation command
        subprocess.run(
            f"obrms -f -m {ref_confo_path} {cluster_reps_sdf} -o {output_sdf} >> {dat_file}",
            shell=True,
            check=True
        )

        # Define the shell command to extract the minimum RMSD
        command = f"awk '{{print $NF}}' {dat_file} | sort -n | head -n 1"

        # Run the command and capture the output
        result = subprocess.check_output(command, shell=True).decode('utf-8').strip()

        # Return the minimum RMSD as a float
        return float(result)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running subprocess commands: {e}")
        return None
    except ValueError as e:
        print(f"Error occurred while converting RMSD result: {e}")
        return None

def split_sdf_into_individual_files(input_sdf_path: str, output_dir: str) -> None:
    """
    Splits a multi-conformer SDF file into individual SDF files, one for each conformer.

    Args:
        input_sdf_path (str): Path to the input multi-conformer SDF file.
        output_dir (str): Directory where individual SDF files will be saved.
    """
    # Create the output directory if it doesn't exist, and clear it if it does
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)

    # Read the input SDF file
    suppl = Chem.SDMolSupplier(input_sdf_path, removeHs=False)

    # Iterate through each molecule (conformer) and save it as a separate SDF file
    for i, mol in enumerate(suppl):
        if mol is None:
            continue
        
        # Sanitize the molecule (fix common errors)
        Chem.SanitizeMol(mol)

        # Define the output filename for the individual conformer
        output_filename = os.path.join(output_dir, f"conformer_{i}.sdf")
        
        # Write the individual conformer to a new SDF file
        writer = Chem.SDWriter(output_filename)
        writer.write(mol)
        writer.close()
    print(f"Split {i+1} conformers from {input_sdf_path} into individual SDF files in {output_dir}")
