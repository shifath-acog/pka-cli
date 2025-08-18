# python dihedral_filter.py 'O=C(O)c1ccccc1O' initial_generated_conformers.sdf filtered_isomers_1.sdf

import sys
import shutil
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
import argparse

# Function to check for cis stereochemistry in a SMILES string
def is_cis_smiles(smiles):
    """
    Checks if a SMILES string contains a cis (Z) double bond.

    Args:
        smiles (str): The SMILES string to check.

    Returns:
        bool: True if a cis double bond is found, False otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Warning: Could not parse SMILES: {smiles}", file=sys.stderr)
        return False
    # Iterate through bonds and check for STEREOZ on double bonds
    return any(bond.GetStereo() == Chem.BondStereo.STEREOZ for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)

# Function to check for trans stereochemistry in a SMILES string
def is_trans_smiles(smiles):
    """
    Checks if a SMILES string contains a trans (E) double bond.

    Args:
        smiles (str): The SMILES string to check.

    Returns:
        bool: True if a trans double bond is found, False otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Warning: Could not parse SMILES: {smiles}", file=sys.stderr)
        return False
    # Iterate through bonds and check for STEREOE on double bonds
    return any(bond.GetStereo() == Chem.BondStereo.STEREOE for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)

# Function to filter isomers based on dihedral angles
def filter_isomers_by_dihedral(input_sdf, output_sdf, stereo_type, angle_target, tolerance=40.0):
    """
    Filters conformers in an SDF file based on the dihedral angle of specified stereocenters.

    Args:
        input_sdf (str): Path to the input SDF file.
        output_sdf (str): Path to the output SDF file for filtered conformers.
        stereo_type (Chem.BondStereo): The RDKit stereochemistry type (STEREOZ for cis, STEREOE for trans).
        angle_target (float): The target dihedral angle in degrees (0.0 for cis, 180.0 for trans).
        tolerance (float): The allowed deviation from the target angle in degrees.

    Returns:
        tuple: A tuple containing (count_pass, count_total) of conformers.
    """
    writer = Chem.SDWriter(output_sdf)
    count_total = 0
    count_pass = 0

    # Iterate through molecules in the input SDF file
    for mol in Chem.SDMolSupplier(input_sdf, removeHs=False):
        if mol is None:
            print(f"Warning: Skipping invalid molecule in {input_sdf}", file=sys.stderr)
            continue
        count_total += 1

        # Ensure the molecule has a conformer
        if not mol.GetNumConformers():
            print(f"Warning: Molecule {count_total} has no conformers. Skipping.", file=sys.stderr)
            continue

        conformer = mol.GetConformer()
        all_dihedrals_ok = True
        dihedral_info = []

        # Check all double bonds for the specified stereochemistry type
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetStereo() == stereo_type:
                begin_atom = bond.GetBeginAtom()
                end_atom = bond.GetEndAtom()

                # Get neighbors for dihedral angle calculation
                begin_neighbors = [a.GetIdx() for a in begin_atom.GetNeighbors() if a.GetIdx() != end_atom.GetIdx()]
                end_neighbors = [a.GetIdx() for a in end_atom.GetNeighbors() if a.GetIdx() != begin_atom.GetIdx()]

                # Need at least one neighbor on each side to define a dihedral
                if not begin_neighbors or not end_neighbors:
                    print(f"Warning: Skipping bond between {begin_atom.GetIdx()}-{end_atom.GetIdx()} due to insufficient neighbors for dihedral.", file=sys.stderr)
                    continue

                # Define the four atoms for dihedral calculation
                atom1 = begin_neighbors[0]
                atom2 = begin_atom.GetIdx()
                atom3 = end_atom.GetIdx()
                atom4 = end_neighbors[0]

                # Calculate the dihedral angle
                angle = rdMolTransforms.GetDihedralDeg(conformer, atom1, atom2, atom3, atom4)
                dihedral_info.append((atom1, atom2, atom3, atom4, angle))

                # Check if the angle is within tolerance
                # Using abs(abs(angle) - angle_target) handles both positive and negative angles
                if abs(abs(angle) - angle_target) > tolerance:
                    all_dihedrals_ok = False
                    break # If one dihedral is out of range, no need to check others for this conformer

        # If all relevant dihedrals are within tolerance and at least one was found
        if all_dihedrals_ok and dihedral_info:
            for i, (a1, a2, a3, a4, angle) in enumerate(dihedral_info):
                mol.SetProp(f"Dihedral_{i+1}", f"{a1}-{a2}-{a3}-{a4} : {angle:.2f}°")
            writer.write(mol)
            count_pass += 1

    writer.close()
    return count_pass, count_total

# Main execution block
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter conformers in an SDF file based on cis/trans stereochemistry and dihedral angles."
    )
    parser.add_argument(
        "input_smiles",
        type=str,
        help="Input SMILES string to determine cis/trans stereochemistry (e.g., 'O=C(O)/C=C/C(=O)O')."
    )
    parser.add_argument(
        "input_sdf",
        type=str,
        help="Path to the input SDF file containing conformers."
    )
    parser.add_argument(
        "output_sdf",
        type=str,
        help="Path to the output SDF file for filtered conformers."
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=40.0,
        help="Tolerance in degrees for dihedral angle filtering (default: 40.0)."
    )

    args = parser.parse_args()

    input_smiles = args.input_smiles
    input_sdf = args.input_sdf
    output_sdf = args.output_sdf
    tolerance = args.tolerance

    mol = Chem.MolFromSmiles(input_smiles)
    if mol is None:
        print(f"Error: Could not parse input SMILES '{input_smiles}'. Exiting.", file=sys.stderr)
        sys.exit(1)

    canon_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)

    if is_cis_smiles(canon_smiles):
        stereo = Chem.BondStereo.STEREOZ
        angle_target = 0.0
        print(f"Cis isomer detected for '{input_smiles}' → filtering around 0° with tolerance {tolerance}°")
        matched, total = filter_isomers_by_dihedral(input_sdf, output_sdf, stereo, angle_target, tolerance=tolerance)
        print(f"Matched {matched} / {total} conformers → '{output_sdf}'")
    elif is_trans_smiles(canon_smiles):
        stereo = Chem.BondStereo.STEREOE
        angle_target = 180.0
        print(f"Trans isomer detected for '{input_smiles}' → filtering around ±180° with tolerance {tolerance}°")
        matched, total = filter_isomers_by_dihedral(input_sdf, output_sdf, stereo, angle_target, tolerance=tolerance)
        print(f"Matched {matched} / {total} conformers → '{output_sdf}'")
    else:
        print(f"Input SMILES '{input_smiles}' does not contain cis/trans stereochemistry. Copying input SDF to output SDF.")
        try:
            shutil.copy(input_sdf, output_sdf)
            print(f"Copied '{input_sdf}' to '{output_sdf}'.")
        except FileNotFoundError:
            print(f"Error: Input SDF file '{input_sdf}' not found. Please ensure it exists.", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(f"An error occurred while copying the file: {e}", file=sys.stderr)
            sys.exit(1)

