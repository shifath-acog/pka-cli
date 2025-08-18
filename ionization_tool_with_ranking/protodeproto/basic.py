
import os
from rdkit import Chem
from rdkit.Chem import AllChem

def protonate_most_negative_nitrogen(sdf_path: str, charge_path: str, output_path: str = None) -> None:
    """
    Protonates the most negative nitrogen in a molecule using Mulliken charges from a .charge file
    and optimizes geometry with UFF. The updated molecule is saved as an SDF.

    Parameters:
        sdf_path (str): Path to input .sdf file
        charge_path (str): Path to Mulliken .charge file
        output_path (str): Path to save protonated and optimized molecule (default: adds '_base_treated.sdf')
    """
    mol = Chem.MolFromMolFile(sdf_path, removeHs=False)
    if mol is None:
        raise ValueError(f"Failed to load molecule from: {sdf_path}")

    base = os.path.splitext(os.path.basename(sdf_path))[0]
    if output_path:
        output_path = os.path.join(output_path, f"{base}_protonated.sdf")
    else:
        output_path = f"{base}_protonated.sdf"

    # Parse .charge file
    nitrogen_charges = []
    with open(charge_path) as f:
        lines = f.readlines()[1:]  # Skip header
        for line in lines:
            parts = line.split()
            if len(parts) >= 3:
                idx = int(parts[0]) - 1  # Convert to 0-based indexing
                symbol = parts[1]
                charge = float(parts[2])
                if symbol == 'N':
                    nitrogen_charges.append((idx, charge))

    if not nitrogen_charges:
        raise ValueError("No nitrogen atoms found in .charge file")

    # Get the nitrogen with most negative charge
    n_idx, n_charge = min(nitrogen_charges, key=lambda x: x[1])
    print(f"Most negative N: atom index {n_idx}, charge {n_charge:.6f}")

    # Modify molecule: add H to nitrogen
    editable = Chem.RWMol(mol)
    h_idx = editable.AddAtom(Chem.Atom(1))  # Add hydrogen atom
    editable.AddBond(n_idx, h_idx, Chem.rdchem.BondType.SINGLE)
    mol_protonated = editable.GetMol()

    # Set nitrogen's formal charge to +1
    mol_protonated.GetAtomWithIdx(n_idx).SetFormalCharge(1)

    # Sanitize and optimize
    Chem.SanitizeMol(mol_protonated)
    AllChem.EmbedMolecule(mol_protonated, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol_protonated)

    # Save output
    writer = Chem.SDWriter(output_path)
    writer.write(mol_protonated)
    writer.close()
    print(f"Protonated and optimized structure saved to '{output_path}'")
