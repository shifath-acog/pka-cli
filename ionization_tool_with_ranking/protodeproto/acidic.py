
import os
from rdkit import Chem
from rdkit.Chem import AllChem

def deprotonate_most_positive_oh_or_sh(sdf_path: str, charge_path: str, output_path: str = None) -> None:
    """
    Deprotonates the most positively charged OH or SH group using Mulliken charges and writes the updated molecule.

    Parameters:
        sdf_path (str): Path to input .sdf file with explicit Hs
        charge_path (str): Path to Mulliken .charge file
        output_path (str): Path to save deprotonated molecule
    """
    mol = Chem.MolFromMolFile(sdf_path, removeHs=False)
    if mol is None:
        raise ValueError(f"Could not load molecule from {sdf_path}")

    base = os.path.splitext(os.path.basename(sdf_path))[0]
    if output_path:
        output_path = os.path.join(output_path, f"{base}_deprotonated.sdf")
    else:
        output_path = f"{base}_deprotonated.sdf"

    # Parse Mulliken charges
    atom_charges = {}
    with open(charge_path) as f:
        lines = f.readlines()[1:]  # Skip header
        for line in lines:
            parts = line.split()
            if len(parts) >= 3:
                idx = int(parts[0]) - 1  # 0-based indexing
                atom_charges[idx] = float(parts[2])

    # Identify O-H and S-H atoms with their charges
    candidates = []
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol in ('O', 'S'):
            # Check if this atom has at least one H neighbor
            has_hydrogen = any(n.GetSymbol() == 'H' for n in atom.GetNeighbors())
            if has_hydrogen:
                idx = atom.GetIdx()
                charge = atom_charges.get(idx, 0.0)
                candidates.append((idx, charge, symbol))


    if not candidates:
        raise ValueError("No O–H or S–H groups found in molecule.")

    # Pick the most positively charged O or S with at least one H
    candidates.sort(key=lambda x: x[1], reverse=True)
    target_idx, target_charge, target_symbol = candidates[0]
    print(f"Deprotonating {target_symbol}-H site at atom index {target_idx + 1} with Mulliken charge {target_charge:.6f}")

    # Remove the bonded H
    editable = Chem.RWMol(mol)
    h_to_remove = None
    for neighbor in mol.GetAtomWithIdx(target_idx).GetNeighbors():
        if neighbor.GetSymbol() == 'H':
            h_to_remove = neighbor.GetIdx()
            break

    if h_to_remove is not None:
        editable.RemoveBond(target_idx, h_to_remove)
        editable.RemoveAtom(h_to_remove)
    else:
        raise ValueError(f"Could not find hydrogen to remove from {target_symbol}.")

    # Finalize and save the molecule
    mol_updated = editable.GetMol()
    mol_updated.GetAtomWithIdx(target_idx).SetFormalCharge(-1) # assign charge
    Chem.SanitizeMol(mol_updated)

    writer = Chem.SDWriter(output_path)
    writer.write(mol_updated)
    writer.close()

    print(f"Deprotonated molecule saved as '{output_path}'")
