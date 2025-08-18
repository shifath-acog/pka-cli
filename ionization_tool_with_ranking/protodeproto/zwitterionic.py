
import os
from rdkit import Chem
from rdkit.Chem import AllChem


def load_mol_and_charges(sdf_path: str, charge_path: str):
    mol = Chem.MolFromMolFile(sdf_path, removeHs=False)
    if mol is None:
        raise ValueError(f"Could not load molecule from {sdf_path}")

    atom_charges = {}
    with open(charge_path) as f:
        lines = f.readlines()[1:]  # skip header
        for line in lines:
            parts = line.split()
            if len(parts) >= 3:
                idx = int(parts[0]) - 1
                charge = float(parts[2])
                atom_charges[idx] = charge
    return mol, atom_charges

def deprotonate_positive_n_with_h(mol: Chem.Mol) -> Chem.Mol:
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N' and atom.GetFormalCharge() == 1:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'H':
                    editable = Chem.RWMol(mol)
                    atom_idx = atom.GetIdx()
                    h_idx = neighbor.GetIdx()

                    editable.RemoveBond(atom_idx, h_idx)
                    editable.RemoveAtom(h_idx)
                    editable.GetAtomWithIdx(atom_idx).SetFormalCharge(0)

                    mol_new = editable.GetMol()
                    Chem.SanitizeMol(mol_new)
                    return mol_new
    raise ValueError("No protonated N with a bound hydrogen found.")


def protonate_most_negative_anion(mol: Chem.Mol, charges: dict) -> Chem.Mol:
    candidates = []
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() == -1 and atom.GetSymbol() in ['O', 'S', 'N']:
            candidates.append((atom.GetIdx(), charges.get(atom.GetIdx(), 0.0)))

    if not candidates:
        raise ValueError("No negatively charged atoms found to protonate.")

    candidates.sort(key=lambda x: x[1])  # most negative first
    target_idx = candidates[0][0]

    editable = Chem.EditableMol(mol)
    h_idx = editable.AddAtom(Chem.Atom(1))  # Hydrogen
    editable.AddBond(target_idx, h_idx, Chem.rdchem.BondType.SINGLE)

    mol_new = editable.GetMol()
    mol_new.GetAtomWithIdx(target_idx).SetFormalCharge(0)
    Chem.SanitizeMol(mol_new)
    AllChem.EmbedMolecule(mol_new, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol_new)
    return mol_new


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
        if atom.GetSymbol() in ('O', 'S'):
            idx = atom.GetIdx()
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'H':
                    charge = atom_charges.get(idx, 0.0)
                    candidates.append((idx, charge, atom.GetSymbol()))
                    break

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

def handle_zwitterionic(sdf_path, charge_path, output_dir):
    mol, charges = load_mol_and_charges(sdf_path, charge_path)

    base = os.path.splitext(os.path.basename(sdf_path))[0]
    base_file = os.path.join(output_dir, f"{base}_deprotonated.sdf")
    acid_file = os.path.join(output_dir, f"{base}_protonated.sdf")

    # Try deprotonate positively charged nitrogen, or fallback to OH/SH deprotonation
    try:
        mol_base = deprotonate_positive_n_with_h(mol)
        writer = Chem.SDWriter(base_file)
        writer.write(mol_base)
        writer.close()
    except ValueError as e:
        print(f"Warning: {e} — attempting OH/SH deprotonation instead.")
        try:
            deprotonate_most_positive_oh_or_sh(sdf_path, charge_path, output_dir)
        except ValueError as e2:
            print(f"Warning: {e2} — no suitable deprotonation site found.")

    # Always attempt protonate negatively charged site
    mol_acid = protonate_most_negative_anion(mol, charges)
    writer = Chem.SDWriter(acid_file)
    writer.write(mol_acid)
    writer.close()

