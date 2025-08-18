
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Dict


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


def deprotonate_most_positive_oh_sh_or_nh(mol: Chem.Mol, charges: Dict[int, float]) -> Chem.Mol:
    candidates = []
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol in ['O', 'S', 'N']:  # Include N now
            idx = atom.GetIdx()
            has_h = any(neighbor.GetSymbol() == 'H' for neighbor in atom.GetNeighbors())
            if has_h:
                candidates.append((idx, charges.get(idx, 0.0), symbol))

    if not candidates:
        raise ValueError("No O-H, S-H or N-H groups found.")

    # Sort by most positive charge to pick the most acidic site
    candidates.sort(key=lambda x: x[1], reverse=True)
    target_idx, target_charge, target_symbol = candidates[0]

    print(f"Deprotonating {target_symbol}-H site at atom index {target_idx + 1} with charge {target_charge:.6f}")

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

    mol_deprot = editable.GetMol()
    mol_deprot.GetAtomWithIdx(target_idx).SetFormalCharge(-1)
    Chem.SanitizeMol(mol_deprot)
    return mol_deprot


def protonate_most_negative_nitrogen(mol: Chem.Mol, charges: Dict[int, float]) -> Chem.Mol:
    candidates = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N' and atom.GetFormalCharge() == 0:
            # Exclude nitrogens that have at least one hydrogen attached
            has_h = any(neighbor.GetSymbol() == 'H' for neighbor in atom.GetNeighbors())
            if not has_h:
                idx = atom.GetIdx()
                candidates.append((idx, charges.get(idx, 0.0)))

    if not candidates:
        raise ValueError("No neutral N atoms without attached H found.")

    candidates.sort(key=lambda x: x[1])  # most negative charge first
    target_idx = candidates[0][0]

    editable = Chem.EditableMol(mol)
    h_idx = editable.AddAtom(Chem.Atom(1))  # Add hydrogen atom
    editable.AddBond(target_idx, h_idx, Chem.rdchem.BondType.SINGLE)

    mol_protonated = editable.GetMol()
    mol_protonated.GetAtomWithIdx(target_idx).SetFormalCharge(1)
    Chem.SanitizeMol(mol_protonated)
    AllChem.EmbedMolecule(mol_protonated, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol_protonated)
    return mol_protonated


def handle_ambiphilic(sdf_path, charge_path, output_dir):
    mol, charges = load_mol_and_charges(sdf_path, charge_path)

    base = os.path.splitext(os.path.basename(sdf_path))[0]
    acid_file = os.path.join(output_dir, f"{base}_deprotonated.sdf")
    base_file = os.path.join(output_dir, f"{base}_protonated.sdf")

    # Use the updated deprotonation function
    mol_acid = deprotonate_most_positive_oh_sh_or_nh(mol, charges)
    writer = Chem.SDWriter(acid_file)
    writer.write(mol_acid)
    writer.close()

    # Protonation logic stays the same
    neutral_nitrogens = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == "N" and atom.GetFormalCharge() == 0]
    if neutral_nitrogens:
        mol_base = protonate_most_negative_nitrogen(mol, charges)
        writer = Chem.SDWriter(base_file)
        writer.write(mol_base)
        writer.close()
    else:
        print("No neutral N atoms found, skipping base protonation.")
