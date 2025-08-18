
import os
from rdkit import Chem


def already_protonated(sdf_path, output_path):
    suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
    mol = suppl[0]
    if mol is None:
        raise ValueError(f"Could not read molecule from {sdf_path}")

    editable = Chem.RWMol(mol)
    modified = False
    atoms_to_remove = []

    # Flag to stop after removing one proton
    proton_removed = False

    for atom in editable.GetAtoms():
        if proton_removed:
            break  # Stop if already removed one proton

        if atom.GetSymbol() in ['O', 'N']:
            valence = sum([b.GetBondTypeAsDouble() for b in atom.GetBonds()])
            expected = 2 if atom.GetSymbol() == 'O' else 3
            if valence > expected:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'H':
                        h_idx = neighbor.GetIdx()
                        editable.RemoveBond(atom.GetIdx(), h_idx)
                        atoms_to_remove.append(h_idx)
                        modified = True
                        proton_removed = True
                        break  # stop checking neighbors once proton removed

    for h_idx in sorted(atoms_to_remove, reverse=True):
        editable.RemoveAtom(h_idx)

    if not modified:
        raise ValueError("No overprotonated O or N atoms found.")

    mol_fixed = editable.GetMol()

    # Clear aromaticity flags to avoid kekulization problems (optional)
    for atom in mol_fixed.GetAtoms():
        atom.SetIsAromatic(False)
    for bond in mol_fixed.GetBonds():
        bond.SetIsAromatic(False)

    # Run UFF optimization
    # try:
    #     AllChem.UFFOptimizeMolecule(mol_fixed)
    # except:
    #     print("Warning: UFF optimization failed.")

    writer = Chem.SDWriter(output_path)
    writer.write(mol_fixed)
    writer.close()


def deprotonate_most_positive_oh_or_sh(sdf_path: str, charge_path: str, output_path: str = None) -> None:
    mol = Chem.MolFromMolFile(sdf_path, removeHs=False)
    if mol is None:
        raise ValueError(f"Could not load molecule from {sdf_path}")

    base = os.path.splitext(os.path.basename(sdf_path))[0]

    # If output_path is None, create filename in current dir
    if output_path is None:
        output_path = f"{base}_deprotonated.sdf"
    else:
        # If output_path is a directory, join with filename
        if os.path.isdir(output_path):
            output_path = os.path.join(output_path, f"{base}_deprotonated.sdf")
        else:
            # Otherwise assume output_path is full filename, use as is
            pass        
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

    # print(f"Deprotonated molecule saved as '{output_path}'")

def handle_protonated(sdf_path: str, charge_path: str, output_dir: str) -> None:
    import os
    from rdkit import Chem

    base = os.path.splitext(os.path.basename(sdf_path))[0]
    output_file = os.path.join(output_dir, f"{base}_deprotonated.sdf")

    # Try to remove an excess proton from overprotonated N/O
    try:
        already_protonated(sdf_path, output_file)
        print(f"Deprotonated overprotonated atom. Saved as '{output_file}'")
        return
    except ValueError as e:
        if "No overprotonated O or N atoms found." in str(e):
            print("No overprotonated N/O atoms found, trying OH/SH instead...")
        else:
            raise

    # If that fails, try to deprotonate OH or SH
    try:
        deprotonate_most_positive_oh_or_sh(sdf_path, charge_path, output_file)
        print(f"Deprotonated OH/SH. Saved as '{output_file}'")
    except ValueError as e:
        print(f"Deprotonation failed: {e}")
