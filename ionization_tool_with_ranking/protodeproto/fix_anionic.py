
import os
from rdkit import Chem
from rdkit.Chem import AllChem

def fix_anionic(sdf_path: str, output_path: str):
    mol = Chem.MolFromMolFile(sdf_path, removeHs=False)
    if mol is None:
        raise ValueError(f"Could not load molecule from {sdf_path}")

    base = os.path.splitext(os.path.basename(sdf_path))[0]
    if output_path:
        output_path = os.path.join(output_path, f"{base}_protonated.sdf")
    else:
        output_path = f"{base}_protonated.sdf"

    editable = Chem.RWMol(mol)
    treated = False

    for atom in editable.GetAtoms():
        if atom.GetFormalCharge() == -1 and atom.GetAtomicNum() in (7, 8):  # N or O
            atom.SetFormalCharge(0)
            h_idx = editable.AddAtom(Chem.Atom(1))  # Add H atom
            editable.AddBond(atom.GetIdx(), h_idx, Chem.BondType.SINGLE)
            treated = True
            break  # Only one center treated

    if not treated:
        print(f"No anionic center found in {sdf_path}")
        return

    mol_treated = editable.GetMol()
    Chem.SanitizeMol(mol_treated)

    mol_treated = Chem.AddHs(mol_treated)
    AllChem.EmbedMolecule(mol_treated, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol_treated)

    Chem.SDWriter(output_path).write(mol_treated)
    print(f"Saved anion-treated molecule to {output_path}")
