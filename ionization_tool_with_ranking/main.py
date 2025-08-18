import os
import sys
from rdkit import Chem
from protodeproto.acidic import deprotonate_most_positive_oh_or_sh
from protodeproto.deprotonate_basic import deprotonate_most_positive_nh
from protodeproto.basic import protonate_most_negative_nitrogen
from protodeproto.ambiphilic import handle_ambiphilic
from protodeproto.zwitterionic import handle_zwitterionic  # NEW IMPORT
from protodeproto.protonated import handle_protonated
from protodeproto.fix_anionic import fix_anionic

def extract_charge_info(sdf_path):
    """Extracts formal charges from M  CHG lines in a .sdf file"""
    with open(sdf_path) as f:
        for line in f:
            if line.startswith("M  CHG"):
                tokens = line.strip().split()
                n_entries = int(tokens[2])
                charges = {}
                for i in range(n_entries):
                    atom_idx = int(tokens[3 + i * 2])
                    charge = int(tokens[4 + i * 2])
                    charges[atom_idx] = charge
                return charges
    return {}


def classify_by_charge(charges, sdf_path):
    if charges:
        pos = [i for i, c in charges.items() if c > 0]
        neg = [i for i, c in charges.items() if c < 0]
        overall_charge = sum(charges.values())

        if overall_charge == 1:
            return "protonated"
        elif overall_charge == -1:
            return "anionic"

        # Load molecule once if needed
        mol = Chem.MolFromMolFile(sdf_path, sanitize=True)
        if mol is None:
            print(f"Could not parse molecule from {sdf_path}")
            return None

        # SMARTS for NO2 group
        no2_smarts = Chem.MolFromSmarts("[N+](=O)[O-]")
        no2_matches = mol.GetSubstructMatches(no2_smarts)
        no2_atoms = set()
        for match in no2_matches:
            no2_atoms.update(match)  # RDKit atom indices (0-based)

        # Convert CHG indices (1-based) to 0-based for comparison
        filtered_pos = [i - 1 for i in pos]
        filtered_neg = [i - 1 for i in neg]

        # Filter out NO2 atoms
        filtered_pos = [i for i in filtered_pos if i not in no2_atoms]
        filtered_neg = [i for i in filtered_neg if i not in no2_atoms]

        if filtered_pos and filtered_neg:
            return "zwitterionic"
        elif filtered_pos:
            return "protonated"
        elif filtered_neg:
            return "anionic"
        else:
            # All charges came from NO2 — fall back to structural logic
            pass

    # No CHG section or only NO2 charges — fall back logic
    with open(sdf_path) as f:
        lines = f.readlines()

    atom_count = int(lines[3].split()[0])
    bond_count = int(lines[3].split()[1])
    atom_lines = lines[4:4 + atom_count]
    bond_lines = lines[4 + atom_count:4 + atom_count + bond_count]

    elements = {}
    adjacency = {}

    for idx, line in enumerate(atom_lines):
        tokens = line.split()
        elements[idx + 1] = tokens[3]
        adjacency[idx + 1] = []

    for line in bond_lines:
        a1 = int(line[0:3])
        a2 = int(line[3:6])
        adjacency[a1].append(a2)
        adjacency[a2].append(a1)

    has_acidic_group = False
    has_acidic_amino_group = False
    has_basic_group = False

    for atom_idx, element in elements.items():
        if element == "N":
            has_H = any(elements[neighbor] == "H" for neighbor in adjacency[atom_idx])
            if has_H:
                has_acidic_amino_group = True
            else:
                has_basic_group = True
        elif element in ("O", "S"):
            for neighbor in adjacency[atom_idx]:
                if elements[neighbor] == "H":
                    has_acidic_group = True
                    break

    #if has_acidic_group and has_basic_group:
    if (has_acidic_group or has_acidic_amino_group) and has_basic_group:
        return "ambiphilic"
    elif has_acidic_group:
        return "acidic"
    elif has_acidic_amino_group:
        return "deprotonate_basic"
    elif has_basic_group:
        return "basic"
    else:
        return None

def process_file(sdf_path, charge_path, output_dir):  # Edit 2
    base = os.path.splitext(os.path.basename(sdf_path))[0]
    charges = extract_charge_info(sdf_path)
    classification = classify_by_charge(charges, sdf_path)

    if classification is None:
        print(f"No titratable groups found in {base}.sdf")
        return

    print(f"{base}.sdf classified as {classification}")
    print(f"Classification for {base}.sdf: {classification}") # Debug print

    if classification == "acidic":
        deprotonate_most_positive_oh_or_sh(sdf_path, charge_path, output_dir)  # Edit 3
    elif classification == "deprotonate_basic":
        deprotonate_most_positive_nh(sdf_path, charge_path, output_dir)  # Edit 3
    elif classification == "basic":
        protonate_most_negative_nitrogen(sdf_path, charge_path, output_dir)  # Edit 3
    elif classification == "ambiphilic":
        handle_ambiphilic(sdf_path, charge_path, output_dir)  # Edit 3
    elif classification == "zwitterionic":
        handle_zwitterionic(sdf_path, charge_path, output_dir)  # Edit 3
    elif classification == "protonated":
        handle_protonated(sdf_path, charge_path, output_dir)  # Edit 3
    elif classification == "anionic":
        fix_anionic(sdf_path, output_dir)

def main():
    if len(sys.argv) != 2:
        print("Usage: python main.py input_folder/")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_dir = os.path.join(input_dir, "treated")
    os.makedirs(output_dir, exist_ok=True)

    for file in os.listdir(input_dir):
        if not file.endswith(".sdf") or file.endswith("_treated.sdf") or "_aug-cc-pVDZ_M06-2X.sdf" in file:
            continue

        sdf_path = os.path.join(input_dir, file)
        charge_path = os.path.join(input_dir, file.replace(".sdf", "_aug-cc-pVDZ_M06-2X.charge"))

        try:
            with open(charge_path, 'r') as f:
                # File exists and is readable
                pass
        except FileNotFoundError:
            print(f"Missing charge file for {file} (FileNotFoundError), skipping.")
            continue
        except Exception as e:
            print(f"Error opening charge file {file}: {e}, skipping.")
            continue

        process_file(sdf_path, charge_path, output_dir)

if __name__ == "__main__":
    main()