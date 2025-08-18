import pickle
import rdkit
from rdkit.Chem.rdchem import Mol
from rdkit.Chem import rdDepictor
from rdkit.Chem import AllChem
from rdkit import Chem
import csv
import sys

pkl_file = sys.argv[1]
#csv_smiles = sys.argv[2]
output_file = sys.argv[2]
objects = []
with (open(pkl_file, "rb")) as openfile:
    while True:
        try:
            objects.append(pickle.load(openfile))
        except EOFError:
            break
conformers_generated = []
for mol in objects:
    for key in mol:
        conformers_generated.append([key, mol[key]])
'''
with open(csv_smiles, "r") as f_cap:
    reader = csv.reader(f_cap)
    data = [tuple(row) for row in reader]

dict_molecules = {}
for x in data:
    dict_molecules[x[0]] = x[2]
'''
num_smiles = 0
for x in conformers_generated:
    output_file_str = output_file + \
        f"_{num_smiles}_conformers.sdf"
    
    # Create SDF writer
    writer = Chem.SDWriter(output_file_str)
    
    attribute_file_str = output_file + \
        f"_{num_smiles}_attributes.csv"
    with open(attribute_file_str, "w") as fa:
        fa.write("dlogp,euclidean_dlogp,mmff_energy,rmsd \n")
    
    counter = 1
    for molecule in x[1]:
        # Write attributes to CSV
        with open(attribute_file_str, "a") as fa:
            fa.write(f'{molecule.dlogp},{molecule.euclidean_dlogp},{molecule.mmff_energy},{molecule.rmsd}\n')
        
        # Clean and prepare the molecule for SDF output
        try:
            # Sanitize the molecule to ensure proper structure
            Chem.SanitizeMol(molecule)
            
            # Set molecule properties as SDF data fields
            molecule.SetProp("_Name", f"conformer_{counter}")
            molecule.SetProp("dlogp", str(float(molecule.dlogp)))
            molecule.SetProp("euclidean_dlogp", str(float(molecule.euclidean_dlogp)))
            molecule.SetProp("mmff_energy", str(float(molecule.mmff_energy)))
            molecule.SetProp("rmsd", str(float(molecule.rmsd)))
            molecule.SetProp("conformer_id", str(counter))
            
            # Write molecule to SDF
            writer.write(molecule)
            
        except Exception as e:
            print(f"Error processing molecule {counter}: {e}")
            # Try to write without sanitization if sanitization fails
            try:
                molecule.SetProp("_Name", f"conformer_{counter}")
                molecule.SetProp("dlogp", str(float(molecule.dlogp)))
                molecule.SetProp("euclidean_dlogp", str(float(molecule.euclidean_dlogp)))
                molecule.SetProp("mmff_energy", str(float(molecule.mmff_energy)))
                molecule.SetProp("rmsd", str(float(molecule.rmsd)))
                molecule.SetProp("conformer_id", str(counter))
                writer.write(molecule)
            except Exception as e2:
                print(f"Failed to write molecule {counter}: {e2}")
                continue
        
        counter += 1
    
    # Close the SDF writer
    writer.close()
    num_smiles += 1