from rdkit import Chem
from pyscf import gto, dft
from gpu4pyscf.dft import rks
from pyscf.geomopt import geometric_solver
import sys
import glob
import time
import re


import os
# Specify which GPUs to use
os.environ["CUDA_VISIBLE_DEVICES"] = "0,1"
os.environ["OMP_NUM_THREADS"] = "25"
os.environ["MKL_NUM_THREADS"] = "25"
os.environ["OPENBLAS_NUM_THREADS"] = "25"


class DFTClass:
    def __init__(self, sdf_path, eps):
        self.sdf_path= sdf_path
        self.eps=eps
    
    def get_atom_list(self):
        molecules = Chem.SDMolSupplier(self.sdf_path, removeHs=False)
        mol = next(molecules)  # Assuming you have one molecule in the SDF

        # Get the first conformer to access 3D coordinates
        conformer = mol.GetConformer()

        # Extract atomic symbols and coordinates
        atom_list = []
        for atom in mol.GetAtoms():
            pos = conformer.GetAtomPosition(atom.GetIdx())
            atom_list.append((atom.GetSymbol(), (pos.x, pos.y, pos.z)))
        return atom_list
       
    def extract_charge_info(self):
        """
        Extracts total formal charge from M  CHG lines in a .sdf file.
        If no M  CHG line is found, assumes formal charge = 0.
        """
        charges = {}
        with open(self.sdf_path) as f:
            for line in f:
              if line.startswith("M  CHG"):
                 tokens = line.strip().split()
                 n_entries = int(tokens[2])
                 for i in range(n_entries):
                    atom_idx = int(tokens[3 + i * 2])
                    charge = int(tokens[4 + i * 2])
                    charges[atom_idx] = charge

        if charges:
            formal_charge = sum(charges.values())
        else:
            formal_charge = 0  # Default if no charge info present

        #print(f"{os.path.basename(sdf_path)}: formal charge = {formal_charge}")
        return formal_charge  
    
    def get_charge(self, cluster_reps_dir):
        """
        Reads the number of atoms directly from an SDF file.

        Args:
            file_path: The path to the SDF file.

        Returns:
            An integer representing the number of atoms, or None if the
            information cannot be extracted or the file cannot be read.
        """
        def extract_charge_info_static(path):
            """
            Extracts total formal charge from M  CHG lines in a .sdf file.
            If no M  CHG line is found, assumes formal charge = 0.
            """
            charges = {}
            with open(path) as f:
                for line in f:
                    if line.startswith("M  CHG"):
                        tokens = line.strip().split()
                        n_entries = int(tokens[2])
                        for i in range(n_entries):
                            atom_idx = int(tokens[3 + i * 2])
                            charge = int(tokens[4 + i * 2])
                            charges[atom_idx] = charge

            if charges:
                formal_charge = sum(charges.values())
            else:
                formal_charge = 0  # Default if no charge info present

            #print(f"{os.path.basename(sdf_path)}: formal charge = {formal_charge}")
            return formal_charge 


        def get_num_atoms(path):
            try:
                with open(path, 'r') as f:
                    # Skip the first three lines of the SDF file
                    for _ in range(3):
                        f.readline()
                    # Read the fourth line which contains the number of atoms
                    fourth_line = f.readline().strip()
                    parts = fourth_line.split()
                    if parts:
                        try:
                            num_atoms = int(parts[0])
                            return num_atoms
                        except ValueError:
                            print("Error: Could not convert atom count to integer.")
                            return None
                    else:
                        print("Error: Fourth line of SDF is empty.")
                        return None
            except FileNotFoundError:
                print(f"Error: File not found at {self.sdf_path}")
                return None
            except Exception as e:
                print(f"An error occurred while reading the file: {e}")
                return None
            
        #Get the charge
        A_atoms = get_num_atoms(self.sdf_path)
        clean_path = re.sub(r'(_protonated|_deprotonated)', '', self.sdf_path)
        HA_atoms = get_num_atoms(clean_path)
        formal_charge_on_HA=extract_charge_info_static(os.path.join(cluster_reps_dir, "rep_of_cluster_1.sdf"))
        # Check for invalid (None) atom counts
        charge=formal_charge_on_HA + (A_atoms-HA_atoms)
        print(f"Charge of {self.sdf_path}: {charge}")
        return charge
     
    def opti_PCM(self, mol, xyz_filename, charge_filename, xc):
        start_time = time.time()
        # Set up the DFT calculation
        mf = rks.RKS(mol).density_fit()  # Use density fitting for efficiency
        mf.xc = xc  # Set the exchange-correlation functional

        # SCF convergence Criteria
        mf.conv_tol = 1e-8  # Energy convergence
        mf.conv_tol_grad = 3e-4  # Gradient convergence
        mf.max_cycle = 70  # Increase max iterations if needed

        # Apply the solvation model
        mf = mf.PCM()  # Initialize solvation model
        mf.grids.atom_grid = (99, 590)
        mf.with_solvent.lebedev_order = 29  # 302 Lebedev grids
        mf.with_solvent.method = 'IEF-PCM'  # Can be C-PCM, SS(V)PE, COSMO
        mf.with_solvent.eps = self.eps  # Set the solvent's dielectric constant

        # Perform geometry optimization
        print("Starting geometry optimization...")
        mol_opt = geometric_solver.optimize(mf, max_steps=200, xtol=1e-8, gtol=3e-4, etol=1e-8)

        # Output optimized geometry
        optimized_atoms = [(atom[0], mol_opt.atom_coords(unit='Angstrom')[i]) for i, atom in enumerate(mol_opt.atom)]

        # Ensure SCF calculation is performed after optimization
        mf_scf = rks.RKS(mol_opt).density_fit()
        mf_scf.xc = xc

        # SCF convergence Criteria
        mf_scf.conv_tol = 1e-8  # Energy convergence
        mf_scf.conv_tol_grad = 3e-4  # Gradient convergence
        mf_scf.max_cycle = 70  # Increase max iterations if needed

        # Apply the solvation model
        mf_scf = mf_scf.PCM()  # Initialize solvation model
        mf_scf.grids.atom_grid = (99, 590)
        mf_scf.with_solvent.lebedev_order = 29  # 302 Lebedev grids
        mf_scf.with_solvent.method = 'IEF-PCM'  # Can be C-PCM, SS(V)PE, COSMO
        mf_scf.with_solvent.eps = self.eps  # Set the solvent's dielectric constant

        #Run the scf
        mf_scf.kernel()

        #Mulliken Charge Analysis
        analysis = mf_scf.analyze()
        mulliken_charges = analysis[0][1]  # Get the Mulliken charges

        # Save captured output to file
        with open(charge_filename, "w") as charge_file:
            charge_file.write("Atom Index  Atom Symbol  Mulliken Charge\n")
            for i, charge in enumerate(mulliken_charges):
                atom_symbol = mol.atom_symbol(i)  # Get the atom symbol
                charge_file.write(f"{i+1}             {atom_symbol}        {charge:.6f}\n")


        # Extract the final energy in Hartree
        final_energy_hartree = mf_scf.e_tot

        # Convert energy from Hartree to kJ/mol
        #hartree_to_kjmol = 2625.5
        hartree_to_kcalmol=627.509
        final_energy_kcalmol = final_energy_hartree * hartree_to_kcalmol

        # Save optimized geometry to XYZ file
        with open(xyz_filename, 'w') as xyz_file:
            xyz_file.write(f"{len(optimized_atoms)}\n")
            xyz_file.write(f"Energy: {final_energy_kcalmol:.6f} kcal/mol\n")
            for symbol, coords in optimized_atoms:
                formatted_coords = ' '.join(f"{coord:.8f}" for coord in coords)
                xyz_file.write(f"{symbol} {formatted_coords}\n")

        print(f"Optimized geometry saved to '{xyz_filename}'.")

        # Print the final energy
        print(f"Final energy: {final_energy_hartree:.8f} Hartree ({final_energy_kcalmol:.6f} kcal/mol)")

        # Record the end time
        opt_time = time.time()

        # Calculate and print the total run time
        total_opt_time = opt_time - start_time
        print(f"\nOPT Time: {total_opt_time:.2f} seconds")

        print("################################################################")
        return mol_opt
    

def run_DFT(cluster_reps_dir, dielectric_value, mode='normal'):
    """
    Main function that processes SDF files with specified basis sets and XC functionals.
    Accepts command-line argument for dielectric constant.
    """
    
    
    # Define basis sets and XC functionals to loop over
    basis_sets = ['aug-cc-pVDZ']
    xc_functionals = ['M06-2X']  # Add more if needed

    # Get list of charged SDF files
    #sdf_files = sorted(glob.glob(os.path.join(cluster_reps_dir, file_pattern)))
    if mode == 'deproto':
        sdf_files = sorted(glob.glob(os.path.join(cluster_reps_dir, "*nated.sdf")))
    else:
        sdf_files = sorted(glob.glob(os.path.join(cluster_reps_dir, "*.sdf")))

    if not sdf_files:
        print(f"No SDF files found in {cluster_reps_dir}")
        return

    print(f"Found {len(sdf_files)} SDF files in {cluster_reps_dir}")
    print(f"Basis sets: {basis_sets}")
    print(f"XC functionals: {xc_functionals}")
    print(f"Dielectric constant: {dielectric_value}")
    print("=" * 70)

    # Loop over each SDF file
    for sdf_path in sdf_files:
        dft_obj = DFTClass(sdf_path, dielectric_value)
        atom_list = dft_obj.get_atom_list()

        for basis in basis_sets:
            for xc in xc_functionals:
                start_time_sub = time.time()

                # Use the chosen charge method
                if mode == 'deproto':
                    charge_dft = dft_obj.get_charge(cluster_reps_dir)
                else:
                    charge_dft = dft_obj.extract_charge_info()

                # Build molecule with given basis
                mol1 = gto.M(
                    atom=atom_list,
                    basis=basis,
                    charge =charge_dft,
                    spin=0,
                    verbose=4,
                )

                # Sanitize names for filename safety
                basis_tag = basis.replace("(", "").replace(")", "").replace("*", "").replace("+", "").replace("/", "-").replace(" ", "")
                xc_tag = xc.replace("(", "").replace(")", "").replace("*", "").replace("+", "").replace("/", "-").replace(" ", "")

                # Output file names
                if mode == 'deproto':
                    xyz_filename = sdf_path.replace('.sdf', '.xyz')
                    charge_filename = sdf_path.replace('.sdf', '.charge')
                else:
                    xyz_filename = sdf_path.replace('.sdf', '.xyz')
                    charge_filename = sdf_path.replace('.sdf', '.charge')

                # Call your optimization/PCM function
                A = dft_obj.opti_PCM(mol1, xyz_filename, charge_filename, xc)

                print(f"Finished {os.path.basename(sdf_path)} | basis: {basis} | xc: {xc} in {time.time() - start_time_sub:.2f} s")

    print(f"\n{'='*70}")
    print("All calculations completed!")


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 4:
        print("Usage: python dft.py DFT <cluster_reps_dir> <dielectric_constant>")
        print("   or: python dft.py DFT_DEPROTO <cluster_reps_dir> <dielectric_constant>")
        sys.exit(1)
    mode = sys.argv[1]
    cluster_reps_dir = sys.argv[2]
    dielectric_value = float(sys.argv[3])
    if mode == "DFT":
        run_DFT(cluster_reps_dir, dielectric_value, mode='normal')
    elif mode == "DFT_DEPROTO":
        run_DFT(cluster_reps_dir, dielectric_value, mode='deproto')
    else:
        print("Unknown mode, must be DFT or DFT_DEPROTO")
        sys.exit(1)