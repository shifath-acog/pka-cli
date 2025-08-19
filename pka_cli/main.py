# import argparse
# import os
# import pandas as pd
# from utils.molconsul_helpers import run_molconsul_workflow

# def main():
#     parser = argparse.ArgumentParser(description="pKa Estimation Pipeline CLI.")
#     subparsers = parser.add_subparsers(dest="command", help="Available commands", required=True)

#     # --- MolConSUL Command ---
#     parser_mc = subparsers.add_parser("molconsul", help="Run the MolConSUL (RDKit) workflow.")
    
#     input_group = parser_mc.add_mutually_exclusive_group(required=True)
#     input_group.add_argument("--smiles", type=str, help="Single SMILES string to process.")
#     input_group.add_argument("--csv", type=str, help="Path to a CSV file with 'smiles' and 'id' columns.")

#     parser_mc.add_argument("--output_dir", type=str, default=".", help="Directory to save output files. Defaults to the current directory.")
#     parser_mc.add_argument("--dielectric", type=float, default=46.826, help="Dielectric constant of the solvent.")
#     parser_mc.add_argument("--num_confs", type=int, default=500, help="Number of conformers to generate.")
#     parser_mc.add_argument("--e_avg_proton", type=float, default=-277.60, help="Average proton energy E_H(solv) in kcal/mol.")
#     parser_mc.add_argument("--pka_exp", type=float, default=0.0, help="Experimental pKa (optional).")

#     args = parser.parse_args()

#     # Define the base path within the container
#     container_base_path = "/work"

#     if args.command == "molconsul":
#         main_output_dir = os.path.join(container_base_path, args.output_dir)
#         os.makedirs(main_output_dir, exist_ok=True)

#         if args.smiles:
#             print(f"--- Running MolConSUL for SMILES: {args.smiles} ---")
#             molecule_output_dir = os.path.join(main_output_dir, args.smiles.replace('/', '_').replace('\\', '_'))
#             run_molconsul_workflow(
#                 smiles=args.smiles,
#                 output_dir=molecule_output_dir,
#                 args=args
#             )
#         elif args.csv:
#             print(f"--- Running MolConSUL batch from: {args.csv} ---")
#             csv_path = os.path.join(container_base_path, args.csv)
#             try:
#                 df = pd.read_csv(csv_path)
#                 if not all(col in df.columns for col in ['smiles', 'id']):
#                     print("[ERROR] CSV file must contain 'smiles' and 'id' columns.")
#                     return

#                 for index, row in df.iterrows():
#                     molecule_id = row['id']
#                     smiles = row['smiles']
#                     print(f"\n--- Processing molecule {molecule_id}: {smiles} ---")
                    
#                     molecule_output_dir = os.path.join(main_output_dir, str(molecule_id))
                    
#                     run_molconsul_workflow(
#                         smiles=smiles,
#                         output_dir=molecule_output_dir,
#                         args=args
#                     )

#             except FileNotFoundError:
#                 print(f"[ERROR] CSV file not found at: {csv_path}")
#             except Exception as e:
#                 print(f"[ERROR] An error occurred while processing the CSV file: {e}")

# if __name__ == "__main__":
#     main()


import argparse
import os
import pandas as pd
import sys
import subprocess

# --- Helper Functions ---
def is_in_gnnis_env():
    """Check if the script is running inside the GNNIS conda environment."""
    # Check if the current executable path contains the GNNImplicitSolvent environment name
    return 'GNNImplicitSolvent' in sys.executable

def relaunch_in_gnnis_env(base_python_executable):
    """Re-launches the script within the GNNIS conda environment."""
    print("--- Not in GNNIS environment. Re-launching... ---")
    try:
        # Construct the command to re-run the script with conda run
        # sys.argv[0] is the script name, sys.argv[1:] are the arguments
        conda_command = ['conda', 'run', '-n', 'GNNImplicitSolvent', 'python', sys.argv[0]] + sys.argv[1:] + ['--base_python_executable', base_python_executable]
        print(f"Executing: {' '.join(conda_command)}")
        # Execute the command and exit, passing stdout/stderr through
        subprocess.run(conda_command, check=True)
        sys.exit(0) # Exit after the re-launched process finishes
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Failed to re-launch in GNNIS environment: {e}")
        sys.exit(1)
    except FileNotFoundError:
        print("[ERROR] 'conda' command not found. Please ensure conda is installed and in your PATH.")
        sys.exit(1)

# --- Main Application Logic ---
def main():
    parser = argparse.ArgumentParser(description="pKa Estimation Pipeline CLI.")
    subparsers = parser.add_subparsers(dest="command", help="Available commands", required=True)

    # --- MolConSUL Command ---
    parser_mc = subparsers.add_parser("molconsul", help="Run the MolConSUL (RDKit) workflow.")
    mc_input_group = parser_mc.add_mutually_exclusive_group(required=True)
    mc_input_group.add_argument("--smiles", type=str, help="Single SMILES string to process.")
    mc_input_group.add_argument("--csv", type=str, help="Path to a CSV file with 'smiles' and 'id' columns.")
    parser_mc.add_argument("--output_dir", type=str, default=".", help="Directory to save output files. Defaults to the current directory.")
    parser_mc.add_argument("--dielectric", type=float, default=46.826, help="Dielectric constant of the solvent.")
    parser_mc.add_argument("--num_confs", type=int, default=500, help="Number of conformers to generate.")
    parser_mc.add_argument("--e_avg_proton", type=float, default=-277.60, help="Average proton energy E_H(solv) in kcal/mol.")
    parser_mc.add_argument("--pka_exp", type=float, default=0.0, help="Experimental pKa (optional).")

    # --- GNNIS Command ---
    parser_gn = subparsers.add_parser("gnnis", help="Run the GNNIS (GNNImplicitSolvent) workflow.")
    gn_input_group = parser_gn.add_mutually_exclusive_group(required=True)
    gn_input_group.add_argument("--smiles", type=str, help="Single SMILES string to process.")
    gn_input_group.add_argument("--csv", type=str, help="Path to a CSV file with 'smiles' and 'id' columns.")
    parser_gn.add_argument("--output_dir", type=str, default=".", help="Directory to save output files. Defaults to the current directory.")
    parser_gn.add_argument("--solvent_name", type=str, default="tip3p", help="Name of the solvent for GNNIS (e.g., 'tip3p', 'DMSO').")
    parser_gn.add_argument("--dielectric", type=float, default=46.826, help="Dielectric constant of the solvent for DFT.")
    parser_gn.add_argument("--num_confs", type=int, default=10, help="Number of conformers to generate.")
    parser_gn.add_argument("--e_avg_proton", type=float, default=-277.60, help="Average proton energy E_H(solv) in kcal/mol.")
    parser_gn.add_argument("--pka_exp", type=float, default=0.0, help="Experimental pKa (optional).")
    # Add base_python_executable to the gnnis sub-parser so it's recognized after re-launch
    parser_gn.add_argument("--base_python_executable", type=str, help=argparse.SUPPRESS) # Hidden argument

    args = parser.parse_args()

    # Determine the base_python_executable before any potential re-launch
    base_python_executable = args.base_python_executable or sys.executable

    # If the command is gnnis, ensure we are in the correct conda environment
    if args.command == 'gnnis' and not is_in_gnnis_env():
        relaunch_in_gnnis_env(base_python_executable)

    # Import workflow functions AFTER potential re-launch and environment check
    # This ensures that GNNIS-related imports only happen when in the correct env
    from utils.molconsul_helpers import run_molconsul_workflow
    # Only import gnnis_helpers if the command is 'gnnis'
    if args.command == 'gnnis':
        from utils.gnnis_helpers import run_gnnis_workflow

    container_base_path = "/work"

    if args.command == "molconsul" or args.command == "gnnis":
        main_output_dir = os.path.join(container_base_path, args.output_dir)
        os.makedirs(main_output_dir, exist_ok=True)

        workflow_func = None
        if args.command == "molconsul":
            workflow_func = run_molconsul_workflow
        elif args.command == "gnnis":
            workflow_func = run_gnnis_workflow

        if args.smiles:
            print(f"--- Running {args.command.upper()} for SMILES: {args.smiles} ---")
            molecule_output_dir = os.path.join(main_output_dir, args.smiles.replace('/', '_').replace('\\', '_'))
            if args.command == 'gnnis':
                workflow_func(
                    smiles=args.smiles,
                    output_dir=molecule_output_dir,
                    args=args,
                    base_python_executable=base_python_executable
                )
            else:
                workflow_func(
                    smiles=args.smiles,
                    output_dir=molecule_output_dir,
                    args=args
                )
        elif args.csv:
            print(f"--- Running {args.command.upper()} batch from: {args.csv} ---")
            csv_path = os.path.join(container_base_path, args.csv)
            try:
                df = pd.read_csv(csv_path)
                if not all(col in df.columns for col in ['smiles', 'id']):
                    print("[ERROR] CSV file must contain 'smiles' and 'id' columns.")
                    return

                for index, row in df.iterrows():
                    molecule_id = row['id']
                    smiles = row['smiles']
                    print(f"\n--- Processing molecule {molecule_id}: {smiles} ---")
                    
                    molecule_output_dir = os.path.join(main_output_dir, str(molecule_id))
                    
                    if args.command == 'gnnis':
                        workflow_func(
                            smiles=smiles,
                            output_dir=molecule_output_dir,
                            args=args,
                            base_python_executable=base_python_executable
                        )
                    else:
                        workflow_func(
                            smiles=smiles,
                            output_dir=molecule_output_dir,
                            args=args
                        )

            except FileNotFoundError:
                print(f"[ERROR] CSV file not found at: {csv_path}")
            except Exception as e:
                print(f"[ERROR] An error occurred while processing the CSV file: {e}")

if __name__ == "__main__":
    main()

