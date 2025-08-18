import argparse
import os
import pandas as pd
from utils.molconsul_helpers import run_molconsul_workflow

def main():
    parser = argparse.ArgumentParser(description="pKa Estimation Pipeline CLI.")
    subparsers = parser.add_subparsers(dest="command", help="Available commands", required=True)

    # --- MolConSUL Command ---
    parser_mc = subparsers.add_parser("molconsul", help="Run the MolConSUL (RDKit) workflow.")
    
    input_group = parser_mc.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--smiles", type=str, help="Single SMILES string to process.")
    input_group.add_argument("--csv", type=str, help="Path to a CSV file with 'smiles' and 'id' columns.")

    parser_mc.add_argument("--output_dir", type=str, default=".", help="Directory to save output files. Defaults to the current directory.")
    parser_mc.add_argument("--dielectric", type=float, default=46.826, help="Dielectric constant of the solvent.")
    parser_mc.add_argument("--num_confs", type=int, default=500, help="Number of conformers to generate.")
    parser_mc.add_argument("--e_avg_proton", type=float, default=-277.60, help="Average proton energy E_H(solv) in kcal/mol.")
    parser_mc.add_argument("--pka_exp", type=float, default=0.0, help="Experimental pKa (optional).")

    args = parser.parse_args()

    # Define the base path within the container
    container_base_path = "/work"

    if args.command == "molconsul":
        main_output_dir = os.path.join(container_base_path, args.output_dir)
        os.makedirs(main_output_dir, exist_ok=True)

        if args.smiles:
            print(f"--- Running MolConSUL for SMILES: {args.smiles} ---")
            molecule_output_dir = os.path.join(main_output_dir, args.smiles.replace('/', '_').replace('\\', '_'))
            run_molconsul_workflow(
                smiles=args.smiles,
                output_dir=molecule_output_dir,
                args=args
            )
        elif args.csv:
            print(f"--- Running MolConSUL batch from: {args.csv} ---")
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
                    
                    run_molconsul_workflow(
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