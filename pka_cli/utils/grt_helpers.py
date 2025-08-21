import os
import subprocess
import pandas as pd
from rdkit import Chem
import shutil

from .common_pka_workflow import run_common_pka_workflow
from atk_conformer_generation_pipeline.utils import (
    generate_conformers,
    mmff_optimize_conformers,
    save_conformers_to_sdf,
    process_conformers,
    calculate_relative_energies,
    convert_conformers_to_smiles
)
from atk_conformer_generation_pipeline.variables import *

def run_td_conformer_generation(smiles, num_confs, output_dir):
    """
    Generates conformers using the Torsional Diffusion method.
    """
    print("--- Running Torsional Diffusion Conformer Generation ---", flush=True)
    td_dir = os.path.join(output_dir, "td_temp")
    os.makedirs(td_dir, exist_ok=True)
    original_cwd = os.getcwd()
    os.chdir(td_dir)

    try:
        with open("isomeric_smi_inp.csv", 'w') as f:
            f.write('smile_str,num_conformers,smile_str\n')
            f.write(f'{smiles},0,{smiles}\n')

        td_command = [
            'conda', 'run', '-n', 'torsional_diffusion', 'python', '/app/torsional-diffusion-master/generate_confs.py',
            '--model_dir', '/app/TD-trained_models/drugs_default/',
            '--out', 'conformers_TD.pkl',
            '--test_csv', 'isomeric_smi_inp.csv',
            '--inference_steps', '20',
            '--confs_per_mol', str(num_confs),
            '--tqdm', '--batch_size', '64', '--post_mmff'
        ]
        print(f"Executing TD...", flush=True)
        subprocess.run(td_command, check=True, capture_output=True, text=True)

        pkl_command = [
            'conda', 'run', '-n', 'torsional_diffusion', 'python', '/app/pkl_to_sdf.py',
            'conformers_TD.pkl',
            'conformers_TD.sdf'
        ]
        print(f"Converting PKL to SDF...", flush=True)
        subprocess.run(pkl_command, check=True, capture_output=True, text=True)

        shutil.move("conformers_TD.sdf_0_conformers.sdf", "conformers_TD.sdf")
        
        with open("conformers_TD.sdf", 'r') as f:
            lines = f.readlines()
        with open("conformers_TD.sdf", 'w') as f:
            if lines:
                f.writelines(lines[1:])

        shutil.copy("conformers_TD.sdf", os.path.join(output_dir, "conformers_TD.sdf"))
        print("--- Torsional Diffusion complete ---", flush=True)

    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Torsional Diffusion failed.", flush=True)
        print(f"STDOUT:\n{e.stdout}", flush=True)
        print(f"STDERR:\n{e.stderr}", flush=True)
        raise
    finally:
        os.chdir(original_cwd)
        shutil.rmtree(td_dir)

def run_gnnis_via_subprocess(smiles, num_confs, solvent, output_dir, python_executable):
    """
    Runs the GNNIS workflow as a subprocess by calling the main CLI script.
    """
    print("--- Calling GNNIS subcommand for conformer generation ---", flush=True)
    gnnis_dir = os.path.join(output_dir, "gnnis_temp")
    
    # The gnnis command creates its own subdirectory based on the SMILES string.
    # We need to account for that to find the output file.
    safe_smiles_name = smiles.replace('/', '_').replace('\\', '_')
    expected_output_file = os.path.join(gnnis_dir, safe_smiles_name, "optimized_generated_conformers.sdf")

    gnnis_command = [
        python_executable,
        "/work/pka_cli/main.py", "gnnis",
        "--smiles", smiles,
        "--output_dir", gnnis_dir,
        "--num_confs", str(num_confs),
        "--solvent_name", solvent
    ]
    
    try:
        subprocess.run(gnnis_command, check=True, capture_output=True, text=True)
        # Copy the final SDF back to the main output directory
        shutil.copy(expected_output_file, os.path.join(output_dir, "conformers_gnnis.sdf"))
        print("--- GNNIS complete ---", flush=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] GNNIS subprocess failed.", flush=True)
        print(f"STDOUT:\n{e.stdout}", flush=True)
        print(f"STDERR:\n{e.stderr}", flush=True)
        raise
    finally:
        if os.path.exists(gnnis_dir):
            shutil.rmtree(gnnis_dir)

def combine_sdfs(sdf_files, output_sdf):
    """
    Combines multiple SDF files into a single file.
    """
    print(f"Combining {len(sdf_files)} SDF files into {output_sdf}...", flush=True)
    with open(output_sdf, 'wb') as wfd:
        for f in sdf_files:
            if os.path.exists(f):
                with open(f, 'rb') as fd:
                    shutil.copyfileobj(fd, wfd)
            else:
                print(f"[Warning] SDF file not found, skipping: {f}", flush=True)

def run_grt_workflow(smiles, output_dir, args, base_python_executable):
    """
    Main workflow for the GRT pKa estimation pipeline.
    """
    print(f"Running GRT workflow for {smiles}")
    print(f"Output will be saved in: {output_dir}")
    
    os.makedirs(output_dir, exist_ok=True)
    original_cwd = os.getcwd()
    os.chdir(output_dir)

    try:
        # --- MODULE 1: Multi-method Conformer Generation ---
        run_td_conformer_generation(smiles, args.num_confs_td, output_dir)
        run_gnnis_via_subprocess(smiles, args.num_confs_gnnis, args.solvent_name, output_dir, base_python_executable)

        print("--- Running RDKit Conformer Generation ---", flush=True)
        rdkit_mol = generate_conformers(smiles, args.num_confs_rdkit)
        rdkit_opt_mol, _ = mmff_optimize_conformers(rdkit_mol)
        save_conformers_to_sdf(rdkit_opt_mol, "conformers_RDKit.sdf")
        print("--- RDKit Generation complete ---", flush=True)

        sdf_to_combine = ["conformers_TD.sdf", "conformers_gnnis.sdf", "conformers_RDKit.sdf"]
        combine_sdfs(sdf_to_combine, init_conf_sdf)

        print("--- Filtering combined conformers ---", flush=True)
        filter_script_path = "/app/dihedral_filter.py"
        filtered_sdf = "filtered_initial_generated_conformers.sdf"
        filter_command = [base_python_executable, filter_script_path, smiles, init_conf_sdf, filtered_sdf]
        subprocess.run(filter_command, check=True, capture_output=True, text=True)

        print("--- Performing final MMFF optimization on combined set ---", flush=True)
        suppl = Chem.SDMolSupplier(filtered_sdf, removeHs=False)
        mols = [m for m in suppl if m is not None]
        if not mols:
            raise ValueError("No valid molecules left after filtering.")
        
        # Re-combine into a single Mol object
        combined_mol = Chem.Mol(mols[0])
        combined_mol.RemoveAllConformers()
        for m in mols:
            for conf in m.GetConformers():
                combined_mol.AddConformer(conf, assignId=True)

        opt_mol, conformer_energies = mmff_optimize_conformers(combined_mol)
        save_conformers_to_sdf(opt_mol, opt_conf_sdf)
        energy_df = pd.DataFrame(list(conformer_energies.items()), columns=['conformer_id', 'energy_in_kcalpermol'])
        
        convert_conformers_to_smiles(opt_conf_sdf, opt_conf_SMILES_file)
        _, energy_df = process_conformers(
            opt_conf_SMILES_file, opt_conf_sdf, feasible_geometries_sdf,
            infeasible_geometries_sdf, similarity_output_csv, infeasible_geometries_csv,
            smiles, opt_mol.GetNumConformers(), energy_df
        )
        
        # --- MODULE 2 & 3 ---
        print("\nStep 2 & 3: Clustering, DFT & pKa Calculation...", flush=True)
        run_common_pka_workflow(
            opt_mol=opt_mol,
            energy_df=energy_df,
            smiles=smiles,
            output_dir=output_dir,
            dielectric_value=args.dielectric,
            e_avg_proton=args.e_avg_proton,
            pka_exp=args.pka_exp,
            opt_conf_sdf=opt_conf_sdf,
            base_python_executable=base_python_executable
        )

    except Exception as e:
        print(f"[ERROR] GRT workflow failed: {e}", flush=True)
        raise
    finally:
        os.chdir(original_cwd)
        print(f"--- GRT Workflow for {smiles} complete ---")