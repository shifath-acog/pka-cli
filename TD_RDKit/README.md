Launch a container using the Docker image:
dockerhub.aganitha.ai:4443/chem/pka_pipeline:ovh3

Inside the container, make sure the following are present in the /work directory:
```
/work
│
├── run.ipynb                          # Main notebook to be executed (**update parameters here**)
├── pka_pipeline_TD.ipynb             # Notebook for generating conformers using Torsional Diffusion
├── full_pipeline_pka_rdkit_TD.ipynb  # Second stage of the pipeline (example name)
└── csv                               # Directory containing the required CSV file with SMILES
```

Change the parameters as you required run.ipynb
```
    inp_smiles = row['smiles']
    num_conf_rdkit: int = 200           # Number of Conformers you want to generate using RDKit
    num_conf_TD: int = 200              # Number of Conformers you want to generate using TD
    E_avg_proton: float=-275.755
    pKa_EXP: float= 24.88
    Solvent: str="DMSO"
    dielectric_value = 46.826
```
