

# Titratable Sites Protonation/Deprotonation Tool

Usage:
python main.py name_of_folder_with_sdf_and_charge_files/

Please make sure that .sdf and .charge files are in the same folder with same name (extension has to be  .charge and .sdf for charge and sdf file, respectively).

Outputs are save in a new folder named "treated" in the input_folder. A specific extention "treated" is also appended to the original name, in order to easily understannd the output file.

## Overview

This tool processes molecules (SDF files) to identify and modify titratable functional groups through protonation and deprotonation steps. It is designed for cheminformatics workflows involving molecular ionization state preparation, such as pKa estimation.

---

## Key Features
- **Titratable site identification**:
  The tool first looks for possible titratable sites by identifying groups like -OH, -SH or amine functionalities etc.

- **Treatment of titratable sites**:
  Once the titratable sites has been identified, this tool deprotonates the acidic functionalities, or protonates the basic functionalities.

- **Ambiphilic Molecule Handling**:  
  Molecules containing both acidic and basic groups are treated in both directions, producing two outputs.

- **Force Field Optimization**:  
  After protonation, the molecule is optimized using RDKit’s UFF force field to refine geometry. This does not make much change in the overall strcutre but it is necessary since too much irregularities in geometry of amine functionality after protonation might leads to the error in DFT calculations.

- **ranking**:  
  Ranking is done based on DFT computed Mulliken charges. Among acidic sites, most positive centre is deprotonated vice-versa for basic sites.
    Protonated molecules are simply deprotonated because that seems to be the most logical operation to calculate pKa
    Anionic compounds are simply protonated at the site of -ve charge
    Zwitter-ionic compounds are deprotonated from the local site where valancy is exceeding provided that site has atleast one hydrogen attached.

ionization_tool_with_ranking
  ├── main.py
├── protodeproto
│   ├── __init__.py
│   ├── __pycache__
│   ├── acidic.py
│   ├── ambiphilic.py
│   ├── basic.py
│   ├── fix_anionic.py
│   └── protonated.py
├── README.md
└── sdf_files
    ├── molecule1
    ├── molecule2
    └── molecule3

