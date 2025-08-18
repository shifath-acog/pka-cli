# pKa Estimation CLI Project Plan

This plan outlines the steps to create a unified command-line interface (CLI) for the pKa estimation pipelines.

## Phase 1: Project Scaffolding

- [ ] Create a new directory for the CLI project (e.g., `pka_cli`).
- [ ] Create the main entry point script (e.g., `pka_cli.py`).
- [ ] Use Python's `argparse` library to set up the main command and sub-commands for each method (`molconsul`, `gnnis`, `md-conformers`, `td-gnnis-rdkit`).
- [ ] Create a `src` directory to hold the refactored Python code from the notebooks.

## Phase 2: Implement the `molconsul` Command

- [ ] Create a module `src/molconsul_pipeline.py`.
- [ ] Refactor the code from `full_pipeline_pka_v0.3.ipynb` into functions within this module.
- [ ] Connect the `molconsul` sub-command in `pka_cli.py` to call these functions.
- [ ] Test the `molconsul` command from the terminal.

## Phase 3: Implement the `gnnis` Command

- [ ] Create a module `src/gnnis_pipeline.py`.
- [ ] Analyze and refactor the `gnnis` notebooks into functions.
- [ ] Connect the `gnnis` sub-command to the new functions.
- [ ] Test the `gnnis` command.

## Phase 4: Implement the `td-gnnis-rdkit` Command

- [ ] Create a module `src/td_gnnis_rdkit_pipeline.py`.
- [ ] Analyze the chain of notebooks for this method to understand the workflow.
- [ ] Refactor the logic into functions.
- [ ] Connect the `td-gnnis-rdkit` sub-command.
- [ ] Test the command.

## Phase 5: Stub out `md-conformers`

- [ ] Since this is in development, add the `md-conformers` sub-command but have it print a "Coming Soon" message.

## Phase 6: Final Touches

- [ ] Add help messages and documentation for all commands.
- [ ] Implement robust error handling.
- [ ] Create a `README.md` for the new CLI tool.
