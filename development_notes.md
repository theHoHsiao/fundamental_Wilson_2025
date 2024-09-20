# Contributing to this repository

To get started:

1. Clone the repository as described in `README.md`
2. If you don't already have one,
   create a Conda environment with Snakemake and pre-commit installed:

       conda create -c conda-forge -n snakemake snakemake
       conda activate snakemake
       pip install pre-commit

3. Install `pre-commit` into the repository,
   so that basic code quality checks will be performed when you commit.

       pre-commit install

## Directory structure

These things are committed:

- Helper tools not used by the workflow go in `tools`
- Tools that form part of the main workflow go in `src`
- The workflow definition and associated files go in `workflow`

These things are *not* committed:

- Raw data go in `raw_data`, with a subdirectory by class of data.
  E.g. `raw_data/flows/Sp4b6.8nAS3mAS-1.45L28T32/out_wflow`
- Metadata go in `metadata`.
- Things to include in the paper go in a relevant subdirecotry in `assets`.
- Things to include in the data release go in `data_assets`.
- Intermediary things not to go anywhere go in `intermediary_data`
