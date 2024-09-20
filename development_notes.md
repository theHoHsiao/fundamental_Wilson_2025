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

## Adding code

- Add one or more Python files to the `src` directory
- Add rules to relevant modules in the `workflow/rules` directory
  - I think we only need 3-4 modules for this work:
    1. `gradient_flow.smk`
    2. `mass_spectrum.smk`
    3. (Possibly `pcac_mass.smk`?)
    4. `output.smk`, to contain all the data presentation layer
    5. Maybe one more module for intermediary fits that lead to the plots?
       That could also be in `output.smk` however.
- Add new Conda environments into `workflow/envs`
  if the code has different requirements to existing rules.
  (I'd anticipate we'll minimally need a `spectrum.yml` and a `tabulate.yml`,
  in addition to the existing `flows.yml`.)

To test individual steps,
you can pass an output file to `snakemake`.
For example,

    snakemake --cores 1 --use-conda intermediary_data/Sp4b6.7nAS3mAS-1.055T48L24/w0_mean.csv

will generate just the $w_0$ computation for the single ensemble specified.

Once we start producing output assets,
these want to be specified as `input:` arguments
to a `rule all:` in the main `workflow/Snakefile`.
