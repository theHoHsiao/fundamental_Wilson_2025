# Meson spectroscopy in the Sp(4) gauge theory with three antisymmetric fermions&mdash;Analysis workflow

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13819431.svg)](https://doi.org/10.5281/zenodo.13819431)

The workflow in this repository performs
the analyses presented in the paper
[Meson spectroscopy in the Sp(4) gauge theory with three antisymmetric fermions][paper].

## Requirements

- Conda, for example, installed from [Miniforge][miniforge]
- [Snakemake][snakemake], which may be installed using Conda
- LaTeX, for example, from [TeX Live][texlive]

### Notes for Apple silicon users

This project uses Julia,
installed via a Conda environment,
which at time of writing is not available on Apple silicon.
To automate the environment setup process for the workflow,
an x86-64 version of Snakemake must be used.

To set this up,
create a new x86-64 Conda environment with Snakemake,
using

``` shellsession
conda create -n snakemake_x86 -c conda-forge -c bioconda snakemake
```

then activate this and install Mamba.
At time of writing,
Snakemake is not compatible with versions of Mamba from 2.0.0 onwards,
so this version must be constrained for the Conda integration to work correctly:

``` shellsession
conda activate snakemake_x86
conda install -c conda-forge 'mamba<2.0.0'
```

With this environment active,
the steps below involving running `snakemake` should work correctly.

## Setup

1. Install the dependencies above.
2. Clone this repository including submodules
   (or download its Zenodo release and `unzip` it)
   and `cd` into it:

   ```shellsession
   git clone --recurse-submodules https://github.com/telos-collaboration/antisymmetric_analysis_2024
   cd antisymmetric_analysis_2024
   ```

3. Either:

    1. Download the `raw_data.zip` file from [the data release][datarelease],
       and extract it into the root of the repository,
       or
    2. Download the `correlators_smear.h5`,
       `correlators_wall.h5`,
       and `flows.h5`
       files from [the data release][datarelease],
       and place them into the `data_assets` directory.
       Instruct Snakemake that these files are up to date
       by running

       ```shellsession
       snakemake --touch data_assets/{correlators_smear,correlators_wall,flows}.h5
       ```

4. Download the `ensemble_metadata.csv` file from [the data release][datarelease],
   and place it into the `metadata` directory.

## Running the workflow

The workflow is run using Snakemake:

``` shellsession
snakemake --cores 1 --use-conda
```

where the number `1`
may be replaced by
the number of CPU cores you wish to allocate to the computation.

Snakemake will automatically download and install
all required Python packages.
This requires an Internet connection;
if you are running in an HPC environment where you would need
to run the workflow without Internet access,
details on how to preinstall the environment
can be found in the [Snakemake documentation][snakemake-conda].

Using `--cores 6` on a MacBook Pro with an M1 Pro processor,
the analysis takes around TODO minutes.

## Output

Output plots, tables, and definitions
are placed in the `assets/plots`, `assets/tables`, and `assets/definitions` directories.

Output data assets are placed into the `data_assets` directory.

Intermediary data are placed in the `intermediary_data` directory.

## Reusability

This workflow is relatively tailored to the data
which it was originally written to analyse.
Additional ensembles may be added to the analysis
by adding relevant files to the `raw_data` directory,
and adding corresponding entries to the files in the `metadata` directory.
Tools present in the `tools` directory may be used
to help determine some inputs for the latter;
for example,
the plateaux positions and lengths.
However,
extending the analysis in this way
has not been as fully tested as the rest of the workflow,
and is not guaranteed to be trivial for someone not already familiar with the code.

[datarelease]: https://doi.org/10.5281/zenodo.13819562
[miniforge]: https://github.com/conda-forge/miniforge
[paper]: TODO
[snakemake]: https://snakemake.github.io
[snakemake-conda]: https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html
[texlive]: https://tug.org/texlive/
