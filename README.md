# %title

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13819431.svg)](https://doi.org/10.5281/zenodo.13819431)

The workflow in this repository performs
the analyses presented in the paper
[%title][paper].

## Requirements

- Conda, for example, installed from [Miniforge][miniforge]
- [Snakemake][snakemake], which may be installed using Conda
- LaTeX, for example, from [TeX Live][texlive]

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
    2. Download the `data.h5` file from [the data release][datarelease],
       and place it into the `data_assets` directory.

4. Download the `metadata.csv` file from [the data release][datarelease],
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

[datarelease]: TODO
[miniforge]: https://github.com/conda-forge/miniforge
[paper]: TODO
[snakemake]: https://snakemake.github.io
[snakemake-conda]: https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html
[texlive]: https://tug.org/texlive/
