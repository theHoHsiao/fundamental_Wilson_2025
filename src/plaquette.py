#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

from flow_analysis.stats.autocorrelation import exp_autocorrelation_fit

import h5py
from more_itertools import pairwise
import pandas as pd
import re

from bootstrap import basic_bootstrap, get_rng


def get_args():
    parser = ArgumentParser(
        description=(
            "Get ensemble plaquettes from HDF5, "
            "compute mean and autocorrelation time, and write a table"
        )
    )

    parser.add_argument("h5file", help="The file to read")
    parser.add_argument("--metadata", default=None, help="CSV of ensemble metadata")
    parser.add_argument(
        "--output_table",
        type=FileType("w"),
        default="-",
        help="Where to output the LaTeX table. (Defaults to stdout.)",
    )
    return parser.parse_args()


def get_volumes(ensemble):
    nt, nx, ny, nz = ensemble["lattice"]
    if nx != ny or ny != nz:
        raise NotImplementedError("This code expects NX == NY == NZ")
    return {"NT": nt, "NS": nx}


def get_mass(ensemble):
    if len(ensemble["quarkmasses"]) != 1:
        raise NotImplementedError("This code expects one fermion mass per ensemble.")
    return ensemble["quarkmasses"][0]


def get_index_separation(indices):
    separation = indices[1] - indices[0]
    for idx1, idx2 in pairwise(indices):
        if idx2 - idx1 != separation:
            raise NotImplementedError(
                "Configurations have non-uniform separation or are out of order."
            )
    return separation


def get_cfg_separation(cfgs):
    indices = [
        int(re.match(".*n([0-9]+)$", filename.decode()).groups()[0])
        for filename in cfgs
    ]
    return get_index_separation(indices)


def get_name(ensemble, metadata):
    if metadata is None:
        return None
    subset = metadata[
        (metadata.Nt == ensemble["NT"])
        & (metadata.Ns == ensemble["NS"])
        & (metadata.beta == ensemble["beta"])
        & (metadata.mAS == ensemble["mAS"])
    ]
    if len(subset) == 0:
        return None
    elif len(subset) > 1:
        raise ValueError("Too many ensembles found")
    else:
        return subset.ensemble_name.values[0]


def process_ensemble(ensemble, metadata):
    result = {}
    result.update(get_volumes(ensemble))
    result["mAS"] = get_mass(ensemble)
    result["beta"] = ensemble["beta"][()]
    result["Ncfg"] = len(ensemble["configurations"])
    result["delta_traj"] = get_cfg_separation(ensemble["configurations"])
    result["name"] = get_name(result, metadata)
    result["avg_plaq"] = basic_bootstrap(ensemble["plaquette"], get_rng(ensemble.name))
    result["plaq_autocorr"] = (
        exp_autocorrelation_fit(ensemble["plaquette"]) * result["delta_traj"]
    )
    return result


def tabulate(results, skip_missing_names=True):
    header = (
        "\\begin{tabular}{ccccccccc}\n\\hline\\hline\n"
        + " & ".join(
            [
                "Ensemble",
                "$N_t \\times N_s^3$",
                r"$\beta$",
                "$am_0$",
                r"$N_{\mathrm{cfg}}$",
                r"$\delta_{\mathrm{traj}}$",
                r"$\langle \mathcal{P} \rangle$",
                r"$\tau_{\mathrm{exp}}^{\mathcal{P}}$",
                "Comment",
            ]
        )
        + " \\\\\n"
    )
    footer = "\n\\hline\\hline\n\\end{tabular}"
    content = []
    previous_beta = None
    for result in sorted(results, key=lambda r: (r["beta"], -r["mAS"])):
        if skip_missing_names and result["name"] is None:
            continue
        if result["beta"] != previous_beta:
            previous_beta = result["beta"]
            content.append(r"\hline")
        content.append(
            " & ".join(
                [
                    result["name"],
                    f"${result['NT']} \\times {result['NS']}^3$",
                    f"{result['beta']}",
                    f"{result['mAS']}",
                    f"{result['Ncfg']}",
                    f"{result['delta_traj']}",
                    f"{result['avg_plaq']:.02uSL}",
                    f"{result['plaq_autocorr']:.02uSL}",
                    "tbc",
                ]
            )
            + r" \\"
        )
    return header + "\n".join(content) + footer


def main():
    args = get_args()
    data = h5py.File(args.h5file, "r")
    if args.metadata is not None:
        metadata = pd.read_csv(args.metadata)
    else:
        metadata = None
    results = [process_ensemble(ensemble, metadata) for ensemble in data.values()]
    print(tabulate(results), file=args.output_table)


if __name__ == "__main__":
    main()
