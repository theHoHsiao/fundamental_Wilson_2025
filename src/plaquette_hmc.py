#!/usr/bin/env python3

from argparse import ArgumentParser

from flow_analysis.stats.autocorrelation import exp_autocorrelation_fit

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

from bootstrap import basic_bootstrap, get_rng
from plaquette import get_index_separation
from plots import save_or_show


def get_skip(plaquettes):
    initial_plaquette = plaquettes[0][1]
    for index, (traj_index, plaquette) in enumerate(plaquettes):
        if plaquette != initial_plaquette:
            return index
    return -1


def get_args():
    parser = ArgumentParser(
        description=(
            "Get ensemble plaquettes from HMC logs, "
            "compute mean and autocorrelation time, and plot"
        )
    )

    parser.add_argument(
        "hmc_filenames",
        nargs="+",
        metavar="hmc_filename",
        help="Filename of HMC log file",
    )
    parser.add_argument(
        "--plot_filename",
        default=None,
        help="Where to place the generated plot. Default is to display on screen.",
    )
    parser.add_argument(
        "--plot_styles",
        default="styles/paperdraft.mplstyle",
        help="Stylesheet to use for plots",
    )
    return parser.parse_args()


def get_plaquette(filename):
    result = {"nAS": 0, "nF": 0, "nADJ": 0, "nS": 0}
    plaquettes = []
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("[SYSTEM][0]MACROS="):
                if "-DREPR_ANTISYMMETRIC" in line:
                    rep = "AS"
                if "-DREPR_SYMMETRIC" in line:
                    rep = "S"
                if "-DREPR_FUNDAMENTAL" in line:
                    rep = "F"
                if "-DREPR_ADJOINT" in line:
                    rep = "ADJ"
                result["rep"] = rep

            if line.startswith("[FLOW][0]Starting a new run from a"):
                if "start" in result:
                    raise ValueError(f"Multiple starts in file {filename}")
                result["start"] = line.split()[6]
            if line.startswith("[GEOMETRY_INIT][0]Global size is"):
                Nt, Nx, Ny, Nz = map(int, line.split()[-1].split("x"))
                result["Nx"] = Nx
                result["Ny"] = Ny
                result["Nz"] = Nz
                result["Nt"] = Nt
            if line.startswith("[ACTION][10]Monomial"):
                if "type = gauge," in line:
                    result["beta"] = float(line.split()[-1])
                elif "type = rhmc" in line or "type = hmc" in line:
                    if "type = hmc," in line:
                        result[f"n{rep}"] += 2
                    if "type = rhmc," in line:
                        result[f"n{rep}"] += 1

                    mass = float(line.split()[10].strip(","))
                    if result.get(f"m{rep}", mass) != mass:
                        raise NotImplementedError(
                            "Non-degenerate masses not currently supported"
                        )
                    result[f"m{rep}"] = mass
                else:
                    raise NotImplementedError(
                        f"Monomial not recognised in file {filename}"
                    )
            if line.startswith("[MAIN][0]Trajectory"):
                trajectory = int(line.split()[1].strip("#.:"))
            if line.startswith("[MAIN][0]Initial plaquette"):
                plaquettes.append((0, float(line.split()[-1])))
            if line.startswith("[MAIN][0]Plaquette:"):
                plaquettes.append((trajectory, float(line.split()[-1])))

    separation = get_index_separation([index for index, plaquette in plaquettes])
    skip = get_skip(plaquettes)
    result["plaq_tau_exp"] = (
        exp_autocorrelation_fit([plaquette for index, plaquette in plaquettes[skip:]])
        * separation
    )
    result["plaq_therm"] = int(10 * result["plaq_tau_exp"].nominal_value)
    plaquettes_subset = [
        plaquette
        for index, plaquette in plaquettes
        if index >= skip + result["plaq_therm"]
        and index % int(result["plaq_tau_exp"].nominal_value) == 0
    ]
    result["avg_plaquette"] = basic_bootstrap(plaquettes_subset, get_rng(filename))
    return result


def plot(data):
    betas = sorted(set([datum["beta"] for datum in data]))
    colours = {"unit": "blue", "random": "red"}
    markers = {"unit": "x", "random": "o"}
    fig, axes = plt.subplots(
        nrows=len(betas),
        layout="constrained",
        figsize=(3.5, 1 + 1.5 * len(betas)),
    )

    for beta, ax in zip(betas, axes):
        subset = [datum for datum in data if datum["beta"] == beta]
        reps = set([datum["rep"] for datum in subset])
        rep = reps.pop()
        if len(reps) != 0:
            raise NotImplementedError("Only one rep per beta currently allowed.")

        ax.text(
            0.95, 0.95, f"$\\beta={beta}$", ha="right", va="top", transform=ax.transAxes
        )
        ax.set_ylabel(r"$\langle \mathcal {P} \rangle$")
        ax.set_xlabel(f"$am_0^{{{rep.lower()}}}$")
        for datum in subset:
            ax.errorbar(
                [datum[f"m{rep}"]],
                [datum["avg_plaquette"].nominal_value],
                yerr=[datum["avg_plaquette"].std_dev],
                color=colours[datum["start"]],
                marker=markers[datum["start"]],
                ls="none",
            )
        ax.xaxis.set_major_locator(ticker.MultipleLocator(0.01))

    no_plot = [np.nan]
    ax.errorbar(
        no_plot,
        no_plot,
        yerr=no_plot,
        ls="none",
        color=colours["unit"],
        marker=markers["unit"],
        label="Cold start",
    )
    ax.errorbar(
        no_plot,
        no_plot,
        yerr=no_plot,
        ls="none",
        color=colours["random"],
        marker=markers["random"],
        label="Hot start",
    )
    fig.legend(loc="outside upper center", ncols=2)
    return fig


def main():
    args = get_args()
    plt.style.use(args.plot_styles)
    data = [get_plaquette(filename) for filename in args.hmc_filenames]
    save_or_show(plot(data), args.plot_filename)


if __name__ == "__main__":
    main()
