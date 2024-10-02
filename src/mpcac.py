#!/usr/bin/env python3

from argparse import ArgumentParser, FileType
import re

import h5py
import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat

from bootstrap import get_rng, sample_bootstrap_1d
from dump import dump_dict


def get_args():
    parser = ArgumentParser(
        description="Compute the PCAC mass from correlators in an HDF5 file"
    )

    parser.add_argument("h5file", help="The file to read")
    # parser.add_argument("--metadata", default=None, help="CSV of ensemble metadata")
    parser.add_argument(
        "--beta",
        type=float,
        default=None,
        help="The beta value of the ensemble to analyse",
    )
    parser.add_argument(
        "--mAS",
        type=float,
        default=None,
        help="The antisymmetric fermion mass of the ensemble to analyse",
    )
    parser.add_argument(
        "--Nt",
        type=int,
        default=None,
        help="The temporal extent of the ensemble to analyse",
    )
    parser.add_argument(
        "--Ns",
        type=int,
        default=None,
        help="The spatial extent of the ensemble to analyse",
    )
    parser.add_argument(
        "--plateau_start",
        type=int,
        default=None,
        help="Time slice at which plateau starts",
    )
    parser.add_argument(
        "--plateau_end", type=int, default=None, help="Time slice at which plateau ends"
    )
    parser.add_argument(
        "--min_trajectory",
        type=int,
        default=None,
        help="Lowest trajectory index to consider",
    )
    parser.add_argument(
        "--max_trajectory",
        type=int,
        default=None,
        help="Highest trajectory index to consider",
    )
    parser.add_argument(
        "--trajectory_step",
        type=int,
        default=1,
        help="Interval of trajectories to consider",
    )
    parser.add_argument(
        "--output_file",
        type=FileType("w"),
        default="-",
        help="Where to output the result. (Defaults to stdout.)",
    )
    parser.add_argument(
        "--effmass_plot_file",
        default=None,
        help="Where to output the effective mass plot. (Skipped if not specified)",
    )
    parser.add_argument(
        "--plot_styles",
        default="styles/paperdraft.mplstyle",
        help="Stylesheet to use for plots",
    )
    return parser.parse_args()


def get_correlators(ensembles, beta=None, mAS=None, Nt=None, Ns=None):
    candidate_ensembles = []
    for ensemble in ensembles.values():
        if beta is not None and ensemble.get("beta", {(): None})[()] != beta:
            continue
        if mAS is not None and (
            len(masses := ensemble.get("quarkmasses", [])) != 1 or masses[0] != mAS
        ):
            continue
        if Nt is not None and ensemble.get("lattice", [None])[0] != Nt:
            continue
        if Ns is not None and tuple(ensemble.get("lattice", [None])[-3:]) != (
            Ns,
            Ns,
            Ns,
        ):
            continue
        candidate_ensembles.append(ensemble)
    if len(candidate_ensembles) > 1:
        raise ValueError("Did not uniquely identify one ensemble.")
    elif len(candidate_ensembles) == 0:
        raise ValueError("No ensembles found.")
    else:
        return candidate_ensembles[0]


def get_g5_eff_mass(sampled_correlator):
    eff_mass_samples = np.arccosh(
        (sampled_correlator[:-2] + sampled_correlator[2:])
        / (2 * sampled_correlator[1:-1])
    )
    return np.mean(eff_mass_samples, axis=1), np.std(eff_mass_samples, axis=1)


def get_eff_mass(ensemble, min_trajectory=None, max_trajectory=None, trajectory_step=1):
    indices = np.asarray(
        [
            int(re.match(".*n([0-9]+)$", filename.decode()).groups()[0])
            for filename in ensemble["configurations"]
        ]
    )
    filtered_indices = (
        ((indices >= min_trajectory) if min_trajectory is not None else True)
        & ((indices <= max_trajectory) if max_trajectory is not None else True)
        & ((indices - (min_trajectory or 0)) % trajectory_step == 0)
    )

    g5 = ensemble["TRIPLET"]["g5"][:, filtered_indices]
    g5_g0g5_re = ensemble["TRIPLET"]["g5_g0g5_re"][:, filtered_indices]

    g5_samples = sample_bootstrap_1d(g5.T, get_rng(ensemble.name)).T
    g5_g0g5_re_samples = sample_bootstrap_1d(g5_g0g5_re.T, get_rng(ensemble.name)).T

    g5_eff_mass = get_g5_eff_mass(g5_samples)[0][:, np.newaxis]

    # The factor g5_eff_mass / sinh(g5_eff_mass) is a correction to
    # make the derivative more closely match the continuum value.
    # See eq. B.14 of 1104.4301.
    eff_mass_samples = (
        g5_eff_mass
        / np.sinh(g5_eff_mass)
        * ((g5_g0g5_re_samples[:-2] - g5_g0g5_re_samples[2:]) / (4 * g5_samples[1:-1]))
    )

    eff_mass_mean = np.mean(eff_mass_samples, axis=1)
    eff_mass_error = np.std(eff_mass_samples, axis=1)

    return [ufloat(mean, error) for mean, error in zip(eff_mass_mean, eff_mass_error)]


def weighted_mean_by_uncertainty(results):
    numerator = sum([result / result.std_dev**2 for result in results])
    denominator = sum([1 / result.std_dev**2 for result in results])
    mean = numerator / denominator
    return mean


def fit_eff_mass(eff_mass, plaq_start, plaq_end):
    return weighted_mean_by_uncertainty(eff_mass[plaq_start:plaq_end])


def plot_eff_mass(eff_mass, fitted_mass, plot_filename):
    fig, ax = plt.subplots(layout="constrained")

    ax.errorbar(
        np.arange(len(eff_mass)),
        [mass.nominal_value for mass in eff_mass],
        [mass.std_dev for mass in eff_mass],
        ls="none",
        marker=".",
        label="Effective mass",
        color="C0",
    )
    ax.plot(
        [0, len(eff_mass) - 1],
        [fitted_mass.nominal_value] * 2,
        color="C1",
        label="Fitted mass",
    )
    ax.fill_between(
        [0, len(eff_mass) - 1],
        [fitted_mass.nominal_value - fitted_mass.std_dev] * 2,
        [fitted_mass.nominal_value + fitted_mass.std_dev] * 2,
        color="C1",
        alpha=0.4,
    )
    ax.set_xlabel("$t / a$")
    ax.set_ylabel("$m_{\\mathrm{PCAC}}$")

    ax.legend(loc="best")
    fig.savefig(plot_filename)


def main():
    args = get_args()
    plt.style.use(args.plot_styles)

    data = h5py.File(args.h5file, "r")
    ensemble = get_correlators(
        data, beta=args.beta, mAS=args.mAS, Nt=args.Nt, Ns=args.Ns
    )
    eff_mass = get_eff_mass(
        ensemble, args.min_trajectory, args.max_trajectory, args.trajectory_step
    )
    fitted_mass = fit_eff_mass(eff_mass, args.plateau_start, args.plateau_end)

    if args.effmass_plot_file:
        plot_eff_mass(eff_mass, fitted_mass, args.effmass_plot_file)

    dump_dict(
        {
            "beta": args.beta,
            "mAS": args.mAS,
            "Nt": args.Nt,
            "Ns": args.Ns,
            "mPCAC": fitted_mass,
        },
        args.output_file,
    )


if __name__ == "__main__":
    main()
