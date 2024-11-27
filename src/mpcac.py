#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

import h5py
import matplotlib.pyplot as plt
import numpy as np

from .bootstrap import get_rng, sample_bootstrap_1d, bootstrap_finalize
from .dump import dump_dict, dump_samples
from .read_hdf5 import get_ensemble, filter_configurations


def get_args():
    parser = ArgumentParser(
        description="Compute the PCAC mass from correlators in an HDF5 file"
    )

    parser.add_argument("h5file", help="The file to read")
    parser.add_argument(
        "--ensemble_name",
        default=None,
        help="Name of the ensemble to analyse. Only used for tagging output.",
    )
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
        "--output_file_mean",
        type=FileType("w"),
        default="-",
        help="Where to output the mean and uncertainty of mPCAC. (Defaults to stdout.)",
    )
    parser.add_argument(
        "--output_file_samples",
        type=FileType("w"),
        default=None,
        help="Where to output the bootstrap samples for mPCAC",
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


def get_g5_eff_mass(sampled_correlator):
    eff_mass_samples = (
        (sampled_correlator[:-2] + sampled_correlator[2:])
        / (2 * sampled_correlator[1:-1])
    ).arccosh()
    return eff_mass_samples.mean, eff_mass_samples.std()


def get_eff_mass_samples(
    ensemble, min_trajectory=None, max_trajectory=None, trajectory_step=1
):
    filtered_indices = filter_configurations(
        ensemble, min_trajectory, max_trajectory, trajectory_step
    )

    g5 = ensemble["TRIPLET"]["g5"][:, filtered_indices]
    g5_g0g5_re = ensemble["TRIPLET"]["g5_g0g5_re"][:, filtered_indices]

    g5_samples = sample_bootstrap_1d(g5.T, get_rng(ensemble.name))
    g5_g0g5_re_samples = sample_bootstrap_1d(g5_g0g5_re.T, get_rng(ensemble.name))

    g5_eff_mass = get_g5_eff_mass(g5_samples)[0]

    # The factor g5_eff_mass / sinh(g5_eff_mass) is a correction to
    # make the derivative more closely match the continuum value.
    # See eq. B.14 of 1104.4301.
    return (
        ((g5_g0g5_re_samples[:-2] - g5_g0g5_re_samples[2:]) / (4 * g5_samples[1:-1]))
        * g5_eff_mass
        / np.sinh(g5_eff_mass)
    )


def fit_eff_mass_samples(eff_mass_samples, plateau_start, plateau_end):
    return eff_mass_samples[plateau_start:plateau_end]


def plot_eff_mass(eff_mass_samples, fitted_mass, plot_filename):
    fig, ax = plt.subplots(layout="constrained")

    eff_mass_value = eff_mass_samples.mean
    num_timeslices = len(eff_mass_value)
    eff_mass_uncertainty = eff_mass_samples.std()

    ax.errorbar(
        np.arange(num_timeslices) + 1,
        eff_mass_value,
        eff_mass_uncertainty,
        ls="none",
        marker=".",
        label="Effective mass",
        color="C0",
    )
    ax.plot(
        [0, num_timeslices],
        [fitted_mass.nominal_value] * 2,
        color="C1",
        label="Fitted mass",
    )
    ax.fill_between(
        [0, num_timeslices],
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
    (ensemble,) = get_ensemble(
        data, beta=args.beta, mAS=args.mAS, Nt=args.Nt, Ns=args.Ns
    )
    eff_mass_samples = get_eff_mass_samples(
        ensemble, args.min_trajectory, args.max_trajectory, args.trajectory_step
    )
    fitted_mass_samples = eff_mass_samples[
        args.plateau_start : args.plateau_end
    ].weighted_mean()
    fitted_mass = bootstrap_finalize(fitted_mass_samples)

    if args.effmass_plot_file:
        plot_eff_mass(eff_mass_samples, fitted_mass, args.effmass_plot_file)

    metadata = {
        "ensemble_name": args.ensemble_name,
        "beta": args.beta,
        "mAS": args.mAS,
        "Nt": args.Nt,
        "Ns": args.Ns,
    }
    dump_dict({**metadata, "mPCAC": fitted_mass}, args.output_file_mean)
    if args.output_file_samples:
        dump_samples(
            {**metadata, "mPCAC": fitted_mass_samples},
            args.output_file_samples,
        )


if __name__ == "__main__":
    main()
