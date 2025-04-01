#!/usr/bin/env python3

from argparse import ArgumentParser, FileType
import numpy as np

from .bootstrap import get_rng, sample_bootstrap_1d
from .read_hdf5 import filter_configurations, get_meson_h5_representation


def get_args():
    parser = ArgumentParser(
        description="Compute the mass and matrix element from correlators in an HDF5 file"
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
        "--mF",
        type=float,
        default=None,
        help="The fundamental fermion mass of the ensemble to analyse",
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
        "--gevp_t0",
        type=int,
        default=None,
        help="Fixed time slice for GEVP",
    )
    parser.add_argument(
        "--E0_plateau_start",
        type=int,
        default=None,
        help="Time slice at which plateau starts",
    )
    parser.add_argument(
        "--E0_plateau_end", type=int, default=None, help="Time slice at which plateau ends"
    )
    parser.add_argument(
        "--E1_plateau_start",
        type=int,
        default=None,
        help="Time slice at which plateau starts",
    )
    parser.add_argument(
        "--E1_plateau_end", type=int, default=None, help="Time slice at which plateau ends"
    )
    parser.add_argument(
        "--E2_plateau_start",
        type=int,
        default=None,
        help="Time slice at which plateau starts",
    )
    parser.add_argument(
        "--E2_plateau_end", type=int, default=None, help="Time slice at which plateau ends"
    )
    parser.add_argument(
        "--E3_plateau_start",
        type=int,
        default=None,
        help="Time slice at which plateau starts",
    )
    parser.add_argument(
        "--E3_plateau_end", type=int, default=None, help="Time slice at which plateau ends"
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
        help="Where to output the mean and uncertainty of meson mass. (Defaults to stdout.)",
    )
    parser.add_argument(
        "--output_file_samples",
        type=FileType("w"),
        default=None,
        help="Where to output the bootstrap samples for meson mass",
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
    parser.add_argument(
        "--channel",
        choices=["f_ps", "f_v", "f_t", "f_av", "f_at", "f_s",
                 "as_ps", "as_v", "as_t", "as_av", "as_at", "as_s",
                 "lambda_even", "sigma_even", "sigmastar_even",
                 "lambda_odd", "sigma_odd", "sigmastar_odd"],
        default=None,
        help="Measuring channel",
    )
    parser.add_argument(
        "--epsilon",
        type=float,
        default=None,
        help="Wuppertal smearing epsilon",
    )
    parser.add_argument(
        "--N_sink",
        type=int,
        default=None,
        help="Optimal smearing level",
    )
    parser.add_argument(
        "--n_smear_min",
        type=int,
        default=None,
        help="min smearing level",
    )
    parser.add_argument(
        "--n_smear_max",
        type=int,
        default=None,
        help="max smearing level",
    )
    parser.add_argument(
        "--n_smear_diff",
        type=int,
        default=None,
        help="diff smearing level",
    )
    parser.add_argument(
        "--num_source",
        type=int,
        default=None,
        help="number of source location used for smearing measurements",
    )
    return parser.parse_args()


def renormalisation_constant(ch):
    return {
        "v": -20.57,
        "av": -15.82,
        "ps": -15.82,
    }.get(ch, ch)


def get_correlator_samples(
    ensemble,
    measurement,
    min_trajectory=None,
    max_trajectory=None,
    trajectory_step=1,
):
    filtered_indices = filter_configurations(
        ensemble, min_trajectory, max_trajectory, trajectory_step
    )

    C = ensemble[measurement][:, filtered_indices]

    return sample_bootstrap_1d(C.T, get_rng(ensemble.name))


def get_channel_tags(ch):
    return {
        "ps": ["g5"],
        "v": ["g1", "g2", "g3"],
        "t": ["g0g1", "g0g2", "g0g3"],
        "av": ["g5g1", "g5g2", "g5g3"],
        "at": ["g0g5g1", "g0g5g2", "g0g5g3"],
        "s": ["id"],
    }.get(ch, ch)


def fold_correlators(C):
    return (C + np.roll(np.flip(C, axis=1), 1, axis=1)) / 2


def fold_correlators_cross(C):
    C_fold = (C - np.roll(np.flip(C, axis=1), 1, axis=1)) / 2

    C_fold[:, 0] = C[:, 0]

    return C_fold


def bin_meson_correlator_samples(
    ensemble,
    measurement,
    Nsource,
    Nsink,
    min_trajectory=None,
    max_trajectory=None,
    trajectory_step=1,
):
    filtered_indices = filter_configurations(ensemble, min_trajectory, max_trajectory, trajectory_step)

    rep = get_meson_h5_representation(measurement.split("_")[0])
    target_channels = get_channel_tags(measurement.split("_")[1])
    C_bin =[]
    for channel in target_channels:

        #C = ensemble[f"source_N{Nsource}_sink_N{Nsink}/{rep} {channel}"][:, filtered_indices]
        C = ensemble[f"source_N{Nsource}_sink_N{Nsink}/{rep} {channel}"][:, :]
        # TO DO: how shall we deal with jumpping configs

        C_bin.append(C)
    
    C_bin = np.array(C_bin)
    C = C_bin.mean(axis=0)
    
    if target_channels[0] == "g5_g0g5_re":
        C_flod = -fold_correlators_cross(C.T)
    else:
        C_flod = fold_correlators(C.T)

    return sample_bootstrap_1d(C_flod, get_rng(ensemble.name))

def get_meson_corr( ensemble, args, Nsource, Nsink, channel):

    try:
        corr = bin_meson_correlator_samples(
            ensemble,
            channel,
            Nsource,
            Nsink,
            args.min_trajectory,
            args.max_trajectory,
            args.trajectory_step,
        )
        
        return corr * args.Ns**3

    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return None
    