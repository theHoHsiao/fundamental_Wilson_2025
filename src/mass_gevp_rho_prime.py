#!/usr/bin/env python3

from argparse import ArgumentParser, FileType
import h5py
import numpy as np
from itertools import product


from .bootstrap import BootstrapSampleSet, BOOTSTRAP_SAMPLE_COUNT, bootstrap_finalize
from .dump import dump_dict, dump_samples
from . import extract, fitting
from .mass_smear import bin_multi_source
from .read_hdf5 import get_ensemble


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
        help="Where to output the mean and uncertainty of m_rhoE1. (Defaults to stdout.)",
    )
    parser.add_argument(
        "--output_file_samples",
        type=FileType("w"),
        default=None,
        help="Where to output the bootstrap samples for m_rhoE1",
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
        "--num_source",
        type=int,
        default=None,
        help="number of source location used for smearing measurements",
    )
    parser.add_argument(
        "--GEVP_t0",
        type=int,
        default=None,
        help="number of source location used for smearing measurements",
    )
    return parser.parse_args()


def fold_correlators(C):
    return (C + np.roll(np.flip(C, axis=1), 1, axis=1)) / 2


def fold_correlators_cross(C):
    C_fold = (C - np.roll(np.flip(C, axis=1), 1, axis=1)) / 2

    C_fold[:, 0] = C[:, 0]

    return C_fold


def get_meson_Cmat_mix_N(ensemble, args, ch1, ch2):
    mixing_channel = [ch1, ch2]

    mat = np.zeros(shape=(BOOTSTRAP_SAMPLE_COUNT, args.Nt, 2, 2))
    mat_mean = np.zeros(shape=(1, args.Nt, 2, 2))

    matrix_size_channel = range(2)
    matrix_size_smearing = range(1)

    for a, b, i, j in product(
        matrix_size_channel,
        matrix_size_channel,
        matrix_size_smearing,
        matrix_size_smearing,
    ):
        if a == b:
            ch = mixing_channel[a]
            corr_set = bin_multi_source(ensemble, ch, args)
            mat[:, :, a + i, b + j] = fold_correlators(corr_set.samples)
            mat_mean[:, :, a + i, b + j] = fold_correlators(corr_set.mean)

        else:
            ch = mixing_channel[a] + "_" + mixing_channel[b] + "_re"
            corr_set = bin_multi_source(ensemble, ch, args)
            mat[:, :, a + i, b + j] = -fold_correlators_cross(corr_set.samples)
            mat_mean[:, :, a + i, b + j] = -fold_correlators_cross(corr_set.mean)

    return mat_mean, mat


def get_Cmat_VTmix(ensemble, args):
    target_channels = [
        ["g1", "g0g1"],
        ["g2", "g0g2"],
        ["g3", "g0g3"],
    ]

    bin_samples = []
    mean_bin = []
    for channels in target_channels:
        mean, samples = get_meson_Cmat_mix_N(ensemble, args, channels[0], channels[1])

        mean_bin.append(mean)
        bin_samples.append(samples)

    return np.array(mean_bin).mean(axis=0), np.array(bin_samples).mean(axis=0)


def main():
    args = get_args()
    if args.plateau_start == 0 and args.plateau_end == 0:
        mass, matrix_element, chi2 = (
            BootstrapSampleSet(np.nan, np.zeros(BOOTSTRAP_SAMPLE_COUNT) * np.nan),
            BootstrapSampleSet(np.nan, np.zeros(BOOTSTRAP_SAMPLE_COUNT) * np.nan),
            0,
        )

    else:
        data = h5py.File(args.h5file, "r")
        ensemble = get_ensemble(
            data,
            beta=args.beta,
            mAS=args.mAS,
            Nt=args.Nt,
            Ns=args.Ns,
            num_source=args.num_source,
            epsilon=args.epsilon,
        )

        Cmat_mean, Cmat = get_Cmat_VTmix(ensemble, args)

        eigenvalues = extract.GEVP_fixT(
            Cmat_mean, Cmat, args.GEVP_t0, args.GEVP_t0 + 1, args.Nt
        )

        mass, matrix_element, chi2 = fitting.fit_exp_bootstrap(
            eigenvalues[1], args.plateau_start, args.plateau_end
        )

    fitted_m = bootstrap_finalize(mass)
    fitted_a = bootstrap_finalize(matrix_element)

    metadata = {
        "ensemble_name": args.ensemble_name,
        "beta": args.beta,
        "mAS": args.mAS,
        "Nt": args.Nt,
        "Ns": args.Ns,
    }
    dump_dict(
        {
            **metadata,
            "smear_rhoE1_chisquare": chi2,
            "smear_rhoE1_mass": fitted_m,
            "smear_rhoE1_matrix_element": fitted_a,
        },
        args.output_file_mean,
    )
    if args.output_file_samples:
        dump_samples(
            {
                **metadata,
                "smear_rhoE1_mass_samples": mass.samples,
                "smear_rhoE1_mass_value": mass.mean,
                "smear_rhoE1_matrix_element_samples": matrix_element.samples,
                "smear_rhoE1_matrix_element_value": matrix_element.mean,
            },
            args.output_file_samples,
        )


if __name__ == "__main__":
    main()
