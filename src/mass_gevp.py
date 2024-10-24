#!/usr/bin/env python3

from argparse import ArgumentParser, FileType
import h5py
import numpy as np
from uncertainties import ufloat

from .bootstrap import BootstrapSampleSet, BOOTSTRAP_SAMPLE_COUNT
from .dump import dump_dict, dump_samples
from . import extract, fitting
from .mass_smear import bin_multi_source, get_correlators


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
        help="Where to output the mean and uncertainty of mPCAC. (Defaults to stdout.)",
    )
    parser.add_argument(
        "--output_file_samples",
        type=FileType("w"),
        default=None,
        help="Where to output the bootstrap samples for mPCAC",
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

    for a in range(2):
        for b in range(2):
            for i in range(1):
                for j in range(1):
                    if a == b:
                        ch = mixing_channel[a]
                        corr_set = bin_multi_source(ensemble, ch, args)
                        mat[:, :, a + i, b + j] = fold_correlators(corr_set.samples)
                        mat_mean[:, :, a + i, b + j] = fold_correlators(corr_set.mean)

                    else:
                        ch = mixing_channel[a] + "_" + mixing_channel[b] + "_re"
                        corr_set = bin_multi_source(ensemble, ch, args)
                        mat[:, :, a + i, b + j] = -fold_correlators_cross(
                            corr_set.samples
                        )
                        mat_mean[:, :, a + i, b + j] = -fold_correlators_cross(
                            corr_set.mean
                        )

    return mat_mean, mat


def get_Cmat_VTmix(ensemble, args):
    CHs = [
        ["g1", "g0g1"],
        ["g2", "g0g2"],
        ["g3", "g0g3"],
    ]

    tmp_bin = []
    mean_bin = []
    for i in range(len(CHs)):
        ch = CHs[i]

        mean, samples = get_meson_Cmat_mix_N(ensemble, args, ch[0], ch[1])

        mean_bin.append(mean)
        tmp_bin.append(samples)

    return np.array(mean_bin).mean(axis=0), np.array(tmp_bin).mean(axis=0)


def main():
    args = get_args()
    # plt.style.use(args.plot_styles)

    if args.plateau_start == 0 and args.plateau_end == 0:
        print("no plateau to fit...")
        m_tmp, a_tmp, chi2 = (
            BootstrapSampleSet(np.nan, np.zeros(BOOTSTRAP_SAMPLE_COUNT) * np.nan),
            BootstrapSampleSet(np.nan, np.zeros(BOOTSTRAP_SAMPLE_COUNT) * np.nan),
            0,
        )

    else:
        data = h5py.File(args.h5file, "r")
        ensemble = get_correlators(
            data,
            beta=args.beta,
            mAS=args.mAS,
            Nt=args.Nt,
            Ns=args.Ns,
            num_source=args.num_source,
            epsilon=args.epsilon,
        )

        Cmat_mean, Cmat = get_Cmat_VTmix(ensemble, args)

        LAM = extract.GEVP_fixT(
            Cmat_mean, Cmat, args.GEVP_t0, args.GEVP_t0 + 1, args.Nt
        )

        E_mean, A_mean, chi2, E_samples, A_samples = fitting.fit_exp_booerr(
            LAM[1], args.plateau_start, args.plateau_end
        )
        m_tmp = BootstrapSampleSet(E_mean, E_samples)
        a_tmp = BootstrapSampleSet(
            A_mean / np.sqrt(E_mean), A_samples / np.sqrt(E_samples)
        )

    fitted_m = ufloat(m_tmp.mean, m_tmp.samples.std())
    fitted_a = ufloat(a_tmp.mean, a_tmp.samples.std())

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
                "smear_rhoE1_mass_samples": m_tmp.samples,
                "smear_rhoE1_mass_value": m_tmp.mean,
                "smear_rhoE1_matrix_element_samples": a_tmp.samples,
                "smear_rhoE1_matrix_element_value": a_tmp.mean,
            },
            args.output_file_samples,
        )


if __name__ == "__main__":
    main()
