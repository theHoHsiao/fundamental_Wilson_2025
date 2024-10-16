#!/usr/bin/env python3

from argparse import ArgumentParser, FileType
import re
import h5py
import numpy as np

from .bootstrap import (
    get_rng,
    sample_bootstrap_1d,
    bootstrap_finalize,
    BOOTSTRAP_SAMPLE_COUNT,
)
from .dump import dump_dict, dump_samples
from . import extract, fitting


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


def get_correlators(
    ensembles, beta=None, mAS=None, Nt=None, Ns=None, num_source=None, epsilon=None
):
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
        if epsilon is not None and ensemble.get("Wuppertal_eps_anti", [])[0] != epsilon:
            continue
        candidate_ensembles.append(ensemble)
    if num_source is not None and len(candidate_ensembles) != num_source:
        print(candidate_ensembles, ensemble.get("Wuppertal_eps_anti", [])[0])
        raise ValueError("Did not uniquely identify one ensemble.")
    elif len(candidate_ensembles) == 0:
        raise ValueError("No ensembles found.")
    else:
        return candidate_ensembles


def fold_correlators(C):
    return (C + np.roll(np.flip(C, axis=1), 1, axis=1)) / 2


def fold_correlators_cross(C):
    C_fold = (C - np.roll(np.flip(C, axis=1), 1, axis=1)) / 2

    C_fold[:, 0] = C[:, 0]

    return C_fold


def bin_multi_source(ensemble, ch, args):
    tmp_bin = []
    for n in range(args.num_source):
        tmp_bin.append(
            get_correlator_samples(
                ensemble[n],
                ch,
                args.N_sink,
                args.min_trajectory,
                args.max_trajectory,
                args.trajectory_step,
            )
        )
    return np.array(tmp_bin).mean(axis=0)


def get_correlator_samples(
    ensemble,
    ch,
    N_sink,
    min_trajectory=None,
    max_trajectory=None,
    trajectory_step=1,
):
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

    C = ensemble[f"source_N100_sink_N{N_sink}"][f"TRIPLET {ch}"][:, filtered_indices]

    return sample_bootstrap_1d(C.T, get_rng(ensemble.name))


def get_meson_Cmat_mix_N(ensemble, args, ch1, ch2):
    mixing_channel = [ch1, ch2]

    mat = np.zeros(shape=(BOOTSTRAP_SAMPLE_COUNT, args.Nt, 2, 2))

    for a in range(2):
        for b in range(2):
            for i in range(1):
                for j in range(1):
                    if a == b:
                        ch = mixing_channel[a]
                        mat[:, :, a + i, b + j] = fold_correlators(
                            bin_multi_source(ensemble, ch, args)
                        )

                    else:
                        ch = mixing_channel[a] + "_" + mixing_channel[b] + "_re"
                        mat[:, :, a + i, b + j] = -fold_correlators_cross(
                            bin_multi_source(ensemble, ch, args)
                        )

    return mat


def get_Cmat_VTmix(ensemble, args):
    CHs = [
        ["g1", "g0g1"],
        ["g2", "g0g2"],
        ["g3", "g0g3"],
    ]

    tmp_bin = []
    for i in range(len(CHs)):
        ch = CHs[i]

        tmp_bin.append(get_meson_Cmat_mix_N(ensemble, args, ch[0], ch[1]))

    return np.array(tmp_bin).mean(axis=0)


def main():
    args = get_args()
    # plt.style.use(args.plot_styles)

    if args.plateau_start == 0 and args.plateau_end == 0:
        print("no plateau to fit...")
        m_tmp, a_tmp, chi2 = (
            np.zeros(BOOTSTRAP_SAMPLE_COUNT) * np.nan,
            np.zeros(BOOTSTRAP_SAMPLE_COUNT) * np.nan,
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

        Cmat = get_Cmat_VTmix(ensemble, args)

        LAM, VEC = extract.GEVP_fixT(Cmat, args.GEVP_t0, args.GEVP_t0 + 1, args.Nt)

        m_tmp, a_tmp, chi2 = fitting.fit_exp_booerr(
            LAM[:, :, 1], args.plateau_start, args.plateau_end
        )

    fitted_m = bootstrap_finalize(m_tmp)
    fitted_a = bootstrap_finalize(a_tmp)

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
            "chi_sqr_dof": chi2,
            "smear_rhoE1_mass": fitted_m,
            "smear_rhoE1_matrix_element": fitted_a,
        },
        args.output_file_mean,
    )
    if args.output_file_samples:
        dump_samples(
            {
                **metadata,
                "smear_rhoE1_mass_samples": m_tmp,
                "smear_rhoE1_mass_value": fitted_m.nominal_value,
                "smear_rhoE1_matrix_element_samples": a_tmp,
                "smear_rhoE1_matrix_element_value": fitted_a.nominal_value,
            },
            args.output_file_samples,
        )


if __name__ == "__main__":
    main()
