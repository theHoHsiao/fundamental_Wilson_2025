#!/usr/bin/env python3

from argparse import ArgumentParser, FileType
import re
import h5py
import numpy as np

from .bootstrap import get_rng, sample_bootstrap_1d, bootstrap_finalize
from .dump import dump_dict, dump_samples
from . import extract


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
        "--channel",
        choices=["ps", "v", "t", "av", "at", "s"],
        default=None,
        help="Measuring channel",
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


def get_correlator_samples(
    ensemble,
    ch,
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

    C = ensemble["TRIPLET"][f"{ch}"][:, filtered_indices]

    return sample_bootstrap_1d(C.T, get_rng(ensemble.name))


def channel_tags(ch):
    return {
        "v": ["g1", "g2", "g3"],
        "t": ["g0g1", "g0g2", "g0g3"],
        "av": ["g5g1", "g5g2", "g5g3"],
        "at": ["g0g5g1", "g0g5g2", "g0g5g3"],
        "s": ["id"],
    }.get(ch, ch)


def ps_extraction(ensemble, args):
    corr_aa = (
        get_correlator_samples(
            ensemble,
            "g5",
            args.min_trajectory,
            args.max_trajectory,
            args.trajectory_step,
        )
        * args.Ns**3
    )

    corr_ab = (
        get_correlator_samples(
            ensemble,
            "g5_g0g5_re",
            args.min_trajectory,
            args.max_trajectory,
            args.trajectory_step,
        )
        * args.Ns**3
    )

    m_tmp, a_tmp, chi2 = extract.meson_decay_sample(
        corr_aa, corr_ab, args.plateau_start, args.plateau_end
    )
    return m_tmp, a_tmp, chi2


def ch_extraction(ensemble, args):
    CHs = channel_tags(args.channel)

    tmp_bin = []
    for j in range(len(CHs)):
        tmp_bin.append(
            get_correlator_samples(
                ensemble,
                CHs[j],
                args.min_trajectory,
                args.max_trajectory,
                args.trajectory_step,
            )
        )

    corr = np.array(tmp_bin).mean(axis=0) * args.Ns**3

    m_tmp, a_tmp, chi2 = extract.meson_mass_sample(
        corr, args.plateau_start, args.plateau_end
    )

    return m_tmp, a_tmp, chi2


def main():
    args = get_args()
    # plt.style.use(args.plot_styles)

    data = h5py.File(args.h5file, "r")
    ensemble = get_correlators(
        data, beta=args.beta, mAS=args.mAS, Nt=args.Nt, Ns=args.Ns
    )

    if args.channel == "ps":
        m_tmp, a_tmp, chi2 = ps_extraction(ensemble, args)

    else:
        m_tmp, a_tmp, chi2 = ch_extraction(ensemble, args)

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
            f"{args.channel}_mass": fitted_m,
            f"{args.channel}_matrix_element": fitted_a,
        },
        args.output_file_mean,
    )
    if args.output_file_samples:
        dump_samples(
            {
                **metadata,
                f"{args.channel}_mass_samples": m_tmp,
                f"{args.channel}_matrix_element_samples": a_tmp,
            },
            args.output_file_samples,
        )


if __name__ == "__main__":
    main()
