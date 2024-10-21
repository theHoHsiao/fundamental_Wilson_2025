#!/usr/bin/env python3

from argparse import ArgumentParser, FileType
import re
import h5py
import numpy as np
from uncertainties import ufloat

from .bootstrap import get_rng, sample_bootstrap_1d, BootstrapSampleSet
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
        raise ValueError("Did not uniquely identify one ensemble.")
    elif len(candidate_ensembles) == 0:
        raise ValueError("No ensembles found.")
    else:
        return candidate_ensembles


def fold_correlators(C):
    return (C + np.roll(np.flip(C, axis=1), 1, axis=1)) / 2


def bin_multi_source(ensemble, ch, args):
    tmp_bin = []
    tmp_bin_mean = []
    for n in range(args.num_source):
        tmp_bin.append(
            get_correlator_samples(
                ensemble[n],
                ch,
                args.N_sink,
                args.min_trajectory,
                args.max_trajectory,
                args.trajectory_step,
            ).samples
        )
        tmp_bin_mean.append(
            get_correlator_samples(
                ensemble[n],
                ch,
                args.N_sink,
                args.min_trajectory,
                args.max_trajectory,
                args.trajectory_step,
            ).mean
        )

    mean = np.zeros(shape=(1, args.Nt))
    mean[0] = np.array(tmp_bin_mean).mean(axis=0)

    return BootstrapSampleSet(mean, np.array(tmp_bin).mean(axis=0))


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

    # return sample_bootstrap_1d(C.T, get_rng(ensemble.name))
    return BootstrapSampleSet(
        C.mean(axis=1), sample_bootstrap_1d(C.T, get_rng(ensemble.name))
    )


def channel_tags(ch):
    return {
        "ps": ["g5"],
        "v": ["g1", "g2", "g3"],
        "t": ["g0g1", "g0g2", "g0g3"],
        "av": ["g5g1", "g5g2", "g5g3"],
        "at": ["g0g5g1", "g0g5g2", "g0g5g3"],
        "s": ["id"],
    }.get(ch, ch)


def ch_extraction(ensemble, args):
    CHs = channel_tags(args.channel)

    tmp_bin = []
    mean = []
    for j in range(len(CHs)):
        tmp_bin.append(
            fold_correlators(bin_multi_source(ensemble, CHs[j], args).samples)
            * args.Ns**3
        )
        mean.append(
            fold_correlators(bin_multi_source(ensemble, CHs[j], args).mean) * args.Ns**3
        )

    corr = BootstrapSampleSet(
        np.array(mean).mean(axis=0), np.array(tmp_bin).mean(axis=0)
    )

    # print(corr)

    m_tmp, a_tmp, chi2 = extract.meson_mass_sample(
        corr, args.plateau_start, args.plateau_end
    )

    return m_tmp, a_tmp, chi2


def main():
    args = get_args()
    # plt.style.use(args.plot_styles)

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

    m_tmp, a_tmp, chi2 = ch_extraction(ensemble, args)

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
            "chi_sqr_dof": chi2,
            f"smear_{args.channel}_mass": fitted_m,
            f"smear_{args.channel}_matrix_element": fitted_a,
        },
        args.output_file_mean,
    )
    if args.output_file_samples:
        dump_samples(
            {
                **metadata,
                f"smear_{args.channel}_mass_samples": m_tmp.samples,
                f"smear_{args.channel}_mass_value": m_tmp.mean,
                f"smear_{args.channel}_matrix_element_samples": a_tmp.samples,
                f"smear_{args.channel}_matrix_element_value": a_tmp.mean,
            },
            args.output_file_samples,
        )


if __name__ == "__main__":
    main()
