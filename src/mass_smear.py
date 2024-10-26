#!/usr/bin/env python3

from argparse import ArgumentParser, FileType
import h5py
import numpy as np
from uncertainties import ufloat

from .bootstrap import BootstrapSampleSet
from .dump import dump_dict, dump_samples
from . import extract
from .mass import (
    get_correlators,
    get_correlator_samples,
    channel_tags,
    fold_correlators,
)


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


def bin_multi_source(ensemble, ch, args):
    tmp_bin = []
    tmp_bin_mean = []
    for n in range(args.num_source):
        tmp_set = get_correlator_samples(
            ensemble[n],
            f"source_N100_sink_N{args.N_sink}/TRIPLET {ch}",
            args.min_trajectory,
            args.max_trajectory,
            args.trajectory_step,
        )

        tmp_bin.append(tmp_set.samples)
        tmp_bin_mean.append(tmp_set.mean)

    mean = np.zeros(shape=(1, args.Nt))
    mean[0] = np.array(tmp_bin_mean).mean(axis=0)

    return BootstrapSampleSet(mean, np.array(tmp_bin).mean(axis=0))


def ch_extraction(ensemble, args):
    CHs = channel_tags(args.channel)

    tmp_bin = []
    mean = []
    for j in range(len(CHs)):
        tmp_set = bin_multi_source(ensemble, CHs[j], args)
        tmp_bin.append(fold_correlators(tmp_set.samples) * args.Ns**3)
        mean.append(fold_correlators(tmp_set.mean) * args.Ns**3)

    corr = BootstrapSampleSet(
        np.array(mean).mean(axis=0), np.array(tmp_bin).mean(axis=0)
    )

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
            f"smear_{args.channel}_chisquare": chi2,
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
