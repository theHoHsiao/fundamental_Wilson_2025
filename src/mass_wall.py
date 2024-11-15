#!/usr/bin/env python3


import h5py
import numpy as np


from .bootstrap import BootstrapSampleSet, bootstrap_finalize
from .dump import dump_dict, dump_samples
from . import extract
from .mass import (
    get_correlators,
    get_correlator_samples,
    channel_tags,
    fold_correlators,
    get_args,
)


def ps_extraction(ensemble, args):
    corr_aa = get_correlator_samples(
        ensemble,
        "TRIPLET/g5",
        args.min_trajectory,
        args.max_trajectory,
        args.trajectory_step,
    )

    aa_mean = np.zeros(shape=(1, args.Nt))
    aa_mean[0] = corr_aa.mean * args.Ns**3
    C_aa = BootstrapSampleSet(aa_mean, corr_aa.samples * args.Ns**3)

    corr_ab = get_correlator_samples(
        ensemble,
        "TRIPLET/g5_g0g5_re",
        args.min_trajectory,
        args.max_trajectory,
        args.trajectory_step,
    )

    ab_mean = np.zeros(shape=(1, args.Nt))
    ab_mean[0] = corr_ab.mean * args.Ns**3
    C_ab = BootstrapSampleSet(ab_mean, corr_ab.samples * args.Ns**3)

    mass, matrix_element, chi2 = extract.meson_decay_constant(
        C_aa, C_ab, args.plateau_start, args.plateau_end
    )
    return mass, matrix_element, chi2


def ch_extraction(ensemble, args):
    target_channels = channel_tags(args.channel)

    bin_samples = []
    bin_mean = []
    for j in range(len(target_channels)):
        tmp_set = get_correlator_samples(
            ensemble,
            f"TRIPLET/{target_channels[j]}",
            args.min_trajectory,
            args.max_trajectory,
            args.trajectory_step,
        )

        bin_samples.append(tmp_set.samples * args.Ns**3)
        bin_mean.append(tmp_set.mean * args.Ns**3)

    mean = np.zeros(shape=(1, args.Nt))
    mean[0] = np.array(bin_mean).mean(axis=0)
    mean = fold_correlators(mean)
    samples = fold_correlators(np.array(bin_samples).mean(axis=0))

    corr = BootstrapSampleSet(mean, samples)

    mass, matrix_element, chi2 = extract.extract_meson_mass(
        corr, args.plateau_start, args.plateau_end
    )

    return mass, matrix_element, chi2


def main():
    args = get_args()

    data = h5py.File(args.h5file, "r")
    ensemble = get_correlators(
        data, beta=args.beta, mAS=args.mAS, Nt=args.Nt, Ns=args.Ns
    )[0]

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
            f"{args.channel}_chisquare": chi2,
            f"{args.channel}_mass": fitted_m,
            f"{args.channel}_matrix_element": fitted_a,
        },
        args.output_file_mean,
    )
    if args.output_file_samples:
        dump_samples(
            {
                **metadata,
                f"{args.channel}_mass_samples": m_tmp.samples,
                f"{args.channel}_mass_value": m_tmp.mean,
                f"{args.channel}_matrix_element_samples": a_tmp.samples,
                f"{args.channel}_matrix_element_value": a_tmp.mean,
            },
            args.output_file_samples,
        )


if __name__ == "__main__":
    main()
