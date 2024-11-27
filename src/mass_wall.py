#!/usr/bin/env python3


import h5py
import numpy as np


from .bootstrap import BootstrapSampleSet, bootstrap_finalize
from .dump import dump_dict, dump_samples
from . import extract
from .mass import (
    get_correlator_samples,
    get_channel_tags,
    fold_correlators,
    get_args,
)
from .read_hdf5 import get_ensemble


def ps_extraction(ensemble, args):
    """
    This function packs two sets of correlators for the fitting:
        corr_aa refers to <A><A> exp(-mt)
        corr_ab refers to <A><B> exp(-mt)
    and returns the mass and the matrix element B for calculating decay constants
    """
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
    target_channels = get_channel_tags(args.channel)

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
    (ensemble,) = get_ensemble(
        data, beta=args.beta, mAS=args.mAS, Nt=args.Nt, Ns=args.Ns
    )

    if args.channel == "ps":
        mass, matrix_element, chi2 = ps_extraction(ensemble, args)

    else:
        mass, matrix_element, chi2 = ch_extraction(ensemble, args)

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
                f"{args.channel}_mass_samples": mass.samples,
                f"{args.channel}_mass_value": mass.mean,
                f"{args.channel}_matrix_element_samples": matrix_element.samples,
                f"{args.channel}_matrix_element_value": matrix_element.mean,
            },
            args.output_file_samples,
        )


if __name__ == "__main__":
    main()
