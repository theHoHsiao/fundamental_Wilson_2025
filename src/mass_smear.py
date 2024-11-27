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


def bin_multi_source(ensemble, ch, args):
    bin_samples = []
    bin_mean = []
    for source_location in ensemble:
        tmp_set = get_correlator_samples(
            source_location,
            f"source_N100_sink_N{args.N_sink}/TRIPLET {ch}",
            args.min_trajectory,
            args.max_trajectory,
            args.trajectory_step,
        )

        bin_samples.append(tmp_set.samples)
        bin_mean.append(tmp_set.mean)

    mean = np.zeros(shape=(1, args.Nt))
    mean[0] = np.array(bin_mean).mean(axis=0)

    return BootstrapSampleSet(mean, np.array(bin_samples).mean(axis=0))


def ch_extraction(ensemble, args):
    target_channels = get_channel_tags(args.channel)

    bin_samples = []
    bin_mean = []
    for j in range(len(target_channels)):
        tmp_set = bin_multi_source(ensemble, target_channels[j], args)
        bin_samples.append(fold_correlators(tmp_set.samples) * args.Ns**3)
        bin_mean.append(fold_correlators(tmp_set.mean) * args.Ns**3)

    corr = BootstrapSampleSet(
        np.array(bin_mean).mean(axis=0), np.array(bin_samples).mean(axis=0)
    )

    mass, matrix_element, chi2 = extract.extract_meson_mass(
        corr, args.plateau_start, args.plateau_end
    )

    return mass, matrix_element, chi2


def main():
    args = get_args()

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
                f"smear_{args.channel}_mass_samples": mass.samples,
                f"smear_{args.channel}_mass_value": mass.mean,
                f"smear_{args.channel}_matrix_element_samples": matrix_element.samples,
                f"smear_{args.channel}_matrix_element_value": matrix_element.mean,
            },
            args.output_file_samples,
        )


if __name__ == "__main__":
    main()
