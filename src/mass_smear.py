#!/usr/bin/env python3


import h5py
import numpy as np

from .bootstrap import bootstrap_finalize, BOOTSTRAP_SAMPLE_COUNT, BootstrapSampleSet
from .dump import dump_dict, dump_samples
from .extract import extract_meson_mass
from .mass import (
    get_meson_corr,
    get_args,
)
from .read_hdf5 import get_ensemble
from .plots_common import plot_meson_smear



def ch_extraction(ensemble, args):

    c_sample = np.zeros(shape=(BOOTSTRAP_SAMPLE_COUNT, args.Nt))
    c_mean = np.zeros(shape=(1, args.Nt))

    corr = get_meson_corr(ensemble, args, int(args.n_smear_source), 0, args.channel)

    c_sample[:,:] = corr.samples
    c_mean[:,:] = corr.mean

    corr_set = BootstrapSampleSet(c_mean, c_sample)

    
    mass, matrix_element, chi2 = extract_meson_mass(
        corr_set, int(args.smear_plateau_start), int(args.smear_plateau_end)
    )

    if args.effmass_plot_file:

        plot_meson_smear(args, corr_set, mass)

    return mass, matrix_element, chi2


def main():
    args = get_args()

    data = h5py.File(args.h5file, "r")
    ensemble, = get_ensemble(
        data,
        beta=args.beta,
        mF=args.mF,
        Nt=args.Nt,
        Ns=args.Ns,
        epsilon=args.epsilon,
    )

    mass, matrix_element, chi2 = ch_extraction(ensemble, args)

    fitted_m = bootstrap_finalize(mass)

    metadata = {
        "ensemble_name": args.ensemble_name,
        "beta": args.beta,
        "mF": args.mF,
        "Nt": args.Nt,
        "Ns": args.Ns,
    }

    

    dump_dict(
        {
            **metadata,
            f"smear_{args.channel}_chisquare": chi2,
            f"smear_{args.channel}_mass": fitted_m,
        },
        args.output_file_mean,
    )
    if args.output_file_samples:
        dump_samples(
            {
                **metadata,
                f"smear_{args.channel}_mass_samples": mass.samples,
                f"smear_{args.channel}_mass_value": mass.mean,
            },
            args.output_file_samples,
        )


if __name__ == "__main__":
    main()
