#!/usr/bin/env python3


import h5py
import numpy as np

from .bootstrap import bootstrap_finalize, BOOTSTRAP_SAMPLE_COUNT
from .dump import dump_dict, dump_samples
from . import fitting
from .mass import (
    get_meson_corr,
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
    corr_ss = get_meson_corr(ensemble, args, args.n_smear_max, args.n_smear_max, args.channel)
    
    rep = args.channel.split("_")[0]
    corr_sp = get_meson_corr(ensemble, args, args.n_smear_max, 0, rep+"_"+"ps-av")

    mass, matrix_element, chi2 = fitting.fit_coshsinh_simultaneous(
        corr_ss,
        corr_sp,
        args.E0_plateau_start,
        args.E0_plateau_end,
        args.Nt
    )

    return mass, matrix_element, chi2

def ch_extraction(ensemble, args):
    """
    This function packs two sets of correlators for the fitting:
        corr_aa refers to <A><A> exp(-mt)
        corr_ab refers to <A><B> exp(-mt)
    and returns the mass and the matrix element B for calculating decay constants
    """
    corr_ss = get_meson_corr(ensemble, args, args.n_smear_max, args.n_smear_max, args.channel)
    corr_sp = get_meson_corr(ensemble, args, args.n_smear_max, 0, args.channel)

    mass, matrix_element, chi2 = fitting.fit_cosh_simultaneous(
        corr_ss,
        corr_sp,
        args.E0_plateau_start,
        args.E0_plateau_end,
        args.Nt
    )
    

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
    
    meson = args.channel.split("_")[1]
    if meson == "ps":
        mass, matrix_element, chisquare = ps_extraction(ensemble, args)

    else:
        mass, matrix_element, chisquare = ch_extraction(ensemble, args)


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
            f"{args.channel}_chisquare": chisquare,
            f"{args.channel}_mass": bootstrap_finalize(mass),
            f"{args.channel}_matrix_element": bootstrap_finalize(matrix_element),
        },
        args.output_file_mean,
    )

    if args.output_file_samples:
        data_to_save = {**metadata}

        data_to_save[f"{args.channel}_chisquare"] = chisquare
        data_to_save[f"{args.channel}_mass"] = mass
        data_to_save[f"{args.channel}_matrix_element"] = matrix_element

        dump_samples(data_to_save,args.output_file_samples)


if __name__ == "__main__":
    main()
