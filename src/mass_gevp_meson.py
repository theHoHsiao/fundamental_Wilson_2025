#!/usr/bin/env python3


import h5py
import numpy as np

from .bootstrap import BootstrapSampleSet, BOOTSTRAP_SAMPLE_COUNT, bootstrap_finalize
from .dump import dump_dict, dump_samples
from . import extract
from .mass import (
    get_meson_corr,
    get_args,
)
from .read_hdf5 import get_ensemble, filter_configurations
from .plots_common import plot_meson_gevp_energy_states

def get_meson_Cmat_single(ensemble, args, Nmin, Nmax, Nd, channel):
    source_Ns = [str(i) for i in np.arange(Nmin, Nmax + 1, Nd)]
    sink_Ns = [str(i) for i in np.arange(Nmin, Nmax + 1, Nd)]

    size = len(source_Ns)
    

    mat_sample = np.zeros(shape=(BOOTSTRAP_SAMPLE_COUNT, args.Nt, size, size))
    mat_mean = np.zeros(shape=(1, args.Nt, size, size))

    # print(mat.shape)

    for i in range(size):
        for j in range(size):

            corr = get_meson_corr(
                ensemble, args, source_Ns[i], sink_Ns[j], channel
                )
            mat_sample[:, :, i, j] = corr.samples
            mat_mean[:, :, i, j] = corr.mean


    return BootstrapSampleSet(mat_mean, mat_sample)


def gevp_meson_extraction(ensemble, args):

    corr_mat = get_meson_Cmat_single(ensemble, args, args.n_smear_min, args.n_smear_max, args.n_smear_diff, args.channel)
    eigenvalues = extract.gevp_fixT(corr_mat.mean, corr_mat.samples, args.gevp_t0, args.gevp_t0+1, args.Nt)

    return eigenvalues


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

    eigenvalues = gevp_meson_extraction(ensemble, args)
    masses, chisquares = extract.extract_energy_states(eigenvalues, args)

    filtered_indices = filter_configurations(ensemble, args.min_trajectory, args.max_trajectory, args.trajectory_step)
    N_cnfg = np.sum(filtered_indices)
    #print(filtered_indices)

    if args.effmass_plot_file:
        plot_meson_gevp_energy_states(args, eigenvalues, masses)

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
            f"gevp_{args.channel}_E0_chisquare": chisquares[0],
            f"gevp_{args.channel}_E0_mass": bootstrap_finalize(masses[0]),
            "N_cnfg": N_cnfg,
            "delta_traj_spec": args.trajectory_step,
        },
        args.output_file_mean,
    )

    if args.output_file_samples:
        data_to_save = {**metadata}

        for n, mass in enumerate(masses):
            data_to_save[f"gevp_{args.channel}_E{n}_chisquare"] = chisquares[n]
            data_to_save[f"gevp_{args.channel}_E{n}_mass"] = mass
        
        dump_samples(data_to_save,args.output_file_samples)


if __name__ == "__main__":
    main()
