#!/usr/bin/env python3

from argparse import ArgumentParser, FileType
import h5py
import numpy as np
from itertools import product

from src.plots_common import plot_meson_gevp_energy_states_only, plot_meson_smear



from .bootstrap import BootstrapSampleSet, BOOTSTRAP_SAMPLE_COUNT, bootstrap_finalize
from .dump import dump_dict, dump_samples
from . import extract, fitting
from .mass import get_meson_corr
from .read_hdf5 import get_ensemble


def get_args():
    parser = ArgumentParser(
        description="Compute the mass and matrix element from correlators in an HDF5 file"
    )

    parser.add_argument("h5file", help="The file to read")
    parser.add_argument(
        "--effmass_plot_file",
        default=None,
        help="Where to output the effective mass plot. (Skipped if not specified)",
    )
    parser.add_argument(
        "--plot_styles",
        default="styles/paperdraft.mplstyle",
        help="Stylesheet to use for plots",
    )
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
        "--mF",
        type=float,
        default=None,
        help="The fundamental fermion mass of the ensemble to analyse",
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
        help="Where to output the mean and uncertainty of m_rhoE1. (Defaults to stdout.)",
    )
    parser.add_argument(
        "--output_file_samples",
        type=FileType("w"),
        default=None,
        help="Where to output the bootstrap samples for m_rhoE1",
    )
    parser.add_argument(
        "--epsilon",
        type=float,
        default=None,
        help="Wuppertal smearing epsilon",
    )
    parser.add_argument(
        "--N_source",
        type=int,
        default=None,
        help="Source smearing level",
    )
    parser.add_argument(
        "--N_sink",
        type=int,
        default=None,
        help="Optimal smearing level",
    )
    parser.add_argument(
        "--n_smear_min",
        type=int,
        default=None,
        help="min smearing level",
    )
    parser.add_argument(
        "--n_smear_max",
        type=int,
        default=None,
        help="max smearing level",
    )
    parser.add_argument(
        "--n_smear_diff",
        type=int,
        default=None,
        help="diff smearing level",
    )
    parser.add_argument(
        "--num_source",
        type=int,
        default=1,
        help="number of source location used for smearing measurements",
    )
    parser.add_argument(
        "--gevp_t0",
        type=int,
        default=None,
        help="Fixed time slice for GEVP",
    )
    parser.add_argument(
        "--bin_size",
        type=int,
        default=1,
        help="Number of consecutive configurations to bin together",
    )
    return parser.parse_args()


def fold_correlators(C):
    return (C + np.roll(np.flip(C, axis=1), 1, axis=1)) / 2


def fold_correlators_cross(C):
    C_fold = (C - np.roll(np.flip(C, axis=1), 1, axis=1)) / 2

    C_fold[:, 0] = C[:, 0]

    return C_fold


def get_meson_Cmat_mix_N(ensemble, args, ch1, ch2):
    mixing_channel = [ch1, ch2]

    smear_Ns = [str(i) for i in np.arange(args.n_smear_min, args.n_smear_max + 1, args.n_smear_diff)]
    N_smear_levels = len(smear_Ns)
    
    matrix_size_channel = range(2)
    matrix_size_smearing = range(N_smear_levels)

    m_size = len(matrix_size_channel) * N_smear_levels

    mat = np.zeros(shape=(BOOTSTRAP_SAMPLE_COUNT, args.Nt, m_size, m_size))
    mat_mean = np.zeros(shape=(1, args.Nt, m_size, m_size))

    for a, b, i, j in product(
        matrix_size_channel,
        matrix_size_channel,
        matrix_size_smearing,
        matrix_size_smearing,
    ):
        if a == b:
            ch = "f_"+mixing_channel[a]
            corr_set = get_meson_corr(
                ensemble, args, smear_Ns[i], smear_Ns[j], ch
                )
            mat[:, :, a + i*2, b + j*2] = corr_set.samples
            mat_mean[:, :, a + i*2, b + j*2] = corr_set.mean

        else:
            ch = "f_"+mixing_channel[a] + "-" + mixing_channel[b]
            corr_set = get_meson_corr(
                ensemble, args, smear_Ns[i], smear_Ns[j], ch
                )
            mat[:, :, a + i*2, b + j*2] = -corr_set.samples
            mat_mean[:, :, a + i*2, b + j*2] = -corr_set.mean
    
    #print(mat_mean[0,2,:,:])
    return mat_mean, mat


def main():
    args = get_args()
    if args.plateau_start == 0 and args.plateau_end == 0:
        mass, matrix_element, chi2 = (
            BootstrapSampleSet(np.nan, np.zeros(BOOTSTRAP_SAMPLE_COUNT) * np.nan),
            BootstrapSampleSet(np.nan, np.zeros(BOOTSTRAP_SAMPLE_COUNT) * np.nan),
            0,
        )

    else:
        data = h5py.File(args.h5file, "r")
        ensemble = get_ensemble(
            data,
            beta=args.beta,
            mF=args.mF,
            Nt=args.Nt,
            Ns=args.Ns,
            num_source=args.num_source,
            epsilon=args.epsilon,
        )[0]
        #print(ensemble[0]["source_N0_sink_N120"].keys())
        Cmat_mean, Cmat = get_meson_Cmat_mix_N(ensemble, args, "v", "t")

        eigenvalues = extract.gevp_fixT(Cmat_mean, Cmat, args.gevp_t0, args.gevp_t0+1, args.Nt/2)
       

        mass, matrix_element, chi2 = fitting.fit_exp_bootstrap(
            eigenvalues[1], args.plateau_start, args.plateau_end
        )
        print(f"chi2 = {chi2}")
        if args.effmass_plot_file is not None:
            plot_meson_gevp_energy_states_only(args, eigenvalues)
            #plot_meson_smear(args, eigenvalues[1], mass, args.plateau_start, args.plateau_end)

            

            
        
        

    fitted_m = bootstrap_finalize(mass)
    fitted_a = bootstrap_finalize(matrix_element)

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
            "gevp_f_rhoprime_chisquare": chi2,
            "gevp_f_rhoprime_E0_mass": fitted_m,
            "gevp_f_rhoprime_matrix_element": fitted_a,
        },
        args.output_file_mean,
    )
    if args.output_file_samples:
        dump_samples(
            {
                **metadata,
                "gevp_f_rhoprime_mass_samples": mass.samples,
                "gevp_f_rhoprime_mass_value": mass.mean,
                "gevp_f_rhoprime_matrix_element_samples": matrix_element.samples,
                "gevp_f_rhoprime_matrix_element_value": matrix_element.mean,
            },
            args.output_file_samples,
        )


if __name__ == "__main__":
    main()
