#!/usr/bin/env python3


import h5py
import numpy as np
import re

from .dump import dump_dict
from flow_analysis.stats.autocorrelation import integrated_autocorrelation_time
from .mass import get_args
from .read_hdf5 import get_ensemble, filter_configurations
from .utils import get_index_separation
from .bootstrap import get_rng, sample_bootstrap_0d
import matplotlib.pyplot as plt


def bin_data(values, bin_size):
    """
    Bin data into blocks of size bin_size by averaging within each block.
    Any leftover data at the end that doesn't fit into a full block is dropped.

    Parameters
    ----------
    values : array_like
        Input data (1D or 2D).
    bin_size : int
        Number of consecutive samples per bin.

    Returns
    -------
    binned_values : ndarray
        Binned data with shape (n_bins, ...) where n_bins = len(values)//bin_size
    """
    values = np.asarray(values)
    n = len(values)
    n_bins = n // bin_size
    trimmed = values[: n_bins * bin_size]

    if values.ndim == 1:
        return trimmed.reshape(n_bins, bin_size).mean(axis=1)
    elif values.ndim == 2:
        return trimmed.reshape(n_bins, bin_size, values.shape[1]).mean(axis=1)
    else:
        raise ValueError("bin_data only supports 1D or 2D arrays")

def check_bin_sizes(args, ensemble, corr_ps_point, corr_ps_smear, max_size=10):
    plt.style.use(args.plot_styles)
    fig, ax = plt.subplots(layout="constrained")
    
    ax2 = ax.twinx()
    
    bootsample_point = sample_bootstrap_0d(corr_ps_point, get_rng(ensemble.name))
    bootsample_smear = sample_bootstrap_0d(corr_ps_smear, get_rng(ensemble.name))
    ax.errorbar(1-0.2, bootsample_point.mean, yerr=bootsample_point.samples.std(), fmt='o', label='Point Source', color='C0')
    ax2.errorbar(1+0.2, bootsample_smear.mean, yerr=bootsample_smear.samples.std(), fmt='o', label='Smeared Source', color='C1')
    
    for bin_size in range(2, max_size+1):
        bined_smear = bin_data(corr_ps_smear, bin_size)
        bined_point = bin_data(corr_ps_point, bin_size)

        bootsample_point = sample_bootstrap_0d(bined_point, get_rng(ensemble.name))
        bootsample_smear = sample_bootstrap_0d(bined_smear, get_rng(ensemble.name))
        ax.errorbar(bin_size-0.2, bootsample_point.mean, yerr=bootsample_point.samples.std(), fmt='o', color='C0')
        ax2.errorbar(bin_size+0.2, bootsample_smear.mean, yerr=bootsample_smear.samples.std(), fmt='o', color='C1')
    
    ax.set_xlabel("Bin Size")
    ax.set_ylabel("Point-Point PS Correlator")
    ax2.set_ylabel("Smeares-Smeared PS Correlator")

    ax.tick_params(axis="y", colors="C0")
    ax2.tick_params(axis="y", colors="C1")
    
    ax.set_xticks(np.arange(1, max_size+1))

    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    fig.legend(lines1 + lines2, labels1 + labels2, loc="outside upper center", ncol=2)


    fig.savefig(args.effmass_plot_file)

   
        

def ps_correlator_autocorrelation(ensemble, args):
    indices = np.asarray(
        [
            int(re.match(".*n([0-9]+)$", filename.decode()).groups()[0])
            for filename in ensemble["configurations"]
        ]
    )
    
    filtered_indices = filter_configurations(ensemble, args.min_trajectory, args.max_trajectory, args.trajectory_step)
    Ncfg = sum(filtered_indices)

    corr_ps_smear = ensemble[f"source_N{args.n_smear_max}_sink_N{args.n_smear_max}/TRIPLET g5"][args.E0_plateau_start, filtered_indices]
    corr_ps_point = ensemble[f"source_N0_sink_N0/TRIPLET g5"][args.E0_plateau_start, filtered_indices]
    
    trajectory_separation = get_index_separation(indices[filtered_indices])

    tau_ps_correlator_smear = integrated_autocorrelation_time(corr_ps_smear) * trajectory_separation
    tau_ps_correlator_point = integrated_autocorrelation_time(corr_ps_point) * trajectory_separation


    if args.effmass_plot_file:
        check_bin_sizes(args, ensemble, corr_ps_point, corr_ps_smear)

    
    return Ncfg, tau_ps_correlator_smear, tau_ps_correlator_point, trajectory_separation


def main():
    args = get_args()

    data = h5py.File(args.h5file, "r")
    (ensemble,) = get_ensemble(
        data, beta=args.beta, mF=args.mF, Nt=args.Nt, Ns=args.Ns
    )

    Ncfg, auto_smear, auto_point, trajectory_separation = ps_correlator_autocorrelation(ensemble, args)

    metadata = {
        "ensemble_name": args.ensemble_name,
        "beta": args.beta,
        "mAS": args.mAS,
        "Nt": args.Nt,
        "Ns": args.Ns,
        "Ncfg": Ncfg,
    }
    dump_dict(
        {**metadata,
         "tau_exp_ps_correlator_smear": auto_smear,
         "tau_exp_ps_correlator_point": auto_point,
         "delta_traj_auto": trajectory_separation},
        args.output_file_mean,
    )


if __name__ == "__main__":
    main()
