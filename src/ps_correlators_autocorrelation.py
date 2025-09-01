#!/usr/bin/env python3


import h5py
import numpy as np
import re

from .dump import dump_dict
from flow_analysis.stats.autocorrelation import exp_autocorrelation_fit
from .mass import get_args
from .read_hdf5 import get_ensemble, filter_configurations
from .utils import get_index_separation


def ps_correlator_autocorrelation(ensemble, args):
    indices = np.asarray(
        [
            int(re.match(".*n([0-9]+)$", filename.decode()).groups()[0])
            for filename in ensemble["configurations"]
        ]
    )
    
    filtered_indices = filter_configurations(ensemble, args.min_trajectory, args.max_trajectory, args.trajectory_step)

    corr_ps_smear = ensemble[f"source_N{args.n_smear_max}_sink_N{args.n_smear_max}/TRIPLET g5"][args.E0_plateau_start, filtered_indices]
    corr_ps_point = ensemble[f"source_N0_sink_N0/TRIPLET g5"][args.E0_plateau_start, filtered_indices]
    
    trajectory_separation = get_index_separation(indices[filtered_indices])

    tau_ps_correlator_smear = exp_autocorrelation_fit(corr_ps_smear) * trajectory_separation
    tau_ps_correlator_point = exp_autocorrelation_fit(corr_ps_point) * trajectory_separation

    return tau_ps_correlator_smear, tau_ps_correlator_point


def main():
    args = get_args()

    data = h5py.File(args.h5file, "r")
    (ensemble,) = get_ensemble(
        data, beta=args.beta, mF=args.mF, Nt=args.Nt, Ns=args.Ns
    )

    auto_smear, auto_point = ps_correlator_autocorrelation(ensemble, args)

    metadata = {
        "ensemble_name": args.ensemble_name,
        "beta": args.beta,
        "mAS": args.mAS,
        "Nt": args.Nt,
        "Ns": args.Ns,
    }
    dump_dict(
        {**metadata,
         "tau_exp_ps_correlator_smear": auto_smear,
         "tau_exp_ps_correlator_point": auto_point},
        args.output_file_mean,
    )


if __name__ == "__main__":
    main()
