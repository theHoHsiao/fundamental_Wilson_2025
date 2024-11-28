#!/usr/bin/env python3


import h5py
import numpy as np
import re

from .dump import dump_dict
from flow_analysis.stats.autocorrelation import exp_autocorrelation_fit
from .mass import get_args
from .read_hdf5 import get_ensemble
from .utils import get_index_separation


def ps_correlator_autocorrelation(ensemble, args):
    indices = np.asarray(
        [
            int(re.match(".*n([0-9]+)$", filename.decode()).groups()[0])
            for filename in ensemble["configurations"]
        ]
    )
    filtered_indices = (
        (indices >= args.min_trajectory) if args.min_trajectory is not None else True
    ) & ((indices <= args.max_trajectory) if args.max_trajectory is not None else True)

    corr_ps = ensemble["TRIPLET/g5"][args.plateau_start, filtered_indices]
    trajectory_separation = get_index_separation(indices)

    tau_ps_correlator = exp_autocorrelation_fit(corr_ps) * trajectory_separation

    return tau_ps_correlator


def main():
    args = get_args()

    data = h5py.File(args.h5file, "r")
    (ensemble,) = get_ensemble(
        data, beta=args.beta, mAS=args.mAS, Nt=args.Nt, Ns=args.Ns
    )

    auto = ps_correlator_autocorrelation(ensemble, args)

    metadata = {
        "ensemble_name": args.ensemble_name,
        "beta": args.beta,
        "mAS": args.mAS,
        "Nt": args.Nt,
        "Ns": args.Ns,
    }
    dump_dict(
        {**metadata, "tau_exp_ps_correlator": auto},
        args.output_file_mean,
    )


if __name__ == "__main__":
    main()
