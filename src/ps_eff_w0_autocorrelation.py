#!/usr/bin/env python3

from argparse import ArgumentParser, FileType
import h5py
from uncertainties import ufloat
import numpy as np
import re

from .dump import dump_dict
from flow_analysis.stats.autocorrelation import exp_autocorrelation_fit
from flow_analysis.readers import readers
from flow_analysis.measurements.scales import bootstrap_ensemble_w0
from flow_analysis.stats.bootstrap import bootstrap_finalize
#from .mass import get_args
from .read_hdf5 import get_ensemble, filter_configurations
from .utils import get_index_separation
from .flow import read_flows


def get_args():
    parser = ArgumentParser()

    parser = ArgumentParser(
        description="Compute the mass and matrix element from correlators in an HDF5 file"
    )
    parser.add_argument("h5file", help="The file to read")
    parser.add_argument("--flow_filename", help="Filename of flow log to analyse")
    parser.add_argument(
        "--filetype",
        choices=list(readers),
        default="hdf5",
        help="How to interpret the input file",
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
        "--W0", type=float, help="Threshold value at which to take t = w0^2"
    )
    parser.add_argument(
        "--operator",
        choices=["plaq", "sym"],
        default="sym",
        help="Flow operator to use",
    )
    parser.add_argument(
        "--output_file_mean",
        type=FileType("w"),
        default="-",
        help="Where to output the mean and uncertainty of w0",
    )
    parser.add_argument(
        "--output_file_samples",
        type=FileType("w"),
        default=None,
        help="Where to output the bootstrap samples for w0",
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
        "--E0_plateau_start",
        type=int,
        default=None,
        help="Time slice at which plateau starts",
    )
    parser.add_argument(
        "--E0_plateau_end", type=int, default=None, help="Time slice at which plateau ends"
    )
    parser.add_argument(
        "--n_smear_max",
        type=int,
        default=None,
        help="max smearing level",
    )

    return parser.parse_args()


def ps_eff_w0_autocorrelation(ensemble, args):
    indices = np.asarray(
        [
            int(re.match(".*n([0-9]+)$", filename.decode()).groups()[0])
            for filename in ensemble["configurations"]
        ]
    )
    filtered_indices = filter_configurations(ensemble, args.min_trajectory, args.max_trajectory, args.trajectory_step)

    #print(ensemble["source_N0_sink_N0"].keys())
    
    corr_ps_smear = ensemble[f"source_N{args.n_smear_max}_sink_N{args.n_smear_max}/TRIPLET g5"][args.E0_plateau_start-1:args.E0_plateau_start+2, filtered_indices]
    corr_ps_point = ensemble["source_N0_sink_N0/TRIPLET g5"][args.E0_plateau_start-1:args.E0_plateau_start+2, filtered_indices]
    
    
    meff_smear = np.log(
        ( corr_ps_smear[0,:]  )
        / ( corr_ps_smear[1,:])
    )
    
    meff_point = np.log(
        ( corr_ps_point[0,:]  )
        / ( corr_ps_point[1,:])
    )
    
    #print("meff shape =", meff.shape)
    #print(meff)

    flows = read_flows(args)
    thinned_flows = flows.thin(
        min_trajectory=args.min_trajectory,
        max_trajectory=args.max_trajectory,
        trajectory_step=args.trajectory_step,
    )
    w0_samples = bootstrap_ensemble_w0(
            thinned_flows,
            args.W0,
            operator=args.operator,
        )
    
    if thinned_flows.trajectories.shape == indices[filtered_indices].shape:
        #print(f"Thinned flows trajectories: {thinned_flows.trajectories}")
        #print(f"Indices: {indices[filtered_indices]}")
        if (thinned_flows.trajectories == indices[filtered_indices]).all():
            w0_mean = bootstrap_finalize(w0_samples)
            flow_time_index = abs(thinned_flows.times - (w0_mean.nominal_value)**2).argmin()
            energy_density = {"sym": thinned_flows.Ecs, "plaq": thinned_flows.Eps}[args.operator]
        
            
            trajectory_separation = get_index_separation(indices[filtered_indices])
            ps_eff_w0_smear = exp_autocorrelation_fit(meff_smear * energy_density[:,flow_time_index]) * trajectory_separation
            ps_eff_w0_point = exp_autocorrelation_fit(meff_point * energy_density[:,flow_time_index]) * trajectory_separation

        else:
            ps_eff_w0_smear = ufloat(np.nan, np.nan)
            ps_eff_w0_point = ufloat(np.nan, np.nan)
    else:
        ps_eff_w0_smear = ufloat(np.nan, np.nan)
        ps_eff_w0_point = ufloat(np.nan, np.nan)

    return ps_eff_w0_smear, ps_eff_w0_point


def main():
    args = get_args()

    data = h5py.File(args.h5file, "r")
    (ensemble,) = get_ensemble(
        data, beta=args.beta, mF=args.mF, Nt=args.Nt, Ns=args.Ns
    )

    auto_smear, auto_point = ps_eff_w0_autocorrelation(ensemble, args)

    metadata = {
        "ensemble_name": args.ensemble_name,
        "beta": args.beta,
        "mF": args.mF,
        "Nt": args.Nt,
        "Ns": args.Ns,
    }
    dump_dict(
        {**metadata,
         "tau_exp_ps_eff_w0_point": auto_point,
         "tau_exp_ps_eff_w0_smear": auto_smear},
        args.output_file_mean,
    )


if __name__ == "__main__":
    main()
