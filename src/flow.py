from argparse import ArgumentParser, FileType

from flow_analysis.readers import readers
from flow_analysis.measurements.scales import bootstrap_ensemble_w0
from flow_analysis.stats.bootstrap import bootstrap_finalize
from flow_analysis.stats.autocorrelation import exp_autocorrelation_fit

import h5py
import numpy as np
from uncertainties import ufloat

from .dump import dump_dict, dump_samples
from .read_hdf5 import get_ensemble
from .utils import get_index_separation


def get_args():
    parser = ArgumentParser()

    parser.add_argument("flow_filename", help="Filename of flow log to analyse")
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
        "W0", type=float, help="Threshold value at which to take t = w0^2"
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

    return parser.parse_args()


def fit_w0_tau_exp(w0, flows, operator="sym"):
    flow_time_index = abs(flows.times - w0**2).argmin()
    energy_density = {"sym": flows.Ecs, "plaq": flows.Eps}[operator]
    raw_tau_exp = exp_autocorrelation_fit(energy_density[:, flow_time_index])
    return raw_tau_exp * get_index_separation(flows.trajectories)


def read_flows(args):
    if args.filetype != "hdf5":
        return readers[args.filetype](args.flow_filename)

    with h5py.File(args.flow_filename, "r") as h5file:
        (ensemble,) = get_ensemble(h5file, args.beta, args.mF, args.Nt, args.Ns)
        ensemble_name = ensemble.name
    return readers["hdf5"](args.flow_filename, group_name=ensemble_name)


def main():
    args = get_args()

    if args.trajectory_step == 0:
        w0_samples = []
        w0_mean = ufloat(np.nan, np.nan)
        trajectories = 0
        trajectory_step = args.trajectory_step
        tau_exp_w0 = ufloat(np.nan, np.nan)
    else:
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
        w0_mean = bootstrap_finalize(w0_samples)
        tau_exp_w0 = fit_w0_tau_exp(
            w0_mean.nominal_value,
            flows,
            operator=args.operator,
        )
        trajectory_step = get_index_separation(thinned_flows.trajectories)
        trajectories = len(thinned_flows.trajectories)

    dump_dict(
        {
            "ensemble_name": args.ensemble_name,
            "w0": w0_mean,
            "tau_exp_w0": tau_exp_w0,
            "delta_traj_w0": trajectory_step,
            "Ncfg_GF": trajectories,
        },
        args.output_file_mean,
    )
    if args.output_file_samples:
        dump_samples(
            {
                "ensemble_name": args.ensemble_name,
                "w0_value": w0_mean.nominal_value,
                "w0_samples": w0_samples,
            },
            args.output_file_samples,
        )


if __name__ == "__main__":
    main()
