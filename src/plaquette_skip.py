#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

from flow_analysis.stats.autocorrelation import integrated_autocorrelation_time

import h5py
import logging
import re

from .bootstrap import sample_bootstrap_0d, bootstrap_finalize, get_rng
from .dump import dump_dict, dump_samples
from .read_hdf5 import get_ensemble, filter_configurations
from .utils import get_index_separation


def get_args():
    parser = ArgumentParser(
        description=(
            "Get ensemble plaquettes from HDF5, "
            "compute mean and autocorrelation time, and write to a file"
        )
    )

    parser.add_argument("h5file", help="The file to read")
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
        "--mF",
        type=float,
        default=None,
        help="The antisymmetric fermion mass of the ensemble to analyse",
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
        "--trajectory_step_auto",
        type=int,
        default=1,
        help="Interval of trajectories to consider for autocorrelation time",
    )
    parser.add_argument(
        "--bin_size",
        type=int,
        default=1,
        help="Number of consecutive configurations to bin together",
    )
    parser.add_argument(
        "--output_file_mean",
        type=FileType("w"),
        default="-",
        help="Where to output the mean values. (Defaults to stdout.)",
    )
    parser.add_argument(
        "--output_file_samples",
        type=FileType("w"),
        default="-",
        help="Where to output the bootstrap samples.",
    )
    return parser.parse_args()


def get_volumes(ensemble):
    nt, nx, ny, nz = ensemble["lattice"]
    if nx != ny or ny != nz:
        raise NotImplementedError("This code expects NX == NY == NZ")
    return {"Nt": nt, "Ns": nx}


def get_mass(ensemble):
    if len(ensemble["quarkmasses"]) != 1:
        raise NotImplementedError("This code expects one fermion mass per ensemble.")
    return ensemble["quarkmasses"][0]


def get_cfg_indices(cfgs, start_cfg=0, end_cfg=None, cfg_step=1):
    trajectory_indices = []
    array_indices = []
    for array_index, filename in enumerate(cfgs):
        trajectory_index = int(re.match(".*n([0-9]+)$", filename.decode()).groups()[0])
        if trajectory_index < start_cfg:
            continue
        if end_cfg is not None and trajectory_index > end_cfg:
            continue
        if (trajectory_index - start_cfg) % cfg_step == 0:
            trajectory_indices.append(trajectory_index)
            array_indices.append(array_index)

    separation = get_index_separation(trajectory_indices)
    if separation != cfg_step and cfg_step > 1:
        message = f"Requested a delta_traj of {cfg_step}, but actually see {separation}"
        logging.warning(message)

    return trajectory_indices, array_indices


def avg_plaquette(ensemble, args):
    result = {}
    result.update(get_volumes(ensemble))
    result["ensemble_name"] = args.ensemble_name
    result["mF"] = get_mass(ensemble)
    result["beta"] = ensemble["beta"][()]

    trajectory_indices, _ = get_cfg_indices(
        ensemble["configurations"], args.min_trajectory, args.max_trajectory, args.trajectory_step_auto
    )
    result["Ncfg"] = len(trajectory_indices)
    result["delta_traj"] = args.trajectory_step_auto


    spectrum_trajectory_indices, _ = get_cfg_indices(
        ensemble["configurations"], args.min_trajectory, args.max_trajectory, args.trajectory_step
    )
    result["Ncfg_spectrum"] = len(spectrum_trajectory_indices)
    result["delta_traj_spec"] = args.trajectory_step
    filtered_indices = filter_configurations(ensemble, args.min_trajectory, args.max_trajectory, args.trajectory_step)


    plaq_auto = ensemble["plaquette"]
    plaq_spec = ensemble["plaquette"][filtered_indices]
    result["plaquette"] = sample_bootstrap_0d(plaq_spec, get_rng(ensemble.name), bin_size=args.bin_size)
    

    result["avg_plaquette"] = bootstrap_finalize(result["plaquette"])
    result["tau_exp_plaq"] = (
        integrated_autocorrelation_time(plaq_auto ) * result["delta_traj"]
    )
    return result


def main():
    args = get_args()

    data = h5py.File(args.h5file, "r")
    ensemble, = get_ensemble(
        data,
        beta=args.beta,
        mF=args.mF,
        Nt=args.Nt,
        Ns=args.Ns,
    )
    #results = sample_bootstrap_0d(plaq, get_rng(ensemble.name))
    result = avg_plaquette(
        ensemble,
        args,
    )

    metadata_fields = [
        "Nt",
        "Ns",
        "mF",
        "beta",
        "Ncfg",
        "delta_traj",
        "Ncfg_spectrum",
        "delta_traj_spec",
        "ensemble_name",
    ]
    dump_dict(
        {k: result[k] for k in [*metadata_fields, "avg_plaquette", "tau_exp_plaq"]},
        args.output_file_mean,
    )
    if args.output_file_samples:
        dump_samples(
            {k: result[k] for k in [*metadata_fields, "plaquette"]},
            args.output_file_samples,
        )


if __name__ == "__main__":
    main()
