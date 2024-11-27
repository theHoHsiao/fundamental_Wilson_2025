#!/usr/bin/env python3

from argparse import ArgumentParser
from collections import defaultdict
from os.path import basename

import h5py
import numpy as np

from .provenance import get_flat_metadata


def get_args():
    parser = ArgumentParser(
        description=(
            "Get ensemble plaquettes and other data from HMC logs, "
            "and package as HDF5."
        )
    )
    parser.add_argument(
        "hmc_filenames",
        nargs="+",
        metavar="hmc_filename",
        help="Filename of HMC log file",
    )
    parser.add_argument(
        "--h5_filename",
        required=True,
        help="Where to place the output HDF5 file.",
    )
    return parser.parse_args()


def is_first_non_none(plaquettes):
    return set(plaquettes[1:]) == {None}


def update(result, data):
    if "polyakov loop" in data and "polyakov loop" not in result:
        if is_first_non_none(result["plaquette"]):
            result["polyakov loop"] = [[None, None, None, None]] * len(
                result["plaquette"]
            )
        else:
            raise NotImplementedError("Inconsistent presence of Polyakov loop")

    if "polyakov loop" in result and "polyakov loop" not in data:
        data["polyakov loop"] = [None, None, None, None]

    for key, collection in result.items():
        if not isinstance(collection, list):
            continue
        if key == "lattice":
            continue
        result[key].append(data[key])


def check_global_metadata(result, metadata):
    if "nt" in result:
        for key, value in metadata.items():
            if result[key] != value:
                raise ValueError(f"Inconsistent value for {key}")
    else:
        result.update(metadata)


def read_hmc(filename):
    result = {
        "original_filename": filename,
        "trajectory": [],
        "plaquette": [],
        "tlen": [],
        "nsteps": [],
        "gentime": [],
        "accept": [],
        "deltaS": [],
        "configurations": [],
    }
    current_metadata = {
        "nAS": 0,
        "nF": 0,
        "nADJ": 0,
        "nS": 0,
    }
    trajectory_data = {}
    monomials_read = set()
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("[SYSTEM][0]MACROS=") or line.startswith(
                "[MAIN][0]Compiled with macros"
            ):
                if "-DREPR_ANTISYMMETRIC" in line:
                    rep = "AS"
                if "-DREPR_SYMMETRIC" in line:
                    rep = "S"
                if "-DREPR_FUNDAMENTAL" in line:
                    rep = "F"
                if "-DREPR_ADJOINT" in line:
                    rep = "ADJ"
                current_metadata["rep"] = rep

            if line.startswith("[FLOW][0]Starting a new run from a"):
                if "start" in result:
                    raise ValueError(f"Multiple starts in file {filename}")
                result["start"] = line.split()[6]
            if line.startswith("[GEOMETRY_INIT][0]Global size is") or line.startswith(
                "[GEOMETRY][0]Global size is"
            ):
                current_metadata["lattice"] = list(
                    map(int, line.split()[-1].split("x"))
                )
            if line.startswith("[MAIN][0]Gauge group:") or line.startswith(
                "[SYSTEM][0]Gauge group:"
            ):
                current_metadata["gauge group"] = line.split()[-1]
            if line.startswith("[IO][0]Configuration"):
                split_line = line.split()
                trajectory_data["configurations"] = basename(split_line[1].strip("[]"))
            if line.startswith("[ACTION][10]Monomial"):
                monomial_id = line.split()[1].strip(":")
                if monomial_id in monomials_read:
                    continue

                monomials_read.add(monomial_id)
                if "type = gauge," in line:
                    current_metadata["beta"] = float(line.split()[-1])
                elif "type = rhmc" in line or "type = hmc" in line:
                    if "type = hmc," in line:
                        current_metadata[f"n{rep}"] += 2
                    if "type = rhmc," in line:
                        current_metadata[f"n{rep}"] += 1

                    mass = float(line.split()[10].strip(","))
                    if result.get(f"m{rep}", mass) != mass:
                        raise NotImplementedError(
                            "Non-degenerate masses not currently supported"
                        )
                    current_metadata[f"m{rep}"] = mass
                else:
                    raise NotImplementedError(
                        f"Monomial not recognised in file {filename}"
                    )
            if line.startswith("[MAIN][0]Trajectory") and line.endswith("...\n"):
                if len(trajectory_data) > 2:
                    check_global_metadata(result, current_metadata)
                    update(result, trajectory_data)
                trajectory_data = defaultdict(lambda: None)
                trajectory_data["trajectory"] = int(line.split()[1].strip("#.:"))

            if line.startswith("[MD_INT][10]MD parameters: level=0"):
                trajectory_data["tlen"] = float(line.split()[3].split("=")[1])
                trajectory_data["nsteps"] = int(line.split()[4].split("=")[1])

            if line.startswith("[HMC][10]Configuration"):
                if line.endswith("accepted.\n"):
                    trajectory_data["accept"] = True
                elif line.endswith("rejected.\n"):
                    trajectory_data["accept"] = False
                else:
                    raise ValueError("Unexpected acceptance result")
            if line.startswith("[HMC][10][DeltaS ="):
                trajectory_data["deltaS"] = float(line.split()[2].split("]")[0])
            if (
                line.startswith("[MAIN][0]Initial plaquette")
                and not result["plaquette"]
            ):
                update(
                    result,
                    {
                        "trajectory": 0,
                        "plaquette": float(line.split()[-1]),
                        "gentime": None,
                        "accept": None,
                        "deltaS": None,
                        "tlen": None,
                        "nsteps": 0,
                        "configurations": trajectory_data.get("configurations"),
                    },
                )
            if line.startswith("[MAIN][0]Trajectory #") and (
                split_line := line.split()
            )[2:4] == ["generated", "in"]:
                trajectory_data["gentime"] = int(split_line[4].strip("[")) + 1e-6 * int(
                    split_line[6]
                )
            if line.startswith("[FUND_POLYAKOV][0]Polyakov direction"):
                split_line = line.split()
                if split_line[2] == "0":
                    trajectory_data["polyakov loop"] = []
                trajectory_data["polyakov loop"].append(
                    float(split_line[4]) + float(split_line[5]) * 1.0j
                )
            if line.startswith("[MAIN][0]Plaquette:"):
                trajectory_data["plaquette"] = float(line.split()[-1])
            if line.startswith("[HMC][0]Memory deallocated."):
                check_global_metadata(result, current_metadata)
                update(result, trajectory_data)
                trajectory_data = {}

    del result["rep"]
    return result


def name_ensemble(datum):
    return "{}x{}x{}x{}b{}m{}_{}".format(
        *datum["lattice"], datum["beta"], datum["mAS"], datum["start"]
    )


def _get_correct_type(data):
    typename = None
    if isinstance(data, str):
        return "S"
    for datum in data:
        if datum is None:
            continue
        if hasattr(datum, "__iter__"):
            typename = _get_correct_type(datum)
            if typename is not None:
                return typename
        else:
            return type(datum)
    return None


def get_correct_type(data):
    result = _get_correct_type(data)
    if result is None:
        raise ValueError("No types available in data.")
    return result


def write_hdf5(data, filename):
    with h5py.File(filename, "w") as h5file:
        for datum in data:
            group = h5file.create_group(f"hmc_{name_ensemble(datum)}")
            for key, value in datum.items():
                if isinstance(value, list) and not isinstance(value[0], int):
                    value = np.asarray(value, dtype=get_correct_type(value))
                group.create_dataset(key, data=value)

        provenance_group = h5file.create_group("_comment")
        for key, value in sorted(get_flat_metadata().items()):
            provenance_group.create_dataset(key, data=value)


def main():
    args = get_args()
    data = [read_hmc(filename) for filename in args.hmc_filenames]
    write_hdf5(data, args.h5_filename)


if __name__ == "__main__":
    main()
