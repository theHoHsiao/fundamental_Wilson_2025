#!/usr/bin/env python3

from argparse import ArgumentParser
import re

from flow_analysis.readers import read_flows_hirep
import h5py
import numpy as np


def get_args():
    parser = ArgumentParser(
        description=(
            "Get gradient flow data from log files,"
            "and collate them into a single HDF5 file."
        )
    )

    parser.add_argument(
        "flow_filenames",
        nargs="+",
        metavar="flow_filename",
        help="Filename of gradient flow log file",
    )
    parser.add_argument(
        "--h5_filename",
        required=True,
        help="Where to place the combined HDF5 file.",
    )
    return parser.parse_args()


def parse_cfg_filename(filename):
    return re.match(
        r".*/(?P<runid>[^/]*)_(?P<NT>[0-9]+)x(?P<NX>[0-9]+)x(?P<NY>[0-9]+)x"
        r"(?P<NZ>[0-9]+)nc(?P<Nc>[0-9]+)(?:r(?P<rep>[A-Z]+))?(?:nf(?P<nf>[0-9]+))?"
        r"(?:b(?P<beta>[0-9]+\.[0-9]+))?(?:m(?P<mass>-?[0-9]+\.[0-9]+))?"
        r"n(?P<cfg_idx>[0-9]+)",
        filename,
    ).groupdict()


def get_filename_metadata(metadata, content):
    if content[0] == "[IO][0]Configuration" and content[2] == "read":
        provisional_metadata = parse_cfg_filename(content[1])
        keys = {
            "NT": int,
            "NX": int,
            "NY": int,
            "NZ": int,
            "Nc": int,
            "rep": str,
            "beta": float,
        }
        rep_suffix = {"FUN": "F", "ASY": "AS", "SYM": "S", "ADJ": "ADJ"}[
            provisional_metadata["rep"]
        ]
        # Hack to account for HiRep sometimes being inconsistent about minus signs
        # Works for this data as we never go to very heavy mass (positive m0)
        if not provisional_metadata["mass"].startswith("-"):
            provisional_metadata["mass"] = "-" + provisional_metadata["mass"]

        provisional_metadata[f"n{rep_suffix}"] = provisional_metadata["nf"]
        provisional_metadata[f"m{rep_suffix}"] = provisional_metadata["mass"]
        keys[f"n{rep_suffix}"] = int
        keys[f"m{rep_suffix}"] = float
        for key, dtype in keys.items():
            value = dtype(provisional_metadata[key])
            if key in metadata and metadata[key] != value:
                message = f"Metadata inconsistent: {metadata[key]} != {value}."
                raise ValueError(message)
            metadata[key] = value


def process_file(flow_filename, h5file):
    flows = read_flows_hirep(flow_filename, metadata_callback=get_filename_metadata)
    group_name = "gflow_{NT}x{NX}x{NY}x{NZ}b{beta}m{mAS}".format(**flows.metadata)
    group = h5file.create_group(group_name)
    group.create_dataset("beta", data=flows.metadata["beta"])
    group.create_dataset("configurations", data=flows.cfg_filenames.astype("S"))
    group.create_dataset("gauge group", data=f"SP({flows.metadata['Nc']})")
    group.create_dataset(
        "lattice",
        data=np.asarray([flows.metadata[key] for key in ["NT", "NX", "NY", "NZ"]]),
    )
    group.create_dataset("plaquette", data=flows.plaquettes)
    group.create_dataset("quarkmasses", data=[flows.metadata["mAS"]])
    group.create_dataset("flow type", data=flows.metadata.get("flow_type"))

    group.create_dataset("flow times", data=flows.times)
    group.create_dataset("topological charge", data=flows.Qs)
    group.create_dataset("energy density plaq", data=flows.Eps)
    group.create_dataset("energy density sym", data=flows.Ecs)


def main():
    args = get_args()
    with h5py.File(args.h5_filename, "w-") as h5file:
        for flow_filename in args.flow_filenames:
            process_file(flow_filename, h5file)


if __name__ == "__main__":
    main()
