from argparse import ArgumentParser, FileType

from flow_analysis.readers import readers
from flow_analysis.measurements.scales import bootstrap_ensemble_w0
from flow_analysis.stats.bootstrap import bootstrap_finalize

import numpy as np
from uncertainties import ufloat

from .dump import dump_dict, dump_samples


def get_args():
    parser = ArgumentParser()

    parser.add_argument("flow_filename", help="Filename of flow log to analyse")
    parser.add_argument(
        "--filetype",
        choices=list(readers),
        default="hirep",
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
    parser.add_argument("--ensemble_name", default=None, help="Name of ensemble")

    return parser.parse_args()


def main():
    args = get_args()

    if args.trajectory_step == 0:
        w0_samples = []
        w0_mean = ufloat(np.nan, np.nan)
        trajectories = 0
    else:
        flows = readers[args.filetype](args.flow_filename).thin(
            min_trajectory=args.min_trajectory,
            max_trajectory=args.max_trajectory,
            trajectory_step=args.trajectory_step,
        )
        w0_samples = bootstrap_ensemble_w0(flows, args.W0, operator=args.operator)
        w0_mean = bootstrap_finalize(w0_samples)
        trajectories = len(flows.trajectories)

    dump_dict(
        {
            "ensemble_name": args.ensemble_name,
            "w0": w0_mean,
            "trajectory_step": args.trajectory_step,
            "num_configs": trajectories,
        },
        args.output_file_mean,
    )
    if args.output_file_samples:
        dump_samples(
            {"ensemble_name": args.ensemble_name, "w0_samples": w0_samples},
            args.output_file_samples,
        )


if __name__ == "__main__":
    main()
