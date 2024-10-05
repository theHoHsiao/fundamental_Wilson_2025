#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

from flow_analysis.readers import readers
from flow_analysis.measurements.Q import Q_fit
from flow_analysis.stats.autocorrelation import exp_autocorrelation_fit

from .dump import dump_dict


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
        "--output_file",
        type=FileType("w"),
        default="-",
        help="Where to output the results",
    )
    parser.add_argument("--ensemble_name", default=None, help="Name of ensemble")

    return parser.parse_args()


def main():
    args = get_args()

    flows = readers[args.filetype](args.flow_filename).thin(
        min_trajectory=args.min_trajectory,
        max_trajectory=args.max_trajectory,
    )
    Q0, sigma_Q = Q_fit(flows)
    tau_exp_Q = exp_autocorrelation_fit(flows.Q_history())

    dump_dict(
        {
            "ensemble_name": args.ensemble_name,
            "Q0": Q0,
            "sigma_Q": sigma_Q,
            "tau_exp_Q": tau_exp_Q,
        },
        args.output_file,
    )


if __name__ == "__main__":
    main()
