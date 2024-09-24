#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

from flow_analysis.readers import readers
from flow_analysis.measurements.Q import Q_fit
from flow_analysis.stats.autocorrelation import exp_autocorrelation_fit


def get_args():
    parser = ArgumentParser()

    parser.add_argument("flow_filename", help="Filename of flow log to analyse")
    parser.add_argument("--filetype", choices=list(readers), default="hirep", help="How to interpret the input file")
    parser.add_argument("--min_trajectory", type=int, default=None, help="Lowest trajectory index to consider")
    parser.add_argument("--max_trajectory", type=int, default=None, help="Highest trajectory index to consider")
    parser.add_argument("--output_file", type=FileType("w"), default="-", help="Where to output the mean and uncertainty of w0")

    return parser.parse_args()


def main():
    args = get_args()

    flows = readers[args.filetype](args.flow_filename).thin(
        min_trajectory=args.min_trajectory,
        max_trajectory=args.max_trajectory,
    )
    Q0, sigma_Q = Q_fit(flows)
    tau_exp_Q = exp_autocorrelation_fit(flows.Q_history())

    print(f"{Q0.nominal_value},{Q0.std_dev},{sigma_Q.nominal_value},{sigma_Q.std_dev},{tau_exp_Q.nominal_value},{tau_exp_Q.std_dev}", file=args.output_file)


if __name__ == "__main__":
    main()
