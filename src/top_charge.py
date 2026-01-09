#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

from flow_analysis.fit_forms import gaussian
from flow_analysis.readers import readers
from flow_analysis.measurements.Q import Q_fit, flat_bin_Qs
from flow_analysis.stats.autocorrelation import integrated_autocorrelation_time

import numpy as np

from .dump import dump_dict
from .flow import read_flows
from .plots_common import ONE_COLUMN
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
        "--output_file",
        type=FileType("w"),
        default="-",
        help="Where to output the results",
    )
    parser.add_argument(
        "--plot_file",
        default=None,
        help=(
            "Where to place a plot of the history and histogram. "
            "If omitted, no plot is generated."
        ),
    )
    parser.add_argument(
        "--plot_styles",
        default="styles/paperdraft.mplstyle",
        help="Stylesheet to use for plots",
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
        "--mF",
        type=float,
        default=None,
        help="The fermion mass of the ensemble to analyse",
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


def compute_stats(flows):
    A, Q0, sigma_Q = Q_fit(flows, with_amplitude=True)
    tau_exp_Q = integrated_autocorrelation_time(flows.Q_history()) * get_index_separation(
        flows.trajectories
    )
    return {
        "Q_amplitude": A,
        "Q0": Q0,
        "sigma_Q": sigma_Q,
        "tau_exp_Q": tau_exp_Q,
    }


def plot(flows, results, plot_filename, plot_styles):
    import matplotlib.pyplot as plt

    plt.style.use(plot_styles)
    fig, (history_ax, histogram_ax) = plt.subplots(
        1,
        2,
        sharey=True,
        layout="constrained",
        gridspec_kw={"width_ratios": [3, 1]},
        figsize=(ONE_COLUMN, 2.5),
    )

    history_ax.step(flows.trajectories, flows.Q_history())

    Q_range, Q_counts = flat_bin_Qs(flows.Q_history())
    histogram_ax.step(Q_counts, Q_range - 0.5, label="Histogram")

    range_min = min(Q_range) - 0.5
    range_max = max(Q_range) + 0.5
    histogram_ax.set_ylim(range_min, range_max)

    smooth_Q_range = np.linspace(range_min, range_max, 1000)
    histogram_ax.plot(
        gaussian(
            smooth_Q_range,
            results["Q_amplitude"].nominal_value,
            results["Q0"].nominal_value,
            results["sigma_Q"].nominal_value,
        ),
        smooth_Q_range,
        label="Fit",
    )

    history_ax.set_ylabel("$Q$")
    history_ax.set_xlabel("Trajectory")
    histogram_ax.set_xlabel("Count")

    fig.savefig(plot_filename)


def main():
    args = get_args()

    flows = read_flows(args).thin(
        min_trajectory=args.min_trajectory,
        max_trajectory=args.max_trajectory,
        trajectory_step=args.trajectory_step,
    )
    results = compute_stats(flows)
    dump_dict(
        {
            "ensemble_name": args.ensemble_name,
            **results,
        },
        args.output_file,
    )
    if args.plot_file:
        plot(flows, results, args.plot_file, args.plot_styles)


if __name__ == "__main__":
    main()
