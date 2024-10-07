#!/usr/bin/env python3

from argparse import ArgumentParser

from flow_analysis.readers import readers
from flow_analysis.measurements.scales import compute_wt_t

import matplotlib.pyplot as plt

from .plots_common import save_or_show


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
        "--plot_file",
        default=None,
        help="Where to place the resulting plot. Default is to output to screen.",
    )
    parser.add_argument(
        "--plot_styles",
        default="styles/paperdraft.mplstyle",
        help="Stylesheet to use for plots",
    )

    return parser.parse_args()


def plot(flows, threshold):
    fig, ax = plt.subplots(figsize=(3.5, 2.5), layout="constrained")

    ax.set_xlabel("$t / a^2$")
    ax.set_ylabel(r"$t \frac{\mathrm{d}}{\mathrm{d}t}\langle t^2 E\rangle$")

    for operator_idx, (operator, label, marker) in enumerate(
        [("plaq", "Plaquette", "o"), ("sym", "Clover", "s")]
    ):
        w_value, w_error = compute_wt_t(flows, operator)
        ax.errorbar(
            flows.times[1:-1],
            w_value,
            yerr=w_error,
            color=f"C{operator_idx}",
            ls="none",
            marker=marker,
            label=label,
            markersize=2,
            capsize=1,
        )

    ax.axhline(threshold, dashes=(2, 3), color="black")
    ax.legend(loc="best")
    return fig


def main():
    args = get_args()
    plt.style.use(args.plot_styles)

    flows = readers[args.filetype](args.flow_filename).thin(
        min_trajectory=args.min_trajectory,
        max_trajectory=args.max_trajectory,
        trajectory_step=args.trajectory_step,
    )
    save_or_show(plot(flows, args.W0), args.plot_file)


if __name__ == "__main__":
    main()
