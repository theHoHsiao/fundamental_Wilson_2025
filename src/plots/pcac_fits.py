#!/usr/bin/env python3

from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np

from ..dump import read_files
from ..plots_common import errorbar_ufloat, save_or_show


def get_args():
    parser = ArgumentParser()

    parser.add_argument(
        "data_filenames",
        nargs="+",
        metavar="data_filename",
        help="Filenames of PCAC mass results",
    )
    parser.add_argument(
        "--plot_filename",
        help="Where to place the generated plot. Default is to display on screen.",
    )
    parser.add_argument(
        "--plot_styles",
        default="styles/paperdraft.mplstyle",
        help="Stylesheet to use for plots",
    )
    return parser.parse_args()


def get_xlim(data):
    data_min = data.mAS.min()
    data_max = data.mAS.max()
    data_range = data_max - data_min
    return data_min - data_range / 10, data_max + data_range / 10


def get_ylim(data):
    return 0, 1.1 * max([value.nominal_value for value in data.mPCAC])


def plot(data):
    fig, ax = plt.subplots(layout="constrained", figsize=(7, 3))

    ax.set_xlabel("$am_0$")
    ax.set_ylabel(r"$am_{\mathrm{PCAC}}$")

    xmin, xmax = get_xlim(data)
    ymin, ymax = get_ylim(data)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    x_range = np.linspace(xmin, xmax, 1000)

    for beta_idx, beta in enumerate(sorted(set(data.beta))):
        subset = data[data.beta == beta]
        colour = f"C{beta_idx}"
        errorbar_ufloat(ax, subset.mAS, subset.mPCAC, color=colour, marker=".")
        if len(subset) < 4:
            continue

        for order, dashes in (1, (5, 3)), (2, (2, 3)):
            fit_result = np.polynomial.polynomial.Polynomial.fit(
                subset.mAS,
                [v.nominal_value for v in subset.mPCAC],
                order,
                w=[1 / v.std_dev for v in subset.mPCAC],
            )
            ax.plot(x_range, fit_result(x_range), color=colour, dashes=dashes)

    return fig


def main():
    args = get_args()
    plt.style.use(args.plot_styles)
    data = read_files(args.data_filenames)
    save_or_show(plot(data), args.plot_filename)


if __name__ == "__main__":
    main()
