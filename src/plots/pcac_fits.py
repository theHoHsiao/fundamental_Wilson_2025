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
        metavar="plot_data_filename",
        help="Filenames of PCAC mass results",
    )
    parser.add_argument(
        "--fit_filenames",
        nargs="+",
        metavar="fit_data_filename",
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


def plot(plot_data, fit_data):
    fig, ax = plt.subplots(layout="constrained", figsize=(7, 3))

    ax.set_xlabel("$am_0$")
    ax.set_ylabel(r"$am_{\mathrm{PCAC}}$")

    xmin, xmax = get_xlim(plot_data)
    ymin, ymax = get_ylim(plot_data)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    x_range = np.linspace(xmin, xmax, 1000)

    for beta_idx, beta in enumerate(sorted(set(plot_data.beta))):
        subset = plot_data[plot_data.beta == beta]
        colour = f"C{beta_idx}"
        errorbar_ufloat(ax, subset.mAS, subset.mPCAC, color=colour, marker=".")

        fit_subset = fit_data[fit_data.beta == beta]
        if len(fit_subset) == 0:
            continue
        for order, dashes in (1, (5, 3)), (2, (2, 3)):
            fit_result = np.polynomial.polynomial.Polynomial.fit(
                fit_subset.mAS,
                [v.nominal_value for v in fit_subset.mPCAC],
                order,
                w=[1 / v.std_dev for v in fit_subset.mPCAC],
            )
            ax.plot(x_range, fit_result(x_range), color=colour, dashes=dashes)

    return fig


def main():
    args = get_args()
    plt.style.use(args.plot_styles)
    plot_data = read_files(args.data_filenames)
    fit_data = read_files(args.fit_filenames)
    save_or_show(plot(plot_data, fit_data), args.plot_filename)


if __name__ == "__main__":
    main()
