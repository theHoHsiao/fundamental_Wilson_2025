#!/usr/bin/env python3

from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np

from ..dump import read_files
from ..plots_common import (
    beta_iterator,
    errorbar_ufloat,
    save_or_show,
    add_figure_legend,
    TWO_COLUMN,
)


def get_args():
    parser = ArgumentParser()

    parser.add_argument(
        "data_filenames",
        nargs="+",
        metavar="plot_data_filename",
        help="Filenames of PCAC mass results",
    )
    parser.add_argument(
        "--linear_fit_filenames",
        nargs="+",
        metavar="fit_data_filename",
        help="Filenames of PCAC mass results to fit linearly",
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


def plot(plot_data, linear_fit_data):
    fig, ax = plt.subplots(layout="constrained", figsize=(TWO_COLUMN, 3))

    ax.set_xlabel("$am_0$")
    ax.set_ylabel(r"$am_{\mathrm{PCAC}}$")

    xmin, xmax = get_xlim(plot_data)
    ymin, ymax = get_ylim(plot_data)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    x_range = np.linspace(xmin, xmax, 1000)

    for beta, colour, marker in beta_iterator(sorted(set(plot_data.beta))):
        subset = plot_data[plot_data.beta == beta]
        errorbar_ufloat(
            ax, subset.mAS, subset.mPCAC, color=colour, marker=marker, label=f"{beta}"
        )

        linear_fit_subset = linear_fit_data[linear_fit_data.beta == beta]
        if len(linear_fit_subset) > 0:
            linear_fit_result = np.polynomial.polynomial.Polynomial.fit(
                linear_fit_subset.mAS,
                [v.nominal_value for v in linear_fit_subset.mPCAC],
                1,
                w=[1 / v.std_dev for v in linear_fit_subset.mPCAC],
            )
            ax.plot(x_range, linear_fit_result(x_range), color=colour, dashes=(5, 3))

        if len(subset) >= 4:
            quadratic_fit_result = np.polynomial.polynomial.Polynomial.fit(
                subset.mAS,
                [v.nominal_value for v in subset.mPCAC],
                2,
                w=[1 / v.std_dev for v in subset.mPCAC],
            )
            ax.plot(x_range, quadratic_fit_result(x_range), color=colour, dashes=(2, 3))

    add_figure_legend(fig)
    return fig


def main():
    args = get_args()
    plt.style.use(args.plot_styles)
    plot_data = read_files(args.data_filenames)
    fit_data = read_files(args.linear_fit_filenames)
    save_or_show(plot(plot_data, fit_data), args.plot_filename)


if __name__ == "__main__":
    main()
