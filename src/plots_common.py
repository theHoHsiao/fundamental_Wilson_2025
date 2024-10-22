#!/usr/bin/env python3

from argparse import ArgumentParser

import matplotlib.pyplot as plt
from uncertainties import UFloat

from .dump import read_sample_files


def save_or_show(fig, filename=None):
    if filename == "/dev/null":
        plt.close(fig)
    elif filename is not None:
        fig.savefig(filename, transparent=True)
        plt.close(fig)
    else:
        plt.show()


def is_ufloat_sequence(seq):
    if hasattr(seq, "values"):
        return isinstance(seq.values[0], UFloat)
    return isinstance(seq[0], UFloat)


def errorbar_ufloat(ax, x, y, *args, **kwargs):
    if is_ufloat_sequence(x):
        x_values = [xi.nominal_value for xi in x]
        x_errors = [xi.std_dev for xi in x]
    else:
        x_values = x
        x_errors = None

    if is_ufloat_sequence(y):
        y_values = [yi.nominal_value for yi in y]
        y_errors = [yi.std_dev for yi in y]
    else:
        y_values = y
        y_errors = None

    ax.errorbar(
        x_values,
        y_values,
        xerr=x_errors,
        yerr=y_errors,
        ls="none",
        *args,
        **kwargs,
    )


def get_standard_plot_args():
    parser = ArgumentParser()

    parser.add_argument(
        "data_filenames",
        nargs="+",
        metavar="sample_filename",
        help="Filenames of sample files containing data to plot",
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


def standard_plot_main(plot_function):
    args = get_standard_plot_args()
    plt.style.use(args.plot_styles)
    data = read_sample_files(args.data_filenames)
    save_or_show(plot_function(data), args.plot_file)


def beta_color(b):
    return {
        6.6: "k",
        6.65: "r",
        6.7: "b",
        6.75: "m",
        6.8: "g",
        6.9: "tab:brown",
    }.get(b, b)


def ch_tag(ch):
    return {
        "rhoE1": r"v^\prime",
    }.get(ch, ch)
