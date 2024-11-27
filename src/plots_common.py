#!/usr/bin/env python3

from argparse import ArgumentParser, SUPPRESS

import matplotlib.pyplot as plt
import pandas as pd
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


def get_standard_plot_args(fit_results=False, external_data=False):
    parser = ArgumentParser()

    parser.add_argument(
        "data_filenames",
        nargs="+",
        metavar="sample_filename",
        help="Filenames of sample files containing data to plot",
    )
    parser.add_argument(
        "--fit_results",
        nargs="+",
        metavar="fit_result",
        default=[],
        help=("Filenames of fit result files to plot" if fit_results else SUPPRESS),
    )
    parser.add_argument(
        "--external_data",
        default=None,
        help=("Filename of any external data to plot" if external_data else SUPPRESS),
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


def standard_plot_main(plot_function, **args_options):
    args = get_standard_plot_args(**args_options)
    plt.style.use(args.plot_styles)
    data = read_sample_files(args.data_filenames)

    external_data = (
        pd.read_csv(args.external_data) if args.external_data is not None else None
    )

    fit_results = read_sample_files(args.fit_results, group_key="beta")

    save_or_show(
        plot_function(data, external_data=external_data, fit_results=fit_results),
        args.plot_file,
    )


def beta_color(b):
    return {
        6.6: "C0",
        6.65: "C1",
        6.7: "C2",
        6.75: "C3",
        6.8: "C4",
        6.9: "C5",
    }.get(b, b)


def channel_color(ch):
    return {
        "ps": "C0",
        "v": "C1",
        "t": "C2",
        "s": "C3",
        "av": "C4",
        "at": "C5",
        "rhoE1": "C6",
    }.get(ch, ch)


def ch_tag(ch):
    return {
        "rhoE1": r"v^\prime",
    }.get(ch, ch)
