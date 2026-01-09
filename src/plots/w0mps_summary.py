#!/usr/bin/env python3

import matplotlib.pyplot as plt
from ..dump import read_sample_files
from ..plots_common import (
    save_or_show,
    ch_tag,
    channel_color,
    add_figure_legend,
    add_figure_legend_axes,
    ONE_COLUMN,
    TWO_COLUMN,
)
from argparse import ArgumentParser
import numpy as np


def get_args():
    parser = ArgumentParser()

    parser.add_argument(
        "data_filenames",
        nargs="+",
        metavar="sample_filename",
        help="Filenames of sample files containing data to plot",
    )
    parser.add_argument(
        "--fit_parameters",
        nargs="+",
        metavar="sample_filename",
        help="Filenames of sample files containing data to plot",
    )
    parser.add_argument(
        "--plot_file_summary",
        default=None,
        help=(
            "Where to place the resulting summary plot (fit results only). "
            "Default is to output to screen."
        ),
    )
    parser.add_argument(
        "--plot_styles",
        default="styles/paperdraft.mplstyle",
        help="Stylesheet to use for plots",
    )
    return parser.parse_args()


def plot_axpb_y(ax, A, L, ch, alpha, color, hatch=None):
    n_fit = 1000
    Yfit = np.zeros(shape=(A.shape[0], n_fit))

    x_min, x_max = ax.get_xlim()
    x_i = 0
    x_f = np.sqrt(x_max)
    x = np.linspace(x_i, x_f, n_fit)

    y_up = np.zeros(n_fit)
    y_dn = np.zeros(n_fit)

    for n in range(A.shape[0]):
        Yfit[n] = A[n] * (1 + L[n] * x**2)

    for i in range(n_fit):
        y_err = Yfit[0:-1, i].std()
        y_up[i] = Yfit[:, i].mean() + y_err
        y_dn[i] = Yfit[:, i].mean() - y_err

    ax.fill_between(
        x**2, y_up, y_dn, alpha=alpha, label=ch, facecolor=color, edgecolor=None, hatch=hatch
    )


def plot(fit_pars):
    summary_fig, summary_axs = plt.subplots(
        1, 2, num="Summary", figsize=(TWO_COLUMN, 4), layout="constrained"
    )
    summary_ax = summary_axs[0]
    summary_ax.plot([0, 1], [0, 1], "--", color="C0")
    summary_ax.set_ylim(0, 1.6)
    summary_ax.set_xlim(0, 0.4)
    summary_ax.set_xlabel(r"$\hat{m}_{\rm PS}^2$")
    summary_ax.set_ylabel(r"$\hat{m}_{\rm M}^2$")

    ax2 = summary_axs[1] #summary_ax.twinx() 
    ax2.set_ylabel(r"$\hat{f}_{\rm M}^2$")
    ax2.set_xlabel(r"$\hat{m}_{\rm PS}^2$")
    ax2.set_xlim(0, 0.4)
    ax2.set_ylim(0, 0.04)

   
    for ch in ["v", "t", "av", "at", "s"]:
        
        for parameter in fit_pars:

            if "M_samples" not in parameter:
                    continue
            
            if parameter["channel"] == ch:

                plot_axpb_y(
                    summary_ax,
                    parameter["M_samples"].samples,
                    parameter["Lm_samples"].samples,
                    r"$ \rm " + ch_tag(ch) + "$",
                    0.8,
                    channel_color(ch),
                )
    

    for ch in ["ps", "v", "av"]:
        
        for parameter in fit_pars:
            if parameter["channel"] == ch:

                plot_axpb_y(
                    ax2,
                    parameter["F_samples"].samples,
                    parameter["Lf_samples"].samples,
                    r"$ \rm " + ch_tag(ch) + "$",
                    0.8,
                    channel_color(ch),
                )


    add_figure_legend_axes(summary_fig, summary_axs, 7, title=None)
    #add_figure_legend(summary_fig, 3, title=None)

    

    return summary_fig


def main():
    args = get_args()
    plt.style.use(args.plot_styles)
    fit_pars = read_sample_files(args.fit_parameters, group_key="channel")
    summary_fig = plot(fit_pars)
    save_or_show(summary_fig, args.plot_file_summary)


if __name__ == "__main__":
    main()
