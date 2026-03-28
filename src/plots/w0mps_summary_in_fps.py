#!/usr/bin/env python3

import matplotlib.pyplot as plt
from ..dump import read_sample_files
from ..plots_common import (
    save_or_show,
    ch_tag,
    channel_color,
    add_figure_legend_axes,
    TWO_COLUMN,
    MPS_CUT,
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


def plot_axpb_y(ax, A, L, ch, alpha, color, hatch=None, scale=None):
    n_fit = 1000
    Yfit = np.zeros(shape=(A.shape[0], n_fit))

    x_i = np.sqrt(MPS_CUT)
    x_f = np.sqrt(0.4)
    x = np.linspace(x_i, x_f, n_fit)

    y_up = np.zeros(n_fit)
    y_dn = np.zeros(n_fit)

    for n in range(A.shape[0]):
        Yfit[n] = A[n] * (1 + L[n] * x**2)
    
    if scale is not None:
        Yfit = Yfit / scale.mean(axis=0)
    
    #print(Yfit[0,:])

    for i in range(n_fit):
        y_err = Yfit[0:-1, i].std()
        y_up[i] = Yfit[:, i].mean() + y_err
        y_dn[i] = Yfit[:, i].mean() - y_err

    ax.fill_between(
        x**2, y_up, y_dn, alpha=alpha, label=ch, facecolor=color, edgecolor=None, hatch=hatch
    )


def get_fps_fit(A, L,):
    n_fit = 1000
    Yfit = np.zeros(shape=(A.shape[0], n_fit))

    x_i = np.sqrt(MPS_CUT)
    x_f = np.sqrt(0.4)
    x = np.linspace(x_i, x_f, n_fit)

    for n in range(A.shape[0]):
        Yfit[n] = A[n] * (1 + L[n] * x**2)
    

    return Yfit


def plot_ps_ths(ax, times_of_mps, label, alpha, scale=None):
    n_fit=1000
    x_i = np.sqrt(MPS_CUT)
    x_f = np.sqrt(0.4)
    m_ps = np.linspace(x_i, x_f, n_fit)

    Y = (times_of_mps * m_ps)**2
    if scale is not None:
        Y = Y / scale.mean(axis=0)

    ax.plot(m_ps**2, Y, "--", color="k", linewidth=1, label=label, alpha=alpha)






def plot(fit_pars):
    summary_fig, summary_axs = plt.subplots(
        1, 2, num="Summary", figsize=(TWO_COLUMN, 4), layout="constrained"
    )
    summary_ax = summary_axs[0]
   

    summary_ax.set_ylim(0, 300)
    summary_ax.set_xlim(MPS_CUT, 0.4)
    summary_ax.set_xlabel(r"$\hat{m}_{\rm PS}^2$")
    summary_ax.set_ylabel(r"$\hat{m}_{\rm M}^2 / \hat{f}_{\rm PS}^2$")

    ax2 = summary_axs[1] #summary_ax.twinx() 
    ax2.set_ylabel(r"$\hat{f}_{\rm M}^2 / \hat{f}_{\rm PS}^2$")
    ax2.set_xlabel(r"$\hat{m}_{\rm PS}^2$")
    ax2.set_xlim(MPS_CUT, 0.4)
    ax2.set_ylim(0, 4)

    for parameter in fit_pars:
        if parameter["channel"] == "ps":
            fm_ps = parameter["F_samples"].samples
            lf_ps = parameter["Lf_samples"].samples
    
    
    fps_scale = get_fps_fit(fm_ps, lf_ps)
    plot_ps_ths(summary_ax, 1, r"$\hat{m}_{\rm PS}$", 0.8, scale=fps_scale)
    plot_ps_ths(summary_ax, 2, r"$2\hat{m}_{\rm PS}$", 0.6, scale=fps_scale)
    plot_ps_ths(summary_ax, 3, r"$3\hat{m}_{\rm PS}$", 0.4, scale=fps_scale)
    plot_ps_ths(summary_ax, 4, r"$4\hat{m}_{\rm PS}$", 0.3, scale=fps_scale)

    #summary_ax.plot([0, 1], [0, 1], "--", color="k", linewidth=1, label=r"$\hat{m}_{\rm PS}$", alpha=0.8)
    #summary_ax.plot([0, 1], [0, 4], "--", color="k", linewidth=1, label=r"$ 2\hat{m}_{\rm PS}$", alpha=0.6)
    #summary_ax.plot([0, 1], [0, 9], "--", color="k", linewidth=1, label=r"$ 3\hat{m}_{\rm PS}$", alpha=0.4)
    
   
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
                    scale=fps_scale,
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
                    scale=fps_scale,
                )


    add_figure_legend_axes(summary_fig, summary_axs, 10, title=None)
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
