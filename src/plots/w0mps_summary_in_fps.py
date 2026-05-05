#!/usr/bin/env python3

import matplotlib.pyplot as plt
from ..dump import read_sample_files
from ..plots_common import (
    save_or_show,
    ch_tag,
    channel_color,
    add_figure_legend_axes,
    TWO_COLUMN,
    MPS_left_CUT,
    get_fps_fit,
    plot_ps_ths,
    plot_axpb_y,
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


def plot(fit_pars):
    summary_fig, mass_axs = plt.subplots(
        1, 2, num="Summary", figsize=(TWO_COLUMN, 4), layout="constrained"
    )
    mass_ax = mass_axs[0]
    decay_ax = mass_axs[1]

    mass_ax.set_xlabel(r"$m_{\rm PS}^2 / f_{\rm PS}^2$")
    mass_ax.set_ylabel(r"$m_{\rm M}^2 / f_{\rm PS}^2$")
    
    decay_ax.set_ylabel(r"$f_{\rm M}^2 / f_{\rm PS}^2$")
    decay_ax.set_xlabel(r"$m_{\rm PS}^2 / f_{\rm PS}^2$")

    for parameter in fit_pars:
        if parameter["channel"] == "ps":
            fm_ps = parameter["F_samples"].samples
            lf_ps = parameter["Lf_samples"].samples
    
    
    fps_scale = get_fps_fit(fm_ps, lf_ps)
    fps_scale_mean = fps_scale.mean(axis=0)

    plot_ps_ths(mass_ax, 1, r"$m_{\rm PS}^2 / f_{\rm PS}^2$", 0.8, scale=fps_scale)
    plot_ps_ths(mass_ax, 2, r"$4m_{\rm PS}^2 / f_{\rm PS}^2$", 0.6, scale=fps_scale)
    plot_ps_ths(mass_ax, 3, r"$9m_{\rm PS}^2 / f_{\rm PS}^2$", 0.4, scale=fps_scale)
    plot_ps_ths(mass_ax, 4, r"$16m_{\rm PS}^2 / f_{\rm PS}^2$", 0.3, scale=fps_scale)

    mass_ax.set_ylim(0, 300)
    mass_ax.set_xlim(MPS_left_CUT/ fps_scale_mean[0], 0.4 / fps_scale_mean[-1])

    decay_ax.set_xlim(MPS_left_CUT/ fps_scale_mean[0], 0.4 / fps_scale_mean[-1])
    decay_ax.set_ylim(0, 6)
    
   
    for ch in ["v", "t", "av", "at", "s"]:
        
        for parameter in fit_pars:

            if "M_samples" not in parameter:
                    continue
            
            if parameter["channel"] == ch:

                plot_axpb_y(
                    mass_ax,
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
                    decay_ax,
                    parameter["F_samples"].samples,
                    parameter["Lf_samples"].samples,
                    r"$ \rm " + ch_tag(ch) + "$",
                    0.8,
                    channel_color(ch),
                    scale=fps_scale,
                )


    add_figure_legend_axes(summary_fig, mass_axs, 10, title=None)
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
