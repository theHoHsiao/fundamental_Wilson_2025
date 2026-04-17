#!/usr/bin/env python3

import matplotlib.pyplot as plt
from ..dump import read_sample_files
from ..plots_common import (
    save_or_show,
    beta_iterator,
    ch_tag,
    add_figure_legend_axes_split,
    TWO_COLUMN,
    MPS_left_CUT,
    MPS_right_CUT,
    MPS2_right_end,
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
        "--plot_file_data",
        default=None,
        help=(
            "Where to place the resulting plot showing data. "
            "Default is to output to screen."
        ),
    )
    parser.add_argument(
        "--plot_styles",
        default="styles/paperdraft.mplstyle",
        help="Stylesheet to use for plots",
    )
    return parser.parse_args()


def plot(data, fit_pars):

    data_fig, data_axes =plt.subplots(
        3, 2, num="Summary", figsize=(TWO_COLUMN, 8), layout="constrained"
    )
    
    

    for datum in data:
        if datum[f"gevp_f_ps_E0_chisquare"] > 1.61:
            print(datum["ensemble_name"], datum["beta"], datum["mF"], datum[f"gevp_f_ps_E0_chisquare"])

    for ax, ch in zip(data_axes.ravel(), ["v", "t", "av", "at", "s", "rhoE1"]):
        
        if ch == "rhoE1":
            n = 1
            ch = "v"
            ax.set_ylabel(r"$(\hat{E}^{{\rm V}}_1)^2$")
        else:
            n = 0
            ax.set_ylabel(r"$\hat{m}_{\mathrm{" + ch_tag(ch) + "}}^2$")

        ax.set_xlabel(r"$\hat{m}_{\mathrm{PS}}^2$")
        
        ax.set_xlim(0.0, MPS2_right_end)
        
        #ax.set_ylabel(r"$\hat{m}_{\mathrm{" + ch_tag(ch) + "}}^2 - W \hat{a}$")

        betas = sorted(set([datum["beta"] for datum in data]))
        for beta, colour, marker in beta_iterator(betas):
            to_plot = []
            to_plot_light_alpha = []
            to_plot_light_alpha.append((np.nan, np.nan, np.nan, np.nan))
            for datum in data:
                
                if datum["beta"] != beta:
                    continue
                
                if "w0_samples" not in datum:
                    continue

                if f"gevp_f_{ch}_E{n}_mass_samples" not in datum:
                    print("Missing data for channel!!", ch, f"E{n}", "in ensemble", datum["ensemble_name"])
                    to_plot.append((np.nan, np.nan, np.nan, np.nan))
                    continue    
                
                #for parameter in fit_pars:
                #    if parameter["channel"] == ch:
                #       
                #        W_fit = parameter["Wm_samples"]
    
                
                X = ( datum["w0_samples"] * datum["gevp_f_ps_E0_mass_samples"]) ** 2
                
                Y = ( datum["w0_samples"] * datum[f"gevp_f_{ch}_E{n}_mass_samples"]) ** 2 #- W_fit / datum["w0_samples"]

                if datum[f"gevp_f_{ch}_E{n}_chisquare"] > 1.8:
                    print(datum["ensemble_name"], datum["beta"], datum["mF"], Y.mean, X.mean, datum[f"gevp_f_{ch}_E{n}_chisquare"], ch, n)
                
                if datum[f"gevp_f_{ch}_E{n}_mass_samples"].mean > 1.0:
                    print("Unphysically heavy mass!!", ch, f"E0={datum[f'gevp_f_{ch}_E0_mass_samples'].mean}", "in ensemble", datum["ensemble_name"])
                    to_plot_light_alpha.append((Y.mean, Y.samples.std(), X.mean, X.samples.std()))
                    continue
                
                to_plot.append((Y.mean, Y.samples.std(), X.mean, X.samples.std()))

            y_values, y_errors, x_values, x_errors = zip(*to_plot)
            ax.errorbar(
                x_values,
                y_values,
                xerr=x_errors,
                yerr=y_errors,
                ls="none",
                alpha=1,
                color=colour,
                marker=marker,
                label=f"{beta}",
            )

            y_values, y_errors, x_values, x_errors = zip(*to_plot_light_alpha)
            ax.errorbar(
                x_values,
                y_values,
                xerr=x_errors,
                yerr=y_errors,
                ls="none",
                alpha=0.6,
                color=colour,
                marker=marker,
                #label=f"{beta}",
            )
        

        for parameter in fit_pars:
            if parameter["channel"] == "rhoE1":
                continue

            if parameter["channel"] == ch and n == 0:

                plot_axpb_y(
                    ax,
                    parameter["M_a2_samples"].samples,
                    parameter["Lm_a2_samples"].samples,
                    r"$~\left. \hat{m}_{M,\,\chi}^2(1+L_{M}^m \hat{m}_{\rm PS}^2 )+W_{M}^{m} \hat{a}+R_{M}^{m} \hat{a}^2 \right \vert_{\hat{a}=0}$",
                    0.5,
                    "C9",
                )

                plot_axpb_y(
                    ax,
                    parameter["M_a2_am2_samples"].samples,
                    parameter["Lm_a2_am2_samples"].samples,
                    r"$~\left. \hat{m}_{M,\,\chi}^2(1+L_{M}^m \hat{m}_{\rm PS}^2 )+W_{M}^{m} \hat{a}+R_{M}^{m} \hat{a}^2 + C_M^m \hat{a}\hat{m}_{\rm PS}^2 \right \vert_{\hat{a}=0}$",
                    0.5,
                    "C7",
                )

        ax.set_ylim(0.0, None)
        _, ymax = ax.get_ylim()
        ax.fill_between(
            [0, MPS_left_CUT], [0, 0], [ymax, ymax], color="C6", alpha=0.2)
        
        #ax.fill_between([MPS_right_CUT, MPS2_right_end], [0, 0], [ymax, ymax], color="C6", alpha=0.2)


    add_figure_legend_axes_split(data_fig, data_axes)

    return data_fig


def main():
    args = get_args()
    plt.style.use(args.plot_styles)
    data = read_sample_files(args.data_filenames)
    fit_pars = read_sample_files(args.fit_parameters, group_key="channel")
    data_fig = plot(data, fit_pars)
    save_or_show(data_fig, args.plot_file_data)


if __name__ == "__main__":
    main()
