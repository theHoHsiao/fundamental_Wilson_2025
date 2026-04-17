#!/usr/bin/env python3

import matplotlib.pyplot as plt
from ..dump import read_sample_files
from ..plots_common import (
    add_figure_legend_axes,
    save_or_show,
    beta_iterator,
    plot_axpb_y,
    plot_am4pb_y,
    add_figure_legend,
    TWO_COLUMN,
    MPS_left_CUT,
)
from argparse import ArgumentParser
import numpy as np
from uncertainties import ufloat


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
    parser.add_argument(
        "--ansatz",
        default="a",
        help="The ansatz to use for the fit"
    )
    return parser.parse_args()


def fit_form(fit_prefix):
    return {
        "a": r"$\hat{f}_{M,\,\chi}^2 (1+L_{M}^f \hat{f}_{\rm PS}^2 )+W_{M}^f \hat{a}~~$   with $\hat{a}=0$",
        "a2": r"$\hat{f}_{M,\,\chi}^2(1+L_{M}^f \hat{f}_{\rm PS}^2 )+W_{M}^f \hat{a}+R_{M}^f \hat{a}^2~~$   with $\hat{a}=0$",
        "m4": r"$\hat{f}_{M,\,\chi}^2(1+L_{M}^f \hat{f}_{\rm PS}^2 + P_M^f \hat{f}_{\rm PS}^4 )+W_{M}^f \hat{a}~~$   with $\hat{a}=0$",
        "a2_m4": r"$\hat{f}_{M,\,\chi}^2(1+L_{M}^f \hat{f}_{\rm PS}^2 + P_M^f \hat{f}_{\rm PS}^4 )+W_{M}^f \hat{a}+R_{M}^f \hat{a}^2~~$   with $\hat{a}=0$",
        "a2_am2": r"$\hat{f}_{M,\,\chi}^2(1+L_{M}^f \hat{f}_{\rm PS}^2 )+W_{M}^f \hat{a}+R_{M}^f \hat{a}^2 + C_M^f \hat{a}\hat{f}_{\rm PS}^2~~$   with $\hat{a}=0$",
        "full": r"$\hat{f}_{M,\,\chi}^2(1+L_{M}^f \hat{f}_{\rm PS}^2 + P_M^f \hat{f}_{\rm PS}^4 )+W_{M}^f \hat{a}+R_{M}^f \hat{a}^2 + C_M^f \hat{a}\hat{f}_{\rm PS}^2~~$   with $\hat{a}=0$",
    }[fit_prefix]  


def plot(data, fit_pars, fit_prefix=""):
    
    fig = plt.figure(layout="constrained", figsize=(TWO_COLUMN, 5))
    gs = fig.add_gridspec(nrows=2, ncols=4)
    ax0 = fig.add_subplot(gs[0, :2])
    ax1 = fig.add_subplot(gs[0, 2:])
    ax2 = fig.add_subplot(gs[1, 1:3])
    axs = [ax0, ax1, ax2]
    
    right_end = 1

    for ax, ch in zip(axs, ["ps", "v", "av"]):

        ax.set_xlabel(r"$\hat{f}_{\rm PS}^2$")
        ax.set_ylabel(r"$\hat{f}_{\rm " + ch + "}^2$")

        ax.set_xlim(0.0, right_end)
        
        #ax.set_ylabel(r"$\hat{f}_{\mathrm{" + ch_tag(ch) + "}}^2 - W \hat{a}$")

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

                if f"f_{ch}_decay_constant_samples" not in datum:
                    print(f"{ch} decay constant not found in " + datum["ensemble_name"])
                    continue
                
                X = (datum["w0_samples"] * datum["gevp_f_ps_E0_mass_samples"]) ** 2
                Y = (datum["w0_samples"] * datum[f"f_{ch}_decay_constant_samples"]) ** 2 
            

                if datum[f"f_{ch}_chisquare"] > 1.61:
                    print(datum["ensemble_name"], datum["beta"], datum["mF"], Y.mean, X.mean, datum[f"f_{ch}_chisquare"])

                if datum[f"gevp_f_{ch}_E0_mass_samples"].mean > 1.0:
                    #print("Unphysically heavy mass!!", ch, f"E{n}={datum[f'gevp_f_{ch}_E{n}_mass_samples'].mean}", "in ensemble", datum["ensemble_name"])
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
                alpha=0.4,
                color=colour,
                marker=marker,
                #label=f"{beta}",
            )

        
        for parameter in fit_pars:
            if parameter["channel"] == "rhoE1":
                continue
            
            
            if parameter["channel"] == ch:
                if fit_prefix in ["m4", "a2_m4", "full"]:
                    plot_am4pb_y(
                        ax,
                        parameter["F_" + fit_prefix + "_samples"].samples,
                        parameter["Lf_" + fit_prefix + "_samples"].samples,
                        parameter["Pf_" + fit_prefix + "_samples"].samples,
                        fit_form(fit_prefix),
                        0.5,
                        "C7",
                    )
                else:

                    plot_axpb_y(
                        ax,
                        parameter["F_" + fit_prefix + "_samples"].samples,
                        parameter["Lf_" + fit_prefix + "_samples"].samples,
                        fit_form(fit_prefix),
                        0.5,
                        "C7",
                    )

                m_limit = ufloat(parameter["F_" + fit_prefix + "_samples"].mean, parameter["F_" + fit_prefix + "_samples"].samples.std())
                ax.text(0.12, 0.002, f'$\\chi^2/\\mathrm{{dof}} = {parameter["chisquare_" + fit_prefix]:.2f},~~\\hat{{f}}^2_\\chi = {m_limit:.02uSL}$')
            

        ax.set_ylim(0.0, None)
        _, ymax = ax.get_ylim()
        ax.fill_between(
            [0, MPS_left_CUT], [0, 0], [ymax, ymax], color="C6", alpha=0.2)
        #ax.fill_between([MPS_right_CUT, right_end], [0, 0], [ymax, ymax], color="C6", alpha=0.2)


    add_figure_legend(fig)
    
    #plt.tight_layout()
    

    return fig


def main():
    args = get_args()
    plt.style.use(args.plot_styles)
    data = read_sample_files(args.data_filenames)
    fit_pars = read_sample_files(args.fit_parameters, group_key="channel")
    
    data_fig = plot(data, fit_pars, fit_prefix=args.ansatz)
    save_or_show(data_fig, args.plot_file_data)


if __name__ == "__main__":
    main()
