#!/usr/bin/env python3

import matplotlib.pyplot as plt
from ..dump import read_sample_files
from ..plots_common import (
    add_figure_legend_axes,
    save_or_show,
    beta_iterator,
    ch_tag,
    get_args_ansatz,
    plot_axpb_y,
    plot_am4pb_y,
    TWO_COLUMN,
    MPS_left_CUT,
    MPS_right_CUT,
)
from argparse import ArgumentParser
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from uncertainties import ufloat


def fit_form(fit_prefix):
    return {
        "a": r"$\hat{m}_{M,\,\chi}^2 (1+L_{M}^m \hat{m}_{\rm PS}^2 )+W_{M}^{m} \hat{a}~~$   with $\hat{a}=0$",
        "a2": r"$\hat{m}_{M,\,\chi}^2(1+L_{M}^m \hat{m}_{\rm PS}^2 )+W_{M}^{m} \hat{a}+R_{M}^{m} \hat{a}^2~~$   with $\hat{a}=0$",
        "m4": r"$\hat{m}_{M,\,\chi}^2(1+L_{M}^m \hat{m}_{\rm PS}^2 + P_M^m \hat{m}_{\rm PS}^4 )+W_{M}^{m} \hat{a}~~$   with $\hat{a}=0$",
        "a2_m4": r"$\hat{m}_{M,\,\chi}^2(1+L_{M}^m \hat{m}_{\rm PS}^2 + P_M^m \hat{m}_{\rm PS}^4 )+W_{M}^{m} \hat{a}+R_{M}^{m} \hat{a}^2~~$   with $\hat{a}=0$",
        "a2_am2": r"$\hat{m}_{M,\,\chi}^2(1+L_{M}^m \hat{m}_{\rm PS}^2 )+W_{M}^{m} \hat{a}+R_{M}^{m} \hat{a}^2 + C_M^m \hat{a}\hat{m}_{\rm PS}^2~~$   with $\hat{a}=0$",
        "full": r"$\hat{m}_{M,\,\chi}^2(1+L_{M}^m \hat{m}_{\rm PS}^2 + P_M^m \hat{m}_{\rm PS}^4 )+W_{M}^{m} \hat{a}+R_{M}^{m} \hat{a}^2 + C_M^m \hat{a}\hat{m}_{\rm PS}^2~~$   with $\hat{a}=0$",
    }[fit_prefix]  


def plot(data, fit_pars, fit_prefix=""):
    
    data_fig, data_axes =plt.subplots(
        3, 2, num="Summary", figsize=(TWO_COLUMN, 8), layout="constrained"
    )
    #data_fig.suptitle(f"Extrapolation with {fit_prefix}")
    
    right_end = 0.45

    for datum in data:
        if datum[f"gevp_f_ps_E0_chisquare"] > 1.61:
            print(datum["ensemble_name"], datum["beta"], datum["mF"], datum[f"gevp_f_ps_E0_chisquare"])

    for ax, ch in zip(data_axes.ravel(), ["v", "t", "av", "at", "s", "rhoE1"]):
        
        if ch == "rhoE1":
            n = 1
            ch = "v"
        else:
            n = 0

        ax.set_xlabel(r"$\hat{m}_{\mathrm{PS}}^2$")
        ax.set_ylabel(r"$\hat{m}_{\mathrm{" + ch_tag(ch) + "}}^2$")
        ax.set_xlim(0.0, right_end)
        
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
                
                X = ( datum["w0_samples"] * datum["gevp_f_ps_E0_mass_samples"]) ** 2
                
                Y = ( datum["w0_samples"] * datum[f"gevp_f_{ch}_E{n}_mass_samples"]) ** 2 #- W_fit / datum["w0_samples"]

                if datum[f"gevp_f_{ch}_E0_chisquare"] > 1.8:
                    print(datum["ensemble_name"], datum["beta"], datum["mF"], Y.mean, X.mean, datum[f"gevp_f_{ch}_E0_chisquare"])
                
                if datum[f"gevp_f_{ch}_E{n}_mass_samples"].mean > 1.0:
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
                alpha=0.6,
                color=colour,
                marker=marker,
                #label=f"{beta}",
            )
        
        for parameter in fit_pars:
            if parameter["channel"] == "rhoE1":
                continue
            
            
            if parameter["channel"] == ch and n == 0:
                if fit_prefix in ["m4", "a2_m4", "full"]:
                    plot_am4pb_y(
                        ax,
                        parameter["M_" + fit_prefix + "_samples"].samples,
                        parameter["Lm_" + fit_prefix + "_samples"].samples,
                        parameter["Pm_" + fit_prefix + "_samples"].samples,
                        fit_form(fit_prefix),
                        0.5,
                        "C7",
                    )
                else:

                    plot_axpb_y(
                        ax,
                        parameter["M_" + fit_prefix + "_samples"].samples,
                        parameter["Lm_" + fit_prefix + "_samples"].samples,
                        fit_form(fit_prefix),
                        0.5,
                        "C7",
                    )

                m_limit = ufloat(parameter["M_" + fit_prefix + "_samples"].mean, parameter["M_" + fit_prefix + "_samples"].samples.std())
                ax.text(0.12, 0.05, f'$\\chi^2/\\mathrm{{dof}} = {parameter["chisquare_" + fit_prefix]:.2f},~~\\hat{{m}}^2_\\chi = {m_limit:.02uSL}$')
            

        ax.set_ylim(0.0, None)
        _, ymax = ax.get_ylim()
        ax.fill_between(
            [0, MPS_left_CUT], [0, 0], [ymax, ymax], color="C6", alpha=0.2)
        ax.fill_between(
            [MPS_right_CUT, right_end], [0, 0], [ymax, ymax], color="C6", alpha=0.2)


    add_figure_legend_axes(data_fig, data_axes[0])    

    return data_fig


def main():
    args = get_args_ansatz()
    plt.style.use(args.plot_styles)
    data = read_sample_files(args.data_filenames)
    fit_pars = read_sample_files(args.fit_parameters, group_key="channel")
    data_fig = plot(data, fit_pars, fit_prefix=args.ansatz)
    save_or_show(data_fig, args.plot_file_data)
    



if __name__ == "__main__":
    main()
