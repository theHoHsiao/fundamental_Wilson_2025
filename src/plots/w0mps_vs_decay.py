#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from ..plots_common import (
    beta_iterator,
    standard_plot_main,
    add_figure_legend_axes_split,
    ch_tag,
    TWO_COLUMN,
    MPS_left_CUT,
    MPS_right_CUT,
    MPS2_right_end,
    plot_axpb_y,
)


def plot(data, fit_results, **kwargs):
    fig = plt.figure(layout="constrained", figsize=(TWO_COLUMN, 5))
    gs = fig.add_gridspec(nrows=2, ncols=4)
    ax0 = fig.add_subplot(gs[0, :2])
    ax1 = fig.add_subplot(gs[0, 2:])
    ax2 = fig.add_subplot(gs[1, 1:3])
    axs = [ax0, ax1, ax2]
    

    for ch, ax in zip(["ps", "v", "av"], axs):

        ax.set_xlabel(r"$\hat{m}_{\rm PS}^2$")
        ax.set_ylabel(r"$\hat{f}_{\rm " + ch_tag(ch) + "}^2$")
        ax.set_xlim(0.0, MPS2_right_end)

        #ax.set_ylabel(r"$\hat{f}_{\rm " + ch + "}^2 - " r"W \hat{a}$")
        
        betas = sorted(set([datum["beta"] for datum in data]))
        for beta, colour, marker in beta_iterator(betas):
            to_plot = []
            to_plot_light_alpha = []
            to_plot_light_alpha.append((np.nan, np.nan, np.nan, np.nan))
            for datum in data:
                if datum["beta"] != beta:
                    continue

                if "w0_samples" not in datum:
                    print("w0_sample not found in "+datum["ensemble_name"])
                    continue
                if "gevp_f_ps_E0_mass_samples" not in datum:
                    print("gevp" + datum["ensemble_name"])
                    continue

                if f"f_{ch}_decay_constant_samples" not in datum:
                    print("decay constant not found in " + datum["ensemble_name"])
                    continue
                
                X = (datum["w0_samples"] * datum["gevp_f_ps_E0_mass_samples"]) ** 2
                Y = (datum["w0_samples"] * datum[f"f_{ch}_decay_constant_samples"]) ** 2 #- W_fit / datum["w0_samples"]

                

                if datum[f"f_{ch}_chisquare"] > 1.61:
                    print(datum["ensemble_name"], datum["beta"], datum["mF"], Y.mean, X.mean, datum[f"f_{ch}_chisquare"])

                if datum[f"gevp_f_{ch}_E0_mass_samples"].mean > 1.0:
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
                alpha=0.4,
                color=colour,
                marker=marker,
                #label=f"{beta}",
            )

            

        for parameter in fit_results:
            if parameter["channel"] == ch:
            
                plot_axpb_y(
                    ax,
                    parameter["F_a2_samples"].samples,
                    parameter["Lf_a2_samples"].samples,
                    r"$~\left. \hat{f}_{M,\,\chi}^2(1+L_{M}^f \hat{m}_{\rm PS}^2 )+W_{M}^f \hat{a}+R_{M}^f \hat{a}^2 \right \vert_{\hat{a}=0} $",
                    0.5,
                    "C9",
                )

                plot_axpb_y(
                    ax,
                    parameter["F_a2_am2_samples"].samples,
                    parameter["Lf_a2_am2_samples"].samples,
                    r"$~\left. \hat{f}_{M,\,\chi}^2(1+L_{M}^f \hat{m}_{\rm PS}^2 )+W_{M}^f \hat{a}+R_{M}^f \hat{a}^2 + C_M^f \hat{a}\hat{m}_{\rm PS}^2 \right \vert_{\hat{a}=0} $",
                    0.5,
                    "C7",
                )
        
        ax.set_ylim(0.0, None)
        _, ymax = ax.get_ylim()
        ax.fill_between(
            [0, MPS_left_CUT], [0, 0], [ymax, ymax], color="C6", alpha=0.2)
        #ax.fill_between([MPS_right_CUT, MPS2_right_end], [0, 0], [ymax, ymax], color="C6", alpha=0.2)

    

    add_figure_legend_axes_split(fig, axs)
    return fig


if __name__ == "__main__":
    standard_plot_main(plot, fit_results=True, group_key="channel")
