#!/usr/bin/env python3

import matplotlib.pyplot as plt
from ..dump import read_sample_files
from ..plots_common import (
    save_or_show,
    get_args_ansatz,
    ch_tag,
    channel_color,
    add_figure_legend_axes,
    TWO_COLUMN,
    MPS_left_CUT,
    plot_axpb_y,
    plot_am4pb_y,
)

def plot(fit_pars, ansatz=""):
    summary_fig, mass_axs = plt.subplots(
        1, 2, num="Summary", figsize=(TWO_COLUMN, 4), layout="constrained"
    )
    mass_ax = mass_axs[0]
    mass_ax.plot([0, 1], [0, 1], "--", color="k", linewidth=1, label=r"$\hat{m}_{\rm PS}^2$", alpha=0.8)
    mass_ax.plot([0, 1], [0, 4], "--", color="k", linewidth=1, label=r"$ 4\hat{m}_{\rm PS}^2$", alpha=0.6)
    mass_ax.plot([0, 1], [0, 9], "--", color="k", linewidth=1, label=r"$ 9\hat{m}_{\rm PS}^2$", alpha=0.4)
    mass_ax.plot([0, 1], [0, 16], "--", color="k", linewidth=1, label=r"$ 16\hat{m}_{\rm PS}^2$", alpha=0.3)

    mass_ax.set_ylim(0, 3.5)
    mass_ax.set_xlim(MPS_left_CUT, 0.4)
    mass_ax.set_xlabel(r"$\hat{m}_{\rm PS}^2$")
    mass_ax.set_ylabel(r"$\hat{m}_{\rm M}^2$")

    decay_ax = mass_axs[1] #mass_ax.twinx() 
    decay_ax.set_ylabel(r"$\hat{f}_{\rm M}^2$")
    decay_ax.set_xlabel(r"$\hat{m}_{\rm PS}^2$")
    decay_ax.set_xlim(MPS_left_CUT, 0.4)
    decay_ax.set_ylim(0, 0.08)

   
    for ch in ["v", "t", "s", "av", "at"]:
        
        for parameter in fit_pars:

            if f"M_{ansatz}_samples" not in parameter:
                    continue
            
            if parameter["channel"] == ch:

                if ansatz in ["m4", "a2_m4", "full"]:
                    plot_am4pb_y(
                        mass_ax,
                        parameter[f"M_{ansatz}_samples"].samples,
                        parameter[f"Lm_{ansatz}_samples"].samples,
                        parameter[f"Pm_{ansatz}_samples"].samples,
                        r"$ \rm " + ch_tag(ch) + "$",
                        0.5,
                        channel_color(ch),
                    )
                else:

                    plot_axpb_y(
                        mass_ax,
                        parameter[f"M_{ansatz}_samples"].samples,
                        parameter[f"Lm_{ansatz}_samples"].samples,
                        r"$ \rm " + ch_tag(ch) + "$",
                        0.5,
                        channel_color(ch),
                    )
    

    for ch in ["ps", "v", "av"]:
        
        for parameter in fit_pars:
            if parameter["channel"] == ch:

                if ansatz in ["m4", "a2_m4", "full"]:

                    plot_am4pb_y(
                        decay_ax,
                        parameter[f"F_{ansatz}_samples"].samples,
                        parameter[f"Lf_{ansatz}_samples"].samples,
                        parameter[f"Pf_{ansatz}_samples"].samples,
                        r"$ \rm " + ch_tag(ch) + "$",
                        0.5,
                        channel_color(ch),
                    )
                else:

                    plot_axpb_y(
                        decay_ax,
                        parameter[f"F_{ansatz}_samples"].samples,
                        parameter[f"Lf_{ansatz}_samples"].samples,
                        r"$ \rm " + ch_tag(ch) + "$",
                        0.5,
                        channel_color(ch),
                    )


    add_figure_legend_axes(summary_fig, mass_axs, 10, title=None)


    return summary_fig


def main():
    args = get_args_ansatz()
    plt.style.use(args.plot_styles)
    fit_pars = read_sample_files(args.fit_parameters, group_key="channel")
    summary_fig = plot(fit_pars, ansatz=args.ansatz)
    save_or_show(summary_fig, args.plot_file_data)


if __name__ == "__main__":
    main()
