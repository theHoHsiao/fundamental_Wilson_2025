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
    MPS_right_CUT,
    plot_axpb_y,
    plot_am4pb_y,
)
from argparse import ArgumentParser
import numpy as np



def get_fps_fit(A, L,):
    n_fit = 1000
    Yfit = np.zeros(shape=(A.shape[0], n_fit))

    x_i = np.sqrt(MPS_left_CUT)
    x_f = np.sqrt(MPS_right_CUT)
    x = np.linspace(x_i, x_f, n_fit)

    for n in range(A.shape[0]):
        Yfit[n] = A[n] * (1 + L[n] * x**2)
    

    return Yfit


def plot_ps_ths(ax, times_of_mps, label, alpha, scale=None):
    n_fit=1000
    x_i = np.sqrt(MPS_left_CUT)
    x_f = np.sqrt(MPS_right_CUT)
    m_ps = np.linspace(x_i, x_f, n_fit)
    X = m_ps**2

    Y = (times_of_mps * m_ps)**2
    if scale is not None:
        Y = Y / scale.mean(axis=0)
        X = X / scale.mean(axis=0)

    ax.plot(X, Y, "--", color="k", linewidth=1, label=label, alpha=alpha)



def plot(fit_pars, ansatz=""):
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
            fm_ps = parameter[f"F_{ansatz}_samples"].samples
            lf_ps = parameter[f"Lf_{ansatz}_samples"].samples
    
    
    fps_scale = get_fps_fit(fm_ps, lf_ps)
    fps_scale_mean = fps_scale.mean(axis=0)

    plot_ps_ths(mass_ax, 1, r"$m_{\rm PS}^2 / f_{\rm PS}^2$", 0.8, scale=fps_scale)
    plot_ps_ths(mass_ax, 2, r"$4m_{\rm PS}^2 / f_{\rm PS}^2$", 0.6, scale=fps_scale)
    plot_ps_ths(mass_ax, 3, r"$9m_{\rm PS}^2 / f_{\rm PS}^2$", 0.4, scale=fps_scale)
    plot_ps_ths(mass_ax, 4, r"$16m_{\rm PS}^2 / f_{\rm PS}^2$", 0.3, scale=fps_scale)

    mass_ax.set_ylim(0, 240)
    mass_ax.set_xlim(MPS_left_CUT/ fps_scale_mean[0], 0.4 / fps_scale_mean[-1])

    decay_ax.set_xlim(MPS_left_CUT/ fps_scale_mean[0], 0.4 / fps_scale_mean[-1])
    decay_ax.set_ylim(0, 7)
    
   
    for ch in ["v", "t", "s","av", "at"]:
        
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
                        scale=fps_scale,
                    )
                else:

                    plot_axpb_y(
                        mass_ax,
                        parameter[f"M_{ansatz}_samples"].samples,
                        parameter[f"Lm_{ansatz}_samples"].samples,
                        r"$ \rm " + ch_tag(ch) + "$",
                        0.5,
                        channel_color(ch),
                        scale=fps_scale,
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
                        scale=fps_scale,
                    )

                else:

                    plot_axpb_y(
                        decay_ax,
                        parameter[f"F_{ansatz}_samples"].samples,
                        parameter[f"Lf_{ansatz}_samples"].samples,
                        r"$ \rm " + ch_tag(ch) + "$",
                        0.5,
                        channel_color(ch),
                        scale=fps_scale,
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
