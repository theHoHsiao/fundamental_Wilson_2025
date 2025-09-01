#!/usr/bin/env python3

import matplotlib.pyplot as plt
from ..dump import read_sample_files
from ..plots_common import (
    save_or_show,
    beta_iterator,
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
        "--plot_file_data",
        default=None,
        help=(
            "Where to place the resulting plot showing data. "
            "Default is to output to screen."
        ),
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


def plot_axpb_y(ax, A, L, ch, alpha, color):
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
        x**2, y_up, y_dn, alpha=alpha, label=ch, facecolor=color, edgecolor=None
    )


def plot(data, fit_pars):
    data_fig, data_axes = plt.subplots(
        3, 2, num="Figure_12", figsize=(TWO_COLUMN, 8), layout="constrained"
    )

    for datum in data:
        print(datum["ensemble_name"], datum["beta"], datum["mF"])


    for ax, ch in zip(data_axes.ravel(), ["v", "t", "av", "at", "s"]):
        
        print(f"~~~~~~~~~~~~~~~~{ch}~~~~~~~~~~~~~~~~")

        ax.set_xlabel(r"$\hat{m}_{\mathrm{ps}}^2$")
        ax.set_ylabel(r"$\hat{m}_{\mathrm{" + ch_tag(ch) + "}}^2$")
        
        #ax.set_ylabel(r"$\hat{m}_{\mathrm{" + ch_tag(ch) + "}}^2 - W \hat{a}$")


        for parameter in fit_pars:
            if parameter["channel"] == ch:
                L_fit = parameter["L_samples"]
                W_fit = parameter["W_samples"]
                M_fit = parameter["M_samples"]

        betas = sorted(set([datum["beta"] for datum in data]))
        for beta, colour, marker in beta_iterator(betas):
            to_plot = []
            for datum in data:
                
                if datum["beta"] != beta:
                    continue
                
                if "w0_samples" not in datum:
                    continue
                
                if datum[f"gevp_f_{ch}_E0_chisquare"] > 1.6:
                    print(f"gevp_f_{ch}_E0_chisquare > 1.6 in {datum['ensemble_name']}")
                
                X = ( datum["w0_samples"] * datum["gevp_f_ps_E0_mass_samples"]) ** 2
                Y = ( datum["w0_samples"] * datum[f"gevp_f_{ch}_E0_mass_samples"]) ** 2 #- W_fit / datum["w0_samples"]

                print(datum["ensemble_name"], datum["beta"], datum["mF"], Y.mean, X.mean, datum[f"gevp_f_{ch}_E0_chisquare"])
                
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
        
        for parameter in fit_pars:
            if parameter["channel"] == ch:
                plot_axpb_y(
                    ax,
                    parameter["M_samples"].samples,
                    parameter["L_samples"].samples,
                    "",
                    0.4,
                    "k",
                )

        ax.set_xlim(0.0, 0.45)

    #add_figure_legend(data_fig)
    add_figure_legend_axes(data_fig, data_axes)

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
