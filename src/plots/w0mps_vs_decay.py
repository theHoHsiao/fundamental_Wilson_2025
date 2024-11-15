#!/usr/bin/env python3

import matplotlib.pyplot as plt
from ..dump import read_sample_files
from ..plots_common import save_or_show, beta_color
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
        "--plot_file",
        default=None,
        help="Where to place the resulting plot. Default is to output to screen.",
    )
    parser.add_argument(
        "--plot_file2",
        default=None,
        help="Where to place the resulting plot. Default is to output to screen.",
    )
    parser.add_argument(
        "--plot_styles",
        default="styles/paperdraft.mplstyle",
        help="Stylesheet to use for plots",
    )
    return parser.parse_args()


def plot_axpb_y(ax, A, L, ch, offset, color, x_i, x_f):
    n_fit = 1000
    Yfit = np.zeros(shape=(A.shape[0], n_fit))

    x = np.linspace(x_i, x_f, n_fit)

    y_up = np.zeros(n_fit)
    y_dn = np.zeros(n_fit)

    for n in range(A.shape[0]):
        Yfit[n] = A[n] * (1 + L[n] * x**2)

    for i in range(n_fit):
        y_err = Yfit[0:-1, i].std()
        y_up[i] = Yfit[:, i].mean() + y_err
        y_dn[i] = Yfit[:, i].mean() - y_err

    # ax.plot(x**2, Yfit[-1], "--", linewidth=0.75, alpha=0.6)
    ax.fill_between(
        x**2, y_up, y_dn, alpha=0.4, label=ch, facecolor=color, edgecolor=None
    )


def plot(data, fit_pars, **kwargs):
    fig, ax1 = plt.subplots(
        1, 2, num="Figure_13_up", figsize=(7, 2.4), layout="constrained"
    )
    fig2, ax2 = plt.subplots(
        1, 1, num="Figure_13_low", figsize=(3.5, 2.4), layout="constrained"
    )

    axs = [ax1[0], ax1[1], ax2]
    subplot_ind = 0

    for ch in ["ps", "v", "av"]:
        ax = ax = axs[subplot_ind]

        ax.set_xlim(0.8, 1.5)
        ax.set_xlabel(r"$\hat{m}_{\rm ps}^2$")
        ax.set_ylabel(r"$\hat{f}_{\rm " + ch + "}^2$")

        betas = sorted(set([datum["beta"] for datum in data]))
        markers = "o^vsx+"
        for beta_idx, (beta, marker) in enumerate(zip(betas, markers)):
            to_plot = []
            for datum in data:
                # print(datum)

                if datum["beta"] != beta:
                    continue

                if "w0_samples" not in datum or "ps_mass_samples" not in datum:
                    continue

                w0_mps = (
                    datum["w0_samples"].samples * datum["ps_mass_samples"].samples
                ) ** 2
                w0_meson = (
                    datum["w0_samples"].samples
                    * datum[f"{ch}_decay_constant_samples"].samples
                ) ** 2

                to_plot.append(
                    (w0_meson.mean(), w0_meson.std(), w0_mps.mean(), w0_mps.std())
                )

            y_values, y_errors, x_values, x_errors = zip(*to_plot)
            ax.errorbar(
                x_values,
                y_values,
                xerr=x_errors,
                yerr=y_errors,
                ls="none",
                alpha=0.7,
                color=beta_color(beta),
                marker=marker,
                label=f"{beta}",
            )

        for parameter in fit_pars:
            if f"F_{ch}_samples" in parameter:
                plot_axpb_y(
                    ax,
                    parameter[f"F_{ch}_samples"].samples,
                    parameter[f"L_{ch}_samples"].samples,
                    "",
                    0,
                    "k",
                    0.894,
                    1.225,
                )

        ax.set_xlim(0.8, 1.5)
        ax.set_ylim(None, None)

        subplot_ind += 1

    handles, labels = fig.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    fig.legend(
        by_label.values(),
        by_label.keys(),
        loc="outside upper center",
        ncol=6,
        borderaxespad=0.2,
    )

    return fig, fig2


def main():
    args = get_args()
    plt.style.use(args.plot_styles)
    data = read_sample_files(args.data_filenames)
    fit_pars = read_sample_files(args.fit_parameters, group_key="channel")
    fig, fig2 = plot(data, fit_pars)
    save_or_show(fig, args.plot_file)
    save_or_show(fig2, args.plot_file2)


if __name__ == "__main__":
    main()
