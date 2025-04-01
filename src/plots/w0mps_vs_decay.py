#!/usr/bin/env python3

import matplotlib.pyplot as plt
from ..plots_common import (
    beta_iterator,
    standard_plot_main,
    add_figure_legend,
    TWO_COLUMN,
)
import numpy as np


def plot_axpb_y(ax, A, L, ch, offset, color):
    n_fit = 1000
    Yfit = np.zeros(shape=(A.shape[0], n_fit))

    x_min, x_max = ax.get_xlim()
    x_i = np.sqrt(x_min)
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
        x**2, y_up, y_dn, alpha=0.4, label=ch, facecolor=color, edgecolor=None
    )


def plot(data, fit_results, **kwargs):
    fig, data_axes = plt.subplots(
        1, 2, num="Figure_12", figsize=(TWO_COLUMN, 4), layout="constrained"
    )

    for ax, ch in zip([data_axes[0], data_axes[1]], ["ps","v"]):

        ax.set_xlabel(r"$\hat{m}_{\rm ps}^2$")
        ax.set_ylabel(r"$\hat{f}_{\rm " + ch + "}^2$")

        betas = sorted(set([datum["beta"] for datum in data]))
        for beta, colour, marker in beta_iterator(betas):
            to_plot = []
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
                Y = (datum["w0_samples"] * datum[f"f_{ch}_decay_constant_samples"]) ** 2

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

        for parameter in fit_results:
            if parameter["channel"] == ch:
                plot_axpb_y(
                    ax,
                    parameter["F_samples"].samples,
                    parameter["L_samples"].samples,
                    "",
                    0,
                    "k",
                )

        #ax.set_xlim(0.8, 1.5)
        ax.set_ylim(None, None)


    add_figure_legend(fig)
    return fig


if __name__ == "__main__":
    standard_plot_main(plot, fit_results=True, group_key="channel")
