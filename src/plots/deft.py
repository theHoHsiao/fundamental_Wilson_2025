#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

from ..bootstrap import BOOTSTRAP_SAMPLE_COUNT
from ..plots_common import (
    standard_plot_main,
    beta_iterator,
    add_figure_legend,
    ONE_COLUMN,
)


def plot_axpb_y_minus(ax, A, B, ch, offset, color):
    n_fit = 1000
    Yfit = np.zeros(shape=(A.shape[0], n_fit))

    x_min, x_max = ax.get_xlim()
    x_i = -np.sqrt(abs(x_min))
    x_f = -np.sqrt(abs(x_max))
    x = np.linspace(x_i, x_f, n_fit)

    y_up = np.zeros(n_fit)
    y_dn = np.zeros(n_fit)

    for n in range(A.shape[0]):
        Yfit[n] = A[n] + (B[n] * -(x**2))

    for i in range(n_fit):
        y_err = Yfit[0:-1, i].std()
        y_up[i] = Yfit[-1, i] + y_err
        y_dn[i] = Yfit[-1, i] - y_err

    ax.fill_between(
        -(x**2), y_up, y_dn, alpha=0.4, label=ch, facecolor=color, edgecolor=None
    )


def plot(data, fit_results, **kwargs):
    fig, ax = plt.subplots(
        1, 1, num="Figure_22", figsize=(ONE_COLUMN, 2.4), layout="constrained"
    )

    ax.set_ylim(0.4, 1.7)
    ax.set_xlim(-7, -4)
    ax.set_xlabel(r"$\log  [ (af_{\rm ps})^2  ] $")
    ax.set_ylabel(r"$\log [(am_{\rm ps})^2 / (am_{\rm PCAC})]$")

    betas = sorted(set([datum["beta"] for datum in data]))

    for beta, colour, marker in beta_iterator(betas):
        to_plot = []
        for datum in data:
            if datum["beta"] != beta:
                continue
            if "ps_mass_samples" not in datum:
                continue

            Y = (datum["ps_mass_samples"]) ** 2 / datum["mPCAC_samples"]

            X = (datum["ps_decay_constant_samples"]) ** 2

            X = np.log(X)
            Y = np.log(Y)

            to_plot.append((Y.mean, Y.samples.std(), X.mean, X.samples.std()))

        if len(to_plot) < 3:
            continue

        for tmp_result in fit_results:
            if tmp_result["beta"] == beta:
                fit_result = tmp_result

        arbitrary_line_width = 0.0005
        plot_axpb_y_minus(
            ax,
            np.random.normal(
                fit_result["A_samples"].mean,
                arbitrary_line_width,
                BOOTSTRAP_SAMPLE_COUNT,
            ),
            np.random.normal(
                fit_result["B_samples"].mean,
                arbitrary_line_width,
                BOOTSTRAP_SAMPLE_COUNT,
            ),  # line with a small band
            "",
            0,
            colour,
        )

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

    add_figure_legend(fig)
    return fig


if __name__ == "__main__":
    standard_plot_main(plot)
