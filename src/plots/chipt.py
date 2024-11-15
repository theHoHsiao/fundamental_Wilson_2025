#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

from ..plots_common import standard_plot_main, beta_color
from ..bootstrap import BOOTSTRAP_SAMPLE_COUNT


def plot_YABX(ax, A, B, ch, offset, color):
    n_fit = 1000
    Yfit = np.zeros(shape=(A.shape[0], n_fit))

    x_min, x_max = ax.get_xlim()
    x_i = np.sqrt(x_min)
    x_f = np.sqrt(x_max)
    x = np.linspace(x_i, x_f, n_fit)

    y_up = np.zeros(n_fit)
    y_dn = np.zeros(n_fit)

    for n in range(A.shape[0]):
        Yfit[n] = A[n] + (B[n] * x**2)

    for i in range(n_fit):
        y_err = Yfit[0:-1, i].std()
        y_up[i] = Yfit[-1, i] + y_err
        y_dn[i] = Yfit[-1, i] - y_err

    ax.fill_between(
        x**2, y_up, y_dn, alpha=0.4, label=ch, facecolor=color, edgecolor=None
    )


def plot(data, fit_results, **kwargs):
    fig, ax = plt.subplots(
        1, 1, num="Figure_21", figsize=(3.5, 2.4), layout="constrained"
    )

    ax.set_ylim(0, 0.018)
    ax.set_xlim(0, 0.7)
    ax.set_xlabel(r"$(am_{\rm ps})^2$")
    ax.set_ylabel(r"$(af_{\rm ps})^2$")

    betas = sorted(set([datum["beta"] for datum in data]))
    markers = "o^vsx+"

    for beta_idx, (beta, marker) in enumerate(zip(betas, markers)):
        to_plot = []
        for datum in data:
            if datum["beta"] != beta:
                continue
            if "ps_mass_samples" not in datum:
                continue

            X = (datum["ps_mass_samples"]) ** 2

            Y = (datum["ps_decay_constant_samples"]) ** 2

            to_plot.append((Y.mean, Y.samples.std(), X.mean, X.samples.std()))

        if len(to_plot) < 3:
            continue

        for tmp_result in fit_results:
            if tmp_result["beta"] == str(beta):
                fit_result = tmp_result

        arbitrary_line_width = 0.00005
        plot_YABX(
            ax,
            np.random.normal(
                fit_result[f"A_{beta}_samples"].mean,
                arbitrary_line_width,
                BOOTSTRAP_SAMPLE_COUNT,
            ),
            np.random.normal(
                fit_result[f"B_{beta}_samples"].mean,
                arbitrary_line_width,
                BOOTSTRAP_SAMPLE_COUNT,
            ),  # Construct a distribution that will give a narrow band that looks like a line
            "",
            0,
            beta_color(beta),
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

    handles, labels = fig.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    fig.legend(
        by_label.values(),
        by_label.keys(),
        loc="outside upper center",
        ncol=5,
        borderaxespad=0.2,
    )

    return fig


if __name__ == "__main__":
    standard_plot_main(plot)
