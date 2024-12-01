#!/usr/bin/env python3

import matplotlib.pyplot as plt
from ..plots_common import (
    standard_plot_main,
    beta_iterator,
    add_figure_legend,
    ONE_COLUMN,
)


def plot(data, **kwargs):
    fig, ax = plt.subplots(
        1, 1, num="Figure_17", figsize=(ONE_COLUMN, 2.4), layout="constrained"
    )

    ax.set_ylim(1, 5)
    ax.set_xlim(0.79, 1.61)
    ax.set_xlabel(r"$\hat{m}_{\rm ps}^2$")
    ax.set_ylabel(r"$w_0 / a$")

    to_plot = []
    betas = sorted(set([datum["beta"] for datum in data]))
    for beta, colour, marker in beta_iterator(betas):
        to_plot = []
        for datum in data:
            if datum["beta"] != beta:
                continue
            if "ps_mass_samples" not in datum:
                continue
            if "w0_samples" not in datum:
                continue

            w0 = datum["w0_samples"]

            X = (datum["ps_mass_samples"] * w0) ** 2

            Y = w0

            to_plot.append((Y.mean, Y.samples.std(), X.mean, X.samples.std()))

        if not to_plot:
            continue

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
