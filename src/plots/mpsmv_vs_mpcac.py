#!/usr/bin/env python3

import matplotlib.pyplot as plt
from ..plots_common import standard_plot_main, beta_iterator, ONE_COLUMN


def plot(data, **kwargs):
    fig, ax = plt.subplots(
        1, 1, num="Figure_6", figsize=(ONE_COLUMN, 2.4), layout="constrained"
    )

    ax.set_ylim(0.78, 0.94)
    ax.set_xlim(0, 0.18)
    ax.set_xlabel(r"$am_{\rm PCAC}$")
    ax.set_ylabel(r"$m_{\rm ps} / m_{\rm v}$")

    betas = sorted(set([datum["beta"] for datum in data]))
    for beta, colour, marker in beta_iterator(betas):
        to_plot = []
        for datum in data:
            if datum["beta"] != beta:
                continue

            if "ps_mass_samples" not in datum:
                continue

            X = datum["mPCAC_samples"]
            Y = datum["ps_mass_samples"] / datum["v_mass_samples"]

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
        )

    return fig


if __name__ == "__main__":
    standard_plot_main(plot)
