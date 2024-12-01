#!/usr/bin/env python3

import matplotlib.pyplot as plt
from ..plots_common import (
    standard_plot_main,
    channel_iterator,
    add_figure_legend,
    TWO_COLUMN,
)


def plot(data, **kwargs):
    fig, ax = plt.subplots(
        1, 2, num="Figure_8", figsize=(TWO_COLUMN, 2.4), layout="constrained"
    )

    ax[0].set_ylim(0, 1.4)
    ax[0].set_xlim(-1.08, -1.01)
    ax[0].set_xlabel(r"$am_0$")
    ax[0].set_ylabel(r"$am_{\rm M}$")

    ax[1].set_ylim(0, 1.4)
    ax[1].set_xlim(0, 0.17)
    ax[1].set_xlabel(r"$am_{\rm PCAC}$")
    ax[1].set_ylabel(r"$am_{\rm M}$")

    channels = ["ps", "v", "av", "at", "s"]
    for channel, colour, marker in channel_iterator(channels):
        to_plot = []
        bare_mass = []
        for datum in data:
            if "mPCAC_samples" not in datum:
                continue

            if f"{channel}_mass_samples" not in datum:
                continue

            X = datum["mPCAC_samples"]
            Y = datum[f"{channel}_mass_samples"]

            bare_mass.append(datum["mAS"])

            to_plot.append((Y.mean, Y.samples.std(), X.mean, X.samples.std()))

        y_values, y_errors, x_values, x_errors = zip(*to_plot)

        ax[0].errorbar(
            bare_mass,
            y_values,
            yerr=y_errors,
            ls="none",
            alpha=1,
            color=colour,
            marker=marker,
        )

        ax[1].errorbar(
            x_values,
            y_values,
            xerr=x_errors,
            yerr=y_errors,
            ls="none",
            alpha=1,
            color=colour,
            label=channel,
            marker=marker,
        )

    add_figure_legend(fig, title=None)
    return fig


if __name__ == "__main__":
    standard_plot_main(plot)
