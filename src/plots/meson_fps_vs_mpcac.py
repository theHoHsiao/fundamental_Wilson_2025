#!/usr/bin/env python3

import matplotlib.pyplot as plt

from ..plots_common import (
    standard_plot_main,
    channel_iterator,
    add_figure_legend,
    ONE_COLUMN,
)


def plot(data, **kwargs):
    fig, ax = plt.subplots(
        1, 1, num="Figure_10", figsize=(ONE_COLUMN, 2.4), layout="constrained"
    )

    ax.set_ylim(5, 11)
    ax.set_xlim(0, 0.18)
    ax.set_xlabel(r"$am_{\rm PCAC}$")
    ax.set_ylabel(r"$m_{\rm M} / f_{\rm ps}$")

    channels = ["ps", "v", "av", "at", "s"]
    for channel, colour, marker in channel_iterator(channels):
        to_plot = []

        for datum in data:
            if "mPCAC_samples" not in datum:
                continue

            if f"{channel}_mass_samples" not in datum:
                continue

            X = datum["mPCAC_samples"]
            Y = datum[f"{channel}_mass_samples"] / datum["ps_decay_constant_samples"]

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
            label=channel,
            marker=marker,
        )

    add_figure_legend(fig, title=None)
    return fig


if __name__ == "__main__":
    standard_plot_main(plot)
