#!/usr/bin/env python3

import matplotlib.pyplot as plt

from ..plots_common import standard_plot_main, channel_color


def plot(data):
    fig, ax = plt.subplots(1, 1, num="Figure_10", figsize=(6, 4), layout="constrained")

    ax.set_ylim(5, 11)
    ax.set_xlim(0, 0.18)
    ax.set_xlabel(r"$am_{\rm PCAC}$")
    ax.set_ylabel(r"$m_{\rm M} / f_{\rm ps}$")

    channels = ["ps", "v", "av", "at", "s"]
    markers = "o^sx+"

    for channel_idx, (channel, marker) in enumerate(zip(channels, markers)):
        to_plot = []

        for datum in data:
            if "mPCAC_samples" not in datum:
                continue

            if f"{channel}_mass_samples" not in datum:
                continue

            X = datum["mPCAC_samples"].samples
            Y = (
                datum[f"{channel}_mass_samples"].samples
                / datum["ps_decay_constant_samples"].samples
            )

            to_plot.append((Y.mean(), Y.std(), X.mean(), X.std()))

        y_values, y_errors, x_values, x_errors = zip(*to_plot)

        ax.errorbar(
            x_values,
            y_values,
            xerr=x_errors,
            yerr=y_errors,
            ls="none",
            alpha=0.7,
            color=channel_color(channel),
            label=channel,
            marker=marker,
        )

    handles, labels = fig.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    fig.legend(
        by_label.values(),
        by_label.keys(),
        loc="outside upper center",
        ncol=6,
        borderaxespad=0.2,
    )

    return fig


if __name__ == "__main__":
    standard_plot_main(plot)
