#!/usr/bin/env python3

import matplotlib.pyplot as plt

from ..plots_common import standard_plot_main


def plot(data, **kwargs):
    fig, ax = plt.subplots(1, 1, num="Figure_6", figsize=(6, 4), layout="constrained")

    ax.set_ylim(5.5, 6.5)
    ax.set_xlim(0, 0.18)
    ax.set_xlabel(r"$am_{\rm PCAC}$")
    ax.set_ylabel(r"$m_{\rm ps} / f_{\rm ps}$")

    betas = sorted(set([datum["beta"] for datum in data]))
    markers = "s"
    for beta_idx, (beta, marker) in enumerate(zip(betas, markers)):
        to_plot = []
        for datum in data:
            if datum["beta"] != beta:
                continue

            if "ps_mass_samples" not in datum:
                continue

            X = datum["mPCAC_samples"].samples
            Y = (
                datum["ps_mass_samples"].samples
                / datum["ps_decay_constant_samples"].samples
            )

            to_plot.append(
                (
                    datum["ps_mass_samples"].mean
                    / datum["ps_decay_constant_samples"].mean,
                    Y.std(),
                    datum["mPCAC_samples"].mean,
                    X.std(),
                )
            )

        y_values, y_errors, x_values, x_errors = zip(*to_plot)

        ax.errorbar(
            x_values,
            y_values,
            xerr=x_errors,
            yerr=y_errors,
            ls="none",
            alpha=0.7,
            color="b",
            marker="s",
        )

    return fig


if __name__ == "__main__":
    standard_plot_main(plot)
