#!/usr/bin/env python3

import matplotlib.pyplot as plt

from ..plots_common import standard_plot_main


def plot(data, **kwargs):
    fig, ax = plt.subplots(
        1, 1, num="Figure_20", figsize=(3.5, 2.4), layout="constrained"
    )

    ax.set_xlim(0, 0.16)
    ax.set_xlabel(r"$am_{\rm PCAC}$")
    ax.set_ylabel(r"$(af_{\rm ps})^2  (am_{\rm ps})^2$")

    to_plot = []

    for datum in data:
        if "ps_mass_samples" not in datum:
            continue
        if "mPCAC_samples" not in datum:
            continue

        X = datum["mPCAC_samples"]
        Y = (datum["ps_mass_samples"] * datum["ps_decay_constant_samples"]) ** 2

        to_plot.append((Y.mean, Y.samples.std(), X.mean, X.samples.std()))

    y_values, y_errors, x_values, x_errors = zip(*to_plot)

    # print(np.array(to_fit_x).shape)

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
