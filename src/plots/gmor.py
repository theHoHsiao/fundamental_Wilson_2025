#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

from ..plots_common import standard_plot_main
from ..mass import C_R


def plot(data):
    fig, ax = plt.subplots(1, 1, num="Figure_20", figsize=(6, 4), layout="constrained")

    ax.set_xlim(0, 0.16)
    ax.set_xlabel(r"$am_{\rm PCAC}$")
    ax.set_ylabel(r"$(af_{\rm ps})^2  (am_{\rm ps})^2$")

    to_plot = []

    for datum in data:
        if "ps_mass_samples" not in datum:
            continue
        if "mPCAC_samples" not in datum:
            continue

        Z_factor = 1 + 2 * (C_R("ps")) * (8 / datum["beta"]) / (
            16 * np.pi**2 * datum["plaquette_samples"].samples
        )

        X = datum["mPCAC_samples"].samples
        Y = (
            datum["ps_mass_samples"].samples
            * datum["ps_matrix_element_samples"].samples
            * Z_factor
        ) ** 2

        to_plot.append((Y.mean(), Y.std(), X.mean(), X.std()))

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
