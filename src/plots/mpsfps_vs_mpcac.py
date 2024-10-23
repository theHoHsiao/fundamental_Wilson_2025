#!/usr/bin/env python3

import matplotlib.pyplot as plt
from ..plots_common import standard_plot_main
from ..mass import C_R


def plot(data):
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

            Z_factor = 1 + 2 * (C_R("PS")) * (8 / beta) / (
                16 * 3.141592653589793**2 * datum["plaquette_samples"].samples
            )

            X = datum["mPCAC_samples"].samples
            Y = (
                datum["ps_mass_samples"].samples
                / datum["ps_matrix_element_samples"].samples
                / Z_factor
            )

            print(datum["ps_matrix_element_samples"].mean)

            to_plot.append(
                (
                    datum["ps_mass_samples"].mean
                    / datum["ps_matrix_element_samples"].mean
                    / Z_factor.mean(),
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
