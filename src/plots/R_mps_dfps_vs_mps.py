#!/usr/bin/env python3

import matplotlib.pyplot as plt
from ..plots_common import standard_plot_main, beta_color
from ..mass import C_R


def plot(data):
    fig, ax = plt.subplots(1, 1, num="Figure_17", figsize=(6, 4), layout="constrained")

    ax.set_xlim(0.79, 1.61)
    ax.set_xlabel(r"$\hat{m}_{\rm ps}^2$")
    ax.set_ylabel(r"$m_{\rm ps} / f_{\rm ps}$")

    to_plot = []
    betas = sorted(set([datum["beta"] for datum in data]))
    markers = "o^vsx+"
    for beta_idx, (beta, marker) in enumerate(zip(betas, markers)):
        to_plot = []
        for datum in data:
            if datum["beta"] != beta:
                continue
            if "ps_mass_samples" not in datum:
                continue
            if "w0_samples" not in datum:
                continue

            w0 = datum["w0_samples"].samples

            Z_factor = 1 + 2 * (C_R("ps")) * (8 / datum["beta"]) / (
                16 * 3.141592653589793**2 * datum["plaquette_samples"].samples
            )

            X = (datum["ps_mass_samples"].samples * w0) ** 2

            Y = (
                datum["ps_mass_samples"].samples
                / datum["ps_matrix_element_samples"].samples
                / Z_factor
            )

            to_plot.append((Y.mean(), Y.std(), X.mean(), X.std()))

        if not to_plot:
            continue

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
        ncol=3,
        borderaxespad=0.2,
    )

    return fig


if __name__ == "__main__":
    standard_plot_main(plot)
