#!/usr/bin/env python3

import matplotlib.pyplot as plt

from ..plots_common import standard_plot_main


def beta_color(b):
    return {
        6.6: "k",
        6.65: "r",
        6.7: "b",
        6.75: "m",
        6.8: "g",
        6.9: "tab:brown",
    }.get(b, b)


def plot(data):
    fig, axs = plt.subplots(
        3, 2, num="Figure_12", figsize=(12, 12), layout="constrained"
    )
    subplot_row = 0
    subplot_col = 0

    for ch in ["v", "t", "s", "av", "at", "rhoE1"]:
        ax = axs[subplot_row, subplot_col]
        ax.set_xlabel(r"$\hat{m}_{\mathrm{ps}}^2$")
        ax.set_ylabel(r"$\hat{m}_{\mathrm{" + ch + "}}^2$")

        betas = sorted(set([datum["beta"] for datum in data]))
        markers = "o^vsx+"
        for beta_idx, (beta, marker) in enumerate(zip(betas, markers)):
            to_plot = []
            for datum in data:
                # print(datum)

                if datum["beta"] != beta:
                    continue

                if "w0_samples" not in datum or "smear_ps_mass_samples" not in datum:
                    continue

                w0_mps = (
                    datum["w0_samples"].samples * datum["smear_ps_mass_samples"].samples
                ) ** 2
                w0_meson = (
                    datum["w0_samples"].samples
                    * datum[f"smear_{ch}_mass_samples"].samples
                ) ** 2
                print(w0_meson.mean(), w0_meson.std(), w0_mps.mean(), w0_mps.std())

                to_plot.append(
                    (w0_meson.mean(), w0_meson.std(), w0_mps.mean(), w0_mps.std())
                )

            y_values, y_errors, x_values, x_errors = zip(*to_plot)
            ax.errorbar(
                x_values,
                y_values,
                xerr=x_errors,
                yerr=y_errors,
                ls="none",
                color=beta_color(beta),
                marker=marker,
                label=f"{beta}",
            )

        ax.set_xlim(0.8, 1.5)
        ax.set_ylim(None, None)

        subplot_row += 1 - subplot_col * 3
        subplot_col += int(subplot_row / 3)
        subplot_row = subplot_row % 3

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
