#!/usr/bin/env python3

import matplotlib.pyplot as plt

from ..plots_common import standard_plot_main


def plot(data, external_data, fit_results):
    fig, ax = plt.subplots(
        1, 1, num="Figure_15", figsize=(3.5, 2.4), layout="constrained"
    )

    ax.set_xlim(0, 1.5)
    ax.set_xlabel(r"$\hat{m}_{\rm ps(PS)}^2$")
    ax.set_ylabel(r"$m_{\rm v(V)} / f_{\rm ps(PS)}$")

    ax.errorbar(
        external_data.get("R_value").values,
        external_data.get("L_value").values * 2**0.5,
        xerr=external_data.get("R_uncertainty").values,
        yerr=external_data.get("L_uncertainty").values * 2**0.5,
        linestyle="",
        marker="o",
        markerfacecolor="none",
        elinewidth=1,
        capthick=1,
        capsize=1,  # zorder=5,
        color="r",
        alpha=0.7,
        label=r"$N_{\rm f}=2$ $Sp(4)$",
    )

    to_plot = []

    for datum in data:
        if "ps_mass_samples" not in datum:
            continue
        if "v_mass_samples" not in datum:
            continue

        w0 = datum["w0_samples"]

        X = (datum["ps_mass_samples"] * w0) ** 2
        Y = (
            datum["v_mass_samples"] / datum["ps_decay_constant_samples"]
            - fit_results[0]["W_vdfps_samples"] / w0
        )

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
        label=r"$N_{\rm as}=3$ $Sp(4)$",
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
