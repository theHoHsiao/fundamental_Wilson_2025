#!/usr/bin/env python3

import matplotlib.pyplot as plt
from ..plots_common import standard_plot_main, ch_tag


def plot(data, **kwargs):
    fig, axs = plt.subplots(
        2, 1, num="Figure_2a", figsize=(3.5, 4.8), layout="constrained"
    )
    if len(set([(datum["beta"], datum["mAS"]) for datum in data])) > 1:
        raise NotImplementedError("Inconsistent physical parameters found.")

    Ns_max = max([datum["Ns"] for datum in data])
    Ns_max_datum = [datum for datum in data if datum["Ns"] == Ns_max]
    if len(Ns_max_datum) > 1:
        raise ValueError("Multiple ensembles have the same largest volume.")
    m_ps_inf = Ns_max_datum[0]["ps_mass_samples"]

    ch_i = 0
    for ch in ["ps", "v"]:
        ax = axs[ch_i]

        ax.set_xlabel(r"$m_{\rm ps}^{\rm inf} L$")
        ax.set_ylabel(r"$am_{\mathrm{" + ch_tag(ch) + "}}$")

        to_plot = []
        for datum in data:
            if "ps_mass_samples" not in datum:
                continue

            X = m_ps_inf * datum["Ns"]
            Y = datum[f"{ch}_mass_samples"]

            to_plot.append((Y.mean, Y.samples.std(), X.mean, X.samples.std()))

        y_values, y_errors, x_values, x_errors = zip(*to_plot)

        ax.errorbar(
            x_values,
            y_values,
            xerr=x_errors,
            yerr=y_errors,
            ls="none",
            alpha=1,
            color="C0",
            marker="s",
        )

        ch_i += 1
    return fig


if __name__ == "__main__":
    standard_plot_main(plot)
