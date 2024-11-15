#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from ..plots_common import standard_plot_main, beta_color


def plot(data, **kwargs):
    fig, ax = plt.subplots(1, 1, num="6", figsize=(3.5, 2.4), layout="constrained")

    ratio_QCD = np.random.normal(1465, 25, 200) / np.random.normal(775.26, 0.25, 200)
    plt.fill_between(
        [0.01, 0.1],
        [ratio_QCD.mean() - ratio_QCD.std(), ratio_QCD.mean() - ratio_QCD.std()],
        [ratio_QCD.mean() + ratio_QCD.std(), ratio_QCD.mean() + ratio_QCD.std()],
        label=r"QCD $\rho(1450)/\rho(770)$",
        alpha=0.7,
        color="dodgerblue",
    )

    # QCD data from https://pdg.lbl.gov/2019/listings/rpp2019-list-rho-770.pdf
    # QCD data from https://pdg.lbl.gov/2014/listings/rpp2014-list-rho-1450.pdf

    ax.set_xlim(-0.01, 1.51)
    ax.set_xlabel(r"$\hat{m}_{\rm ps}^2$")
    ax.set_ylabel(r"$m_{\rm v^\prime} / m_{\rm v}$")

    to_plot = []
    betas = sorted(set([datum["beta"] for datum in data]))
    markers = "o^vsx+"
    for beta_idx, (beta, marker) in enumerate(zip(betas, markers)):
        to_plot = []
        for datum in data:
            if datum["beta"] != beta:
                continue
            if "smear_ps_mass_samples" not in datum:
                continue
            if "smear_v_mass_samples" not in datum:
                continue
            if "smear_rhoE1_mass_samples" not in datum or np.isnan(
                (datum["smear_rhoE1_mass_samples"].samples).mean()
            ):
                continue

            w0 = datum["w0_samples"].samples

            X = (datum["smear_ps_mass_samples"].samples * w0) ** 2

            Y = (
                datum["smear_rhoE1_mass_samples"].samples
                / datum["smear_v_mass_samples"].samples
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
