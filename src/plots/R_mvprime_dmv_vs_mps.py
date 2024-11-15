#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from ..plots_common import standard_plot_main, beta_color
from ..bootstrap import BOOTSTRAP_SAMPLE_COUNT


def plot(data, **kwargs):
    fig, ax = plt.subplots(1, 1, num="6", figsize=(3.5, 2.4), layout="constrained")

    # from https://pdg.lbl.gov/2019/listings/rpp2019-list-rho-770.pdf
    rho_770_mass = np.random.normal(775.26, 0.25, BOOTSTRAP_SAMPLE_COUNT)

    # from https://pdg.lbl.gov/2014/listings/rpp2014-list-rho-1450.pdf
    rho_1450_mass = np.random.normal(1465, 25, BOOTSTRAP_SAMPLE_COUNT)

    R_ratio_QCD = rho_1450_mass / rho_770_mass

    plt.fill_between(
        [0.01, 0.1],
        [
            R_ratio_QCD.mean() - R_ratio_QCD.std(),
            R_ratio_QCD.mean() - R_ratio_QCD.std(),
        ],
        [
            R_ratio_QCD.mean() + R_ratio_QCD.std(),
            R_ratio_QCD.mean() + R_ratio_QCD.std(),
        ],
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

            w0 = datum["w0_samples"]

            X = (datum["smear_ps_mass_samples"] * w0) ** 2

            Y = datum["smear_rhoE1_mass_samples"] / datum["smear_v_mass_samples"]

            to_plot.append((Y.mean, Y.samples.std(), X.mean, X.samples.std()))

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
