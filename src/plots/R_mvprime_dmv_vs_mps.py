#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat

from ..plots_common import (
    standard_plot_main,
    beta_iterator,
    add_figure_legend,
    ONE_COLUMN,
)


def plot(data, **kwargs):
    fig, ax = plt.subplots(
        1, 1, num="6", figsize=(ONE_COLUMN, 2.4), layout="constrained"
    )

    ax.set_xlim(-0.01, 1.51)
    ax.set_xlabel(r"$\hat{m}_{\rm ps}^2$")
    ax.set_ylabel(r"$m_{\rm v^\prime} / m_{\rm v}$")

    to_plot = []
    betas = sorted(set([datum["beta"] for datum in data]))
    for beta, colour, marker in beta_iterator(betas):
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
            alpha=1,
            color=colour,
            marker=marker,
            label=f"$\\beta={beta}$",
        )

    # Add QCD result
    # from https://pdg.lbl.gov/2019/listings/rpp2019-list-rho-770.pdf
    rho_770_mass = ufloat(775.26, 0.25)

    # from https://pdg.lbl.gov/2014/listings/rpp2014-list-rho-1450.pdf
    rho_1450_mass = ufloat(1465, 25)

    R_ratio_QCD = rho_1450_mass / rho_770_mass
    R_ratio_QCD_lower = R_ratio_QCD.nominal_value - R_ratio_QCD.std_dev
    R_ratio_QCD_upper = R_ratio_QCD.nominal_value + R_ratio_QCD.std_dev

    plt.fill_between(
        [0.01, 0.1],
        [R_ratio_QCD_lower, R_ratio_QCD_lower],
        [R_ratio_QCD_upper, R_ratio_QCD_upper],
        label=r"QCD $\rho(1450)/\rho(770)$",
        alpha=1,
        color="dodgerblue",
    )

    add_figure_legend(fig, ncol=3, title=None)
    return fig


if __name__ == "__main__":
    standard_plot_main(plot)
