#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

from ..plots_common import standard_plot_main


def plot_poly_M4(ax, A, B, C, ch, color):
    n_fit = 1000
    Yfit = np.zeros(shape=(A.shape[0], n_fit))

    x_min, x_max = ax.get_xlim()
    x_i = np.sqrt(x_min)
    x_f = np.sqrt(x_max)
    x = np.linspace(x_i, x_f, n_fit)

    y_up = np.zeros(n_fit)
    y_dn = np.zeros(n_fit)

    for n in range(A.shape[0]):
        Yfit[n] = A[n] + (B[n] * x**2) + (C[n] * x**4)

    for i in range(n_fit):
        y_err = Yfit[0:-1, i].std()
        y_up[i] = Yfit[-1, i] + y_err
        y_dn[i] = Yfit[-1, i] - y_err

    ax.fill_between(
        x**2, y_up, y_dn, alpha=0.4, label=ch, facecolor=color, edgecolor=None
    )


def plot(data, external_data, fit_results):
    fig, ax = plt.subplots(num="Figure_11", figsize=(3.5, 2.4), layout="constrained")

    ax.set_ylim(0, 1.8)
    ax.set_xlim(0, 1.4)
    ax.set_xlabel(r"$\hat{m}^2_{\rm ps},\, \hat{m}^2_{\rm PS}$")
    ax.set_ylabel(r"$\hat{m}^2_{\rm v},\, \hat{m}^2_{\rm V} $")

    # plot fund. rep. results

    ax.errorbar(
        external_data.w0m_ps_sq_avg.values,
        external_data.w0m_v_sq_avg.values,
        xerr=external_data.w0m_ps_sq_err.values,
        yerr=external_data.w0m_v_sq_err.values,
        linestyle="",
        marker="o",
        markerfacecolor="none",
        elinewidth=1,
        capthick=1,
        capsize=1,
        color="C3",
        alpha=1,
        label=r"$N_{\rm f}=2$ $Sp(4)$",
    )

    to_plot = []
    for datum in data:
        if "ps_mass_samples" not in datum:
            continue
        if "v_mass_samples" not in datum:
            continue

        X = datum["ps_mass_samples"] ** 2
        Y = datum["v_mass_samples"] ** 2

        w0 = datum["w0_samples"] ** 2

        X_w0 = X * w0
        Y_w0 = Y * w0

        to_plot.append((Y_w0.mean, Y_w0.samples.std(), X_w0.mean, X_w0.samples.std()))

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
