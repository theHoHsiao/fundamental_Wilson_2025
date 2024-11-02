#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

from ..plots_common import standard_plot_main
from ..fitting import meson_beta_quad


def plot_poly_M4(ax, A, B, C, ch, offset, color, x_i, x_f):
    n_fit = 1000
    Yfit = np.zeros(shape=(A.shape[0], n_fit))

    x = np.linspace(x_i, x_f, n_fit)

    y_up = np.zeros(n_fit)
    y_dn = np.zeros(n_fit)

    for n in range(A.shape[0]):
        Yfit[n] = A[n] + (B[n] * x**2) + (C[n] * x**4)

    for i in range(n_fit):
        y_err = Yfit[0:-1, i].std()
        y_up[i] = Yfit[-1, i] + y_err
        y_dn[i] = Yfit[-1, i] - y_err

    # ax.plot(x**2, Yfit[-1], "--", linewidth=0.75, alpha=0.6)
    ax.fill_between(
        x**2, y_up, y_dn, alpha=0.4, label=ch, facecolor=color, edgecolor=None
    )  # color=plt.gca().lines[-1].get_color()


def plot(data, external_data, fit_results):
    fig, ax = plt.subplots(2, 1, num="Figure_11a", figsize=(6, 8), layout="constrained")

    ax[0].set_ylim(0, 0.8)
    ax[0].set_xlim(0, 0.7)
    ax[0].set_xlabel(r"$(am_{\rm ps(PS)})^2$")
    ax[0].set_ylabel(r"$(am_{\rm v(V)})^2$")

    ax[1].set_ylim(0, 1.8)
    ax[1].set_xlim(0, 1.4)
    ax[1].set_xlabel(r"$\hat{m}_{\rm ps(PS)}^2$")
    ax[1].set_ylabel(r"$\hat{m}_{\rm v(V)}^2$")

    # plot fund. rep. results

    ax[0].errorbar(
        external_data.am_ps_sq_avg.values,
        external_data.am_v_sq_avg.values,
        xerr=external_data.am_ps_sq_err.values,
        yerr=external_data.am_v_sq_err.values,
        linestyle="",
        marker="o",
        markerfacecolor="none",
        elinewidth=1,
        capthick=1,
        capsize=1,  # zorder=5,
        color="r",
        alpha=0.8,
        label=r"$N_f=2$ (f) $Sp(4)$",
    )

    print(fit_results["m_V_vs_m_PS"])

    plot_poly_M4(
        ax[0],
        np.array(fit_results["m_V_vs_m_PS"]),
        np.array(fit_results["L"]),
        np.array(fit_results["W"]),
        "",
        0,
        "r",
        0,
        0.8366600265340756,
    )

    ax[1].errorbar(
        external_data.w0m_ps_sq_avg.values,
        external_data.w0m_v_sq_avg.values,
        xerr=external_data.w0m_ps_sq_err.values,
        yerr=external_data.w0m_v_sq_err.values,
        linestyle="",
        marker="o",
        markerfacecolor="none",
        elinewidth=1,
        capthick=1,
        capsize=1,  # zorder=5,
        color="r",
        alpha=0.8,
        label=r"$N_{\rm f}=2$ $Sp(4)$",
    )

    to_plot = []
    to_plot2 = []
    to_fit_x = []
    to_fit_y = []
    for datum in data:
        if "ps_mass_samples" not in datum:
            continue
        if "v_mass_samples" not in datum:
            continue

        X = datum["ps_mass_samples"].samples ** 2
        Y = datum["v_mass_samples"].samples ** 2

        w0 = datum["w0_samples"].samples ** 2

        X_w0 = X * w0
        Y_w0 = Y * w0

        to_plot.append((Y.mean(), Y.std(), X.mean(), X.std()))
        to_plot2.append((Y_w0.mean(), Y_w0.std(), X_w0.mean(), X_w0.std()))

        to_fit_x.append(np.append(X, datum["ps_mass_samples"].mean ** 2))
        to_fit_y.append(np.append(Y, datum["v_mass_samples"].mean ** 2))

    y_values, y_errors, x_values, x_errors = zip(*to_plot)
    y2_values, y2_errors, x2_values, x2_errors = zip(*to_plot2)

    # print(np.array(to_fit_x).shape)
    fit_val, X2 = meson_beta_quad(np.array(to_fit_x), np.array(to_fit_y))
    plot_poly_M4(
        ax[0], fit_val[0], fit_val[1], fit_val[2], "", 0, "b", 0, 0.8366600265340756
    )

    ax[0].errorbar(
        x_values,
        y_values,
        xerr=x_errors,
        yerr=y_errors,
        ls="none",
        alpha=0.7,
        color="b",
        marker="s",
    )

    ax[1].errorbar(
        x2_values,
        y2_values,
        xerr=x2_errors,
        yerr=y2_errors,
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
