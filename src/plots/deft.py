#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

from ..plots_common import standard_plot_main, beta_color


def plot_axpb_y_minus(ax, A, B, ch, offset, color, x_i, x_f):
    n_fit = 1000
    Yfit = np.zeros(shape=(A.shape[0], n_fit))

    x = np.linspace(x_i, x_f, n_fit)

    y_up = np.zeros(n_fit)
    y_dn = np.zeros(n_fit)

    for n in range(A.shape[0]):
        Yfit[n] = A[n] + (B[n] * -(x**2))

    for i in range(n_fit):
        y_err = Yfit[0:-1, i].std()
        y_up[i] = Yfit[-1, i] + y_err
        y_dn[i] = Yfit[-1, i] - y_err

    # ax.plot(x**2, Yfit[-1], "--", linewidth=0.75, alpha=0.6)
    ax.fill_between(
        -(x**2), y_up, y_dn, alpha=0.4, label=ch, facecolor=color, edgecolor=None
    )  # color=plt.gca().lines[-1].get_color()


def plot(data, fit_results, **kwargs):
    fig, ax = plt.subplots(
        1, 1, num="Figure_22", figsize=(3.5, 2.4), layout="constrained"
    )

    ax.set_ylim(0.4, 1.7)
    ax.set_xlim(-7, -4)
    ax.set_xlabel(r"$\log  [ (af_{\rm ps})^2  ] $")
    ax.set_ylabel(r"$\log [(am_{\rm ps})^2 / (am_{\rm PCAC})]$")

    betas = sorted(set([datum["beta"] for datum in data]))
    markers = "o^vsx+"

    for beta_idx, (beta, marker) in enumerate(zip(betas, markers)):
        to_plot = []
        for datum in data:
            if datum["beta"] != beta:
                continue
            if "ps_mass_samples" not in datum:
                continue

            Y = (datum["ps_mass_samples"]) ** 2 / datum["mPCAC_samples"]

            X = (datum["ps_decay_constant_samples"]) ** 2

            X = np.log(X)
            Y = np.log(Y)

            to_plot.append((Y.mean, Y.samples.std(), X.mean, X.samples.std()))

        if len(to_plot) < 3:
            continue

        for tmp_result in fit_results:
            if tmp_result["beta"] == str(beta):
                fit_result = tmp_result

        # print(fit_results)
        plot_axpb_y_minus(
            ax,
            np.random.normal(fit_result[f"A_b{beta}_samples"].mean, 0.0005, 200),
            np.random.normal(
                fit_result[f"B_b{beta}_samples"].mean, 0.0005, 200
            ),  # line with a small band
            "",
            0,
            beta_color(beta),
            2.645,
            2,
        )

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
        ncol=5,
        borderaxespad=0.2,
    )

    return fig


if __name__ == "__main__":
    standard_plot_main(plot)
