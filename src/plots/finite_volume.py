#!/usr/bin/env python3

import matplotlib.pyplot as plt
from ..plots_common import standard_plot_main, ch_tag
from argparse import ArgumentParser


def get_args():
    parser = ArgumentParser()

    parser.add_argument(
        "data_filenames",
        nargs="+",
        metavar="sample_filename",
        help="Filenames of sample files containing data to plot",
    )
    parser.add_argument(
        "--plot_file",
        default=None,
        help="Where to place the resulting plot. Default is to output to screen.",
    )
    parser.add_argument(
        "--plot_file2",
        default=None,
        help="Where to place the resulting plot. Default is to output to screen.",
    )
    parser.add_argument(
        "--plot_styles",
        default="styles/paperdraft.mplstyle",
        help="Stylesheet to use for plots",
    )
    return parser.parse_args()


def plot(data):
    fig, axs = plt.subplots(2, 1, num="Figure_2a", figsize=(6, 8), layout="constrained")

    ch_i = 0

    for ch in ["ps", "v"]:
        ax = axs[ch_i]

        ax.set_xlabel(r"$m_{\rm ps}^{\rm inf} L$")
        ax.set_ylabel(r"$am_{\mathrm{" + ch_tag(ch) + "}}$")

        betas = sorted(set([datum["beta"] for datum in data]))
        markers = "s"
        for beta_idx, (beta, marker) in enumerate(zip(betas, markers)):
            to_plot = []
            for datum in data:
                # print(datum)
                if datum["Ns"] == 24 and datum["ensemble_name"] == "ASB4M3":
                    m_ps_inf = datum["ps_mass_samples"]
                if datum["beta"] != beta:
                    continue

                if "ps_mass_samples" not in datum:
                    continue

                X = m_ps_inf.samples * datum["Ns"]
                Y = datum[f"{ch}_mass_samples"].samples

                to_plot.append((Y.mean(), Y.std(), X.mean(), X.std()))

            y_values, y_errors, x_values, x_errors = zip(*to_plot)

            ax.errorbar(
                x_values,
                y_values,
                xerr=x_errors,
                yerr=y_errors,
                ls="none",
                alpha=0.7,
                color="b",
                marker="s",
                label=f"{beta}",
            )

        # ax.set_xlim(2, 11)
        ch_i += 1

    return fig


"""
def main():
    args = get_args()
    plt.style.use(args.plot_styles)
    data = read_sample_files(args.data_filenames)
    fig = plot(data)
    save_or_show(fig, args.plot_file)
"""


if __name__ == "__main__":
    standard_plot_main(plot)
