#!/usr/bin/env python3

import matplotlib.pyplot as plt
from ..dump import read_sample_files
from ..plots_common import save_or_show, ch_tag
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
        "--plot_styles",
        default="styles/paperdraft.mplstyle",
        help="Stylesheet to use for plots",
    )
    return parser.parse_args()


def plot(data):
    fig, ax = plt.subplots(1, 1, num="Figure_6", figsize=(6, 4), layout="constrained")

    ch = "ps"

    ax.set_xlabel(r"$m_{\rm ps}^{\rm inf} L$")
    ax.set_ylabel(r"$am_{\mathrm{" + ch_tag(ch) + "}}$")

    betas = sorted(set([datum["beta"] for datum in data]))
    markers = "s"
    for beta_idx, (beta, marker) in enumerate(zip(betas, markers)):
        to_plot = []
        for datum in data:
            if datum["beta"] != beta:
                continue

            if "ps_mass_samples" not in datum:
                continue

            X = datum["mpcac_samples"].samples
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

    return fig


def main():
    args = get_args()
    plt.style.use(args.plot_styles)
    data = read_sample_files(args.data_filenames)
    fig = plot(data)
    save_or_show(fig, args.plot_file)


if __name__ == "__main__":
    main()
