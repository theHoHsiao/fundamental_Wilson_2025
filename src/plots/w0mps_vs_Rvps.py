#!/usr/bin/env python3

import matplotlib.pyplot as plt
from ..dump import read_sample_files
from ..plots_common import (
    save_or_show,
    beta_iterator,
    ch_tag,
    channel_color,
    add_figure_legend,
    ONE_COLUMN,
    MPS2_right_end,
)
from argparse import ArgumentParser
import numpy as np


def get_args():
    parser = ArgumentParser()

    parser.add_argument(
        "data_filenames",
        nargs="+",
        metavar="sample_filename",
        help="Filenames of sample files containing data to plot",
    )
    parser.add_argument(
        "--fit_parameters",
        nargs="+",
        metavar="sample_filename",
        help="Filenames of sample files containing data to plot",
    )
    parser.add_argument(
        "--plot_file_data",
        default=None,
        help=(
            "Where to place the resulting plot showing data. "
            "Default is to output to screen."
        ),
    )
    parser.add_argument(
        "--plot_file_summary",
        default=None,
        help=(
            "Where to place the resulting summary plot (fit results only). "
            "Default is to output to screen."
        ),
    )
    parser.add_argument(
        "--plot_styles",
        default="styles/paperdraft.mplstyle",
        help="Stylesheet to use for plots",
    )
    return parser.parse_args()



def plot(data):
    data_fig, data_axes = plt.subplots(
        1, 1, num="Figure_12", figsize=(ONE_COLUMN, 3.5), layout="constrained"
    )

    for ax, ch in zip([data_axes], ["v"]):

        ax.set_xlabel(r"$\hat{m}_{\mathrm{PS}}^2$")
        ax.set_ylabel(r"$  m_{\mathrm{" + ch_tag(ch) + r"}} / m_{\mathrm{PS}} $")
        ax.set_xlim(0, MPS2_right_end)

        betas = sorted(set([datum["beta"] for datum in data]))
        for beta, colour, marker in beta_iterator(betas):
            to_plot = []
            for datum in data:
                
                if datum["beta"] != beta:
                    continue
                
                if "w0_samples" not in datum:
                    continue
                
                if "gevp_f_ps_E0_mass_samples" not in datum:
                    continue


                X = ( datum["w0_samples"] * datum["gevp_f_ps_E0_mass_samples"]) ** 2
                Y = (  datum["gevp_f_v_E0_mass_samples"] / datum[f"gevp_f_ps_E0_mass_samples"]) 

                to_plot.append((Y.mean, Y.samples.std(), X.mean, X.samples.std()))

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
                label=f"{beta}",
            )

    add_figure_legend(data_fig)

    return data_fig


def main():
    args = get_args()
    plt.style.use(args.plot_styles)
    data = read_sample_files(args.data_filenames)
    data_fig = plot(data)
    save_or_show(data_fig, args.plot_file_data)


if __name__ == "__main__":
    main()
