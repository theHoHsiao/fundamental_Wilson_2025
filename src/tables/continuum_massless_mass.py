#!/usr/bin/env python3


from argparse import ArgumentParser, FileType

from ..dump import read_files_single_channel
from ..plots_common import ch_tag


def format_table(df):
    header = (
        "\\begin{tabular}{|c|c|c|c|c|} \n \\hline\\hline \n"
        "M & $\\hat{m}_{M,\\,\\chi}^2$ & $L_{\\rm M}^m$ & $W_M^m$ & $\\chi^2/{\\rm N_{d.o.f}}$  \\\\ \n"
        "\\hline \n"
    )
    footer = "\\hline \n \\end{tabular}"
    content = []

    for row in df.itertuples():
        content.append(
            (
                "${}$ & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02f}$ \\\\\n"
            ).format(
                "\\rm " + ch_tag(row.channel.split("_")[1]),
                row.M,
                row.L,
                row.W,
                row.chi_sqr_dof,
            )
        )

    return header + "".join(content) + footer


def get_args():
    parser = ArgumentParser(
        description="Take chipt extrapolation results with each beta value"
    )
    parser.add_argument(
        "data_filenames", nargs="+", help="Filenames of plaquette result files"
    )
    parser.add_argument("--output_file", type=FileType("w"), default="-")
    return parser.parse_args()


def main():
    args = get_args()
    data = read_files_single_channel(args.data_filenames)
    print(format_table(data), file=args.output_file)


if __name__ == "__main__":
    main()
