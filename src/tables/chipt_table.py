#!/usr/bin/env python3


from argparse import ArgumentParser, FileType

from ..dump import read_files_single_beta


def format_table(df):
    header = (
        "\\begin{tabular}{|c|c|c|c|c|} \n \\hline\\hline \n"
        "~~$\\beta$~~ & fit range ($am_0$) & $\\chi^2/N_{\\rm d.o.f.}$ & ~~~$A$~~~ & ~~~$B$~~~ \\\\ \\hline \n"
    )
    footer = "\\hline \n \\end{tabular}"
    content = []

    for row in df.itertuples():
        content.append(
            ("{} & ${}$ & ${:.02f}$ & ${:.02uSL}$ & ${:.02uSL}$ \\\\\n").format(
                row.beta,
                row.bare_mass_range,
                row.chi_sqr_dof,
                row.A,
                row.B,
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
    data = read_files_single_beta(args.data_filenames)
    print(format_table(data), file=args.output_file)


if __name__ == "__main__":
    main()
