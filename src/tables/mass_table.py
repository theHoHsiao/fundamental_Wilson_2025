#!/usr/bin/env python3

# from uncertainties import ufloat, UFloat

# from .bootstrap import BootstrapSampleSet

from argparse import ArgumentParser, FileType
from ..dump import read_files
from ..tables_common import by_ensemble_name


def get_args():
    parser = ArgumentParser()

    parser.add_argument(
        "data_filenames", nargs="+", help="Filenames of meson extraction result files"
    )
    parser.add_argument("--output_file", type=FileType("w"), default="-")
    return parser.parse_args()


def format_table(df):
    header = (
        "\\begin{tabular}{|c|c|c|c|c|c|c|}\n"
        "\\hline\\hline\n"
        r"Ensemble & $am_{\mathrm{PCAC}}$ & $am_{\mathrm{ps}}$ & $af_{\mathrm{ps}}$ & "
        r"$am_{\mathrm{s}}$ & $m_{\mathrm{ps}}L$ & $f_{\mathrm{ps}}L$ \\"
        "\n\\hline"
    )
    footer = "\\hline\\hline\n\\end{tabular}"
    content = []
    previous_prefix = None
    for row in df.itertuples():
        if (next_prefix := row.ensemble_name[:4]) != previous_prefix:
            previous_prefix = next_prefix
            content.append("\\hline\n")

        content.append(
            (
                "{} & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ & "
                "${:.02uSL}$ & ${:.02uSL}$ \\\\\n"
            ).format(
                row.ensemble_name,
                row.mPCAC,
                row.ps_mass,
                row.ps_matrix_element,  # TODO
                row.s_mass,
                row.ps_mass * row.Ns,
                row.ps_matrix_element * row.Ns,  # TODO
            )
        )

    return header + "".join(content) + footer


def main():
    args = get_args()
    data = read_files(args.data_filenames)
    print(
        format_table(data.sort_values(by="ensemble_name", key=by_ensemble_name)),
        file=args.output_file,
    )


if __name__ == "__main__":
    main()
