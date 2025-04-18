#!/usr/bin/env python3

# from uncertainties import ufloat, UFloat

# from .bootstrap import BootstrapSampleSet

from ..tables_common import ensemble_table_main
#import numpy as np


def format_table(df):
    header = (
        "\\begin{tabular}{|c|c|c|c|c|c|c|}\n"
        "\\hline\\hline\n"
        r"Ensemble & $\beta$ & $m_0$ & $Nt$ & $Ns$ & $w_0$ & $<P>$ \\"
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
                "{} & {} & {} & {} & {:-1} & ${:.02uSL}$ & ${:.02uSL}$ "
                "\\\\\n"
            ).format(
                row.ensemble_name,
                row.beta,
                row.mF,
                int(row.Nt),
                int(row.Ns),
                row.w0,
                row.avg_plaquette,
            )
        )

    return header + "".join(content) + footer


if __name__ == "__main__":
    ensemble_table_main(format_table)
