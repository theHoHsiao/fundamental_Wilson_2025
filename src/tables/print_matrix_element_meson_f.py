#!/usr/bin/env python3

# from uncertainties import ufloat, UFloat

# from .bootstrap import BootstrapSampleSet

from ..tables_common import ensemble_table_main
import numpy as np


def format_table(df):
    header = (
        "\\begin{tabular}{|c|c|c|c|c|c|c|}\n"
        "\\hline\\hline\n"
        r"Ensemble & $am_{\mathrm{ps}} $ & $af_{\mathrm{ps}} $ & $\chi^2$ / {\textrm{d.o.f.}} &"
        r" $am_{\mathrm{v}} $ & $af_{\mathrm{v}} $ & $\chi^2$ / {\textrm{d.o.f.}} \\"
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
                "{} & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02}$ & ${:.02uSL}$ & "
                "${:.02uSL}$ & ${:.02}$\\\\\n"
            ).format(
                row.ensemble_name,
                row.f_ps_mass,
                row.f_ps_decay_constant,
                row.f_ps_chisquare,
                row.f_v_mass,
                row.f_v_decay_constant,
                row.f_v_chisquare,
            )
        )

    return header + "".join(content) + footer


if __name__ == "__main__":
    ensemble_table_main(format_table)
