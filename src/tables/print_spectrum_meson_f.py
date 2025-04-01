#!/usr/bin/env python3

# from uncertainties import ufloat, UFloat

# from .bootstrap import BootstrapSampleSet

from ..tables_common import ensemble_table_main
import numpy as np


def format_table(df):
    header = (
        "\\begin{tabular}{|c|c|c|c|c|c|c|c|}\n"
        "\\hline\\hline\n"
        r"Ensemble & $w_0$ & $aE^0_{\mathrm{ps}}$ & $am_{\mathrm{ps}} $ & $af_{\mathrm{ps}} $ &"
        r" $aE^0_{\mathrm{v}}$ & $am_{\mathrm{v}} $ & $af_{\mathrm{v}} $ \\"
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
                "{} & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ & "
                "${:.02uSL}$ & ${:.02uSL}$\\\\\n"
            ).format(
                row.ensemble_name,
                row.w0,
                row.gevp_f_ps_E0_mass,
                row.f_ps_mass,
                row.f_ps_decay_constant,
                row.gevp_f_v_E0_mass,
                row.f_v_mass,
                row.f_v_decay_constant,
            )
        )

    return header + "".join(content) + footer


if __name__ == "__main__":
    ensemble_table_main(format_table)
