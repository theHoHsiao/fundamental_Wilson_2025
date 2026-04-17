#!/usr/bin/env python3

# from uncertainties import ufloat, UFloat

# from .bootstrap import BootstrapSampleSet

from ..tables_common import ensemble_table_main
#import numpy as np


def format_table(df):
    header = (
        "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|}\n"
        "\\hline\\hline\n"
        r"Ensemble & $\beta$ & $m_0$ & $\varepsilon_{\rm W}$ & $N^{\mathrm{diff}}_{\rm W}$ & $N^{\mathrm{max}}_{\rm W}$ & "
        r" $am_{\mathrm{PCAC}}$ & $am_{\mathrm{PS}}$  & $af_{\mathrm{PS}} $ &"
        r" $m_{\mathrm{PS}}L$ & $f_{\mathrm{PS}}L$\\"
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
                "{} & {} & {} & {} & {} & {} & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ "
                "\\\\\n"
            ).format(
                row.ensemble_name,
                row.beta,
                row.mF,
                row.epsilon,
                row.n_smear_diff,
                row.n_smear_max,
                row.mPCAC,
                row.gevp_f_ps_E0_mass,
                row.f_ps_decay_constant,
                row.gevp_f_ps_E0_mass * row.Ns,
                row.f_ps_decay_constant * row.Ns,
            )
        )

    return header + "".join(content) + footer


if __name__ == "__main__":
    ensemble_table_main(format_table)
