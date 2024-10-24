#!/usr/bin/env python3

# from uncertainties import ufloat, UFloat

# from .bootstrap import BootstrapSampleSet

from ..tables_common import ensemble_table_main


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


if __name__ == "__main__":
    ensemble_table_main(format_table)
