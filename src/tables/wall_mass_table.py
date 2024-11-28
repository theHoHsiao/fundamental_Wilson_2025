#!/usr/bin/env python3


from ..tables_common import common_table_main


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
                row.ps_decay_constant,
                row.s_mass,
                row.ps_mass * row.Ns,
                row.ps_decay_constant * row.Ns,
            )
        )

    return header + "".join(content) + footer


if __name__ == "__main__":
    common_table_main(format_table)
