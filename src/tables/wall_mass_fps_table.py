#!/usr/bin/env python3

# from uncertainties import ufloat, UFloat

# from .bootstrap import BootstrapSampleSet

from ..tables_common import ensemble_table_main


def format_table(df):
    header = (
        "\\begin{tabular}{|c|c|c|c|c|c|c|}\n"
        "\\hline\\hline\n"
        r"Ensemble & ~~~$m_{\mathrm{ps}} / f_{\mathrm{ps}}$~~~ & ~~~$m_{\mathrm{s}}/ f_{\mathrm{ps}}$~~~ & "
        r"~~~$m_{\mathrm{v}}/ f_{\mathrm{ps}}$~~~ & ~~~$m_{\mathrm{t}}/ f_{\mathrm{ps}}$~~~ & ~~~$m_{\mathrm{av}}/ f_{\mathrm{ps}}$~~~ &"
        r" ~~~$m_{\mathrm{at}}/ f_{\mathrm{ps}}$~~~\\"
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
                row.ps_Rfps,
                row.s_Rfps,
                row.v_Rfps,
                row.t_Rfps,
                row.av_Rfps,
                row.at_Rfps,
            )
        )

    return header + "".join(content) + footer


if __name__ == "__main__":
    ensemble_table_main(format_table)
