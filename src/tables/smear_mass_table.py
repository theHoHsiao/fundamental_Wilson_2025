#!/usr/bin/env python3

# from uncertainties import ufloat, UFloat

# from .bootstrap import BootstrapSampleSet

from ..tables_common import ensemble_table_main


def format_table(df):
    header = (
        "\\begin{tabular}{|c|c|c|c|c|c|c|c}\n"
        "\\hline\\hline\n"
        r"Ensemble & ~~~~~~$am_{\mathrm{ps}}$~~~~~~ & ~~~~~~$am_{\mathrm{s}}$~~~~~~ & ~~~~~~$m_{\mathrm{v}}$~~~~~~ & "
        r"~~~~~~$am_{\mathrm{t}}$~~~~~~ & ~~~~~~$m_{\mathrm{av}}$~~~~~~ & ~~~~~~$m_{\mathrm{at}}$~~~~~~ & ~~~~~~$m_{\mathrm{v}^\prime}$~~~~~~ \\"
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
                "${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ \\\\\n"
            ).format(
                row.ensemble_name,
                row.smear_ps_mass,
                row.smear_s_mass,
                row.smear_v_mass,
                row.smear_t_mass,
                row.smear_av_mass,
                row.smear_at_mass,
                row.smear_rhoE1_mass,
            )
        )

    return header + "".join(content) + footer


if __name__ == "__main__":
    ensemble_table_main(format_table)
