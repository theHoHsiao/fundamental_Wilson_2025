#!/usr/bin/env python3

from ..tables_common import common_table_main
import numpy as np


def format_table(df):
    header = (
        "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|}\n"
        "\\hline\\hline\n"
        r"Ensemble & ~~~$m_{\mathrm{ps}} / f_{\mathrm{ps}}$~~~ & ~~~$m_{\mathrm{s}}/ f_{\mathrm{ps}}$~~~ & "
        r"~~~$m_{\mathrm{v}}/ f_{\mathrm{ps}}$~~~ & ~~~$m_{\mathrm{t}}/ f_{\mathrm{ps}}$~~~ & ~~~$m_{\mathrm{av}}/ f_{\mathrm{ps}}$~~~ &"
        r" ~~~$m_{\mathrm{at}}/ f_{\mathrm{ps}}$~~~ & ~~~$m_{\mathrm{v}^\prime} / f_{\mathrm{ps}} $~~~ & ~~~$m_{\rm v^\prime} / m_{\rm v}$~~~ \\"
        "\n\\hline"
    )
    footer = "\\hline\\hline\n\\end{tabular}"
    content = []
    previous_prefix = None
    for row in df.itertuples():
        if (next_prefix := row.ensemble_name[:4]) != previous_prefix:
            previous_prefix = next_prefix
            content.append("\\hline\n")

        if np.isnan(row.smear_rhoE1_Rfps.nominal_value) or np.isnan(
            row.smear_rhoE1_Rfps.std_dev
        ):
            smear_rhoE1_Rfps = r"\cdots"
            smear_rhoE1_Rmv = r"\cdots"
        else:
            smear_rhoE1_Rfps = f"{row.smear_rhoE1_Rfps:.02uSL}"
            smear_rhoE1_Rmv = f"{row.smear_rhoE1_Rmv:.02uSL}"

        content.append(
            (
                "{} & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ & "
                "${:.02uSL}$ & ${:.02uSL}$ & ${}$ & ${}$ \\\\\n"
            ).format(
                row.ensemble_name,
                row.smear_ps_Rfps,
                row.smear_s_Rfps,
                row.smear_v_Rfps,
                row.smear_t_Rfps,
                row.smear_av_Rfps,
                row.smear_at_Rfps,
                smear_rhoE1_Rfps,
                smear_rhoE1_Rmv,
            )
        )

    return header + "".join(content) + footer


if __name__ == "__main__":
    common_table_main(format_table)
