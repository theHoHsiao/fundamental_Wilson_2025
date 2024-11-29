#!/usr/bin/env python3

from ..tables_common import common_table_main, get_header, get_footer
import numpy as np


def format_table(df):
    header = get_header(
        [
            r"Ensemble",
            r"$am_{\mathrm{ps}}$",
            r"$am_{\mathrm{s}}$",
            r"$m_{\mathrm{v}}$",
            r"$am_{\mathrm{t}}$",
            r"$m_{\mathrm{av}}$",
            r"$m_{\mathrm{at}}$",
            r"$m_{\mathrm{v}^\prime}$",
        ],
        column_separation="0.5em",
        hlines=1,
    )
    footer = get_footer()
    content = []
    previous_prefix = None
    for row in df.itertuples():
        if (next_prefix := row.ensemble_name[:4]) != previous_prefix:
            previous_prefix = next_prefix
            content.append("\\hline\n")

        if np.isnan(row.smear_rhoE1_mass.nominal_value) or np.isnan(
            row.smear_rhoE1_mass.std_dev
        ):
            smear_rhoE1_mass = r"\cdots"

        else:
            smear_rhoE1_mass = f"{row.smear_rhoE1_mass:.02uSL}"

        content.append(
            (
                "{} & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ & "
                "${:.02uSL}$ & ${:.02uSL}$ & ${}$ \\\\\n"
            ).format(
                row.ensemble_name,
                row.smear_ps_mass,
                row.smear_s_mass,
                row.smear_v_mass,
                row.smear_t_mass,
                row.smear_av_mass,
                row.smear_at_mass,
                smear_rhoE1_mass,
            )
        )

    return header + "".join(content) + footer


if __name__ == "__main__":
    common_table_main(format_table)
