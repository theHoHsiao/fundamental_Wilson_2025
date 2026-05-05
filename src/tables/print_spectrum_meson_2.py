#!/usr/bin/env python3

# from uncertainties import ufloat, UFloat

# from .bootstrap import BootstrapSampleSet

from ..tables_common import ensemble_table_main
import numpy as np


def format_table(df):
    header = (
        "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|}\n"
        "\\hline\\hline\n"
        r"Ensemble & $\beta$ & $m_0$ & $am_{\mathrm{V}}$ & $af_{\mathrm{V}} $ & $am_{\mathrm{T}} $ & $am_{\mathrm{AV}}$ "
        r" & $af_{\mathrm{AV}} $ & $am_{\mathrm{AT}}$  & $am_{\mathrm{S}}$ & $aE^{\mathrm{V}}_1$ \\"
        "\n\\hline"
    )
    footer = "\\hline\\hline\n\\end{tabular}"
    content = []
    previous_prefix = None
    for row in df.itertuples():
        if (next_prefix := row.ensemble_name[:4]) != previous_prefix:
            previous_prefix = next_prefix
            content.append("\\hline\n")

        if np.isnan(row.gevp_f_v_E1_mass.nominal_value):
            formatted_v_E1_mass = r"$\cdots$"
        else:
            formatted_v_E1_mass = "{:.02uSL}".format(row.gevp_f_v_E1_mass)


        content.append(
            (
                "{} & {} & {} & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ & "
                "${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ & {}\\\\\n"
            ).format(
                row.ensemble_name,
                row.beta,
                row.mF,
                row.gevp_f_v_E0_mass,
                row.f_v_decay_constant,
                row.gevp_f_t_E0_mass,
                row.gevp_f_av_E0_mass,
                row.f_av_decay_constant,
                row.gevp_f_at_E0_mass,
                row.gevp_f_s_E0_mass,
                formatted_v_E1_mass,
            )
        )

    return header + "".join(content) + footer


if __name__ == "__main__":
    ensemble_table_main(format_table)
