#!/usr/bin/env python3

import numpy as np

from ..tables_common import common_table_main


def format_table(df):
    header = (
        "\\begin{tabular}{|c|c|c|c|c|c|c|c|}\n"
        "\\hline\\hline\n"
        r"Ensemble & $Q_0$ & $\sigma_Q$ & $\tau_{\mathrm{exp}}^Q$ & "
        r"$N_{\mathrm{traj}}^{\mathrm{GF}}$ & $\delta_{\mathrm{traj}}^{\mathrm{GF}}$ & "
        r"$w_0 / a$ & $\tau_{\mathrm{exp}}^{w_0}$ \\\\\n"
        r"\hline"
    )
    footer = "\\hline\\hline\n\\end{tabular}"
    content = []
    previous_prefix = None
    for row in df.itertuples():
        if (next_prefix := row.ensemble_name[:4]) != previous_prefix:
            previous_prefix = next_prefix
            content.append("\\hline\n")

        if np.isnan(row.w0.nominal_value) or np.isnan(row.w0.std_dev):
            w0 = r"\cdots"
            tau_exp_w0 = r"\cdots"
            num_configs = r"$\cdots$"
            trajectory_step = r"$\cdots$"
        else:
            w0 = f"{row.w0:.02uSL}"
            tau_exp_w0 = f"{row.tau_exp_w0:.02uSL}"
            num_configs = row.num_configs
            trajectory_step = row.trajectory_step

        content.append(
            (
                "{} & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ "
                "& {} & {} & ${}$ & ${}$ \\\\\n"
            ).format(
                row.ensemble_name,
                row.Q0,
                row.sigma_Q,
                row.tau_exp_Q,
                num_configs,
                trajectory_step,
                w0,
                tau_exp_w0,
            )
        )

    return header + "".join(content) + footer


if __name__ == "__main__":
    common_table_main(format_table)
