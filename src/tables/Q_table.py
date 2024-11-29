#!/usr/bin/env python3

import numpy as np

from ..tables_common import common_table_main, get_header, get_footer


def format_table(df):
    header = get_header(
        [
            r"Ensemble",
            r"$N_{\mathrm{traj}}^{\mathrm{GF}}$",
            r"$w_0 / a$",
            r"$Q_0$",
            r"$\sigma_Q$ ",
        ],
        hlines=1,
    )
    footer = get_footer()
    content = []
    previous_prefix = None
    for row in df.itertuples():
        if (next_prefix := row.ensemble_name[:4]) != previous_prefix:
            previous_prefix = next_prefix
            content.append("\\hline\n")

        if np.isnan(row.w0.nominal_value) or np.isnan(row.w0.std_dev):
            w0 = r"\cdots"
            num_configs = r"$\cdots$"
        else:
            w0 = f"{row.w0:.02uSL}"
            num_configs = row.num_configs

        content.append(
            "{} & {} & ${}$ & ${:.02uSL}$ & ${:.02uSL}$ \\\\\n".format(
                row.ensemble_name,
                num_configs,
                w0,
                row.Q0,
                row.sigma_Q,
            )
        )

    return header + "".join(content) + footer


if __name__ == "__main__":
    common_table_main(format_table)
