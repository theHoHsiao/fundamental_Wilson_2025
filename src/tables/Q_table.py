#!/usr/bin/env python3

import numpy as np

from ..definitions_common import format_definitions
from ..tables_common import common_table_main, get_header, get_footer


def format_list(ensembles):
    if len(ensembles) == 0:
        return ""
    elif len(ensembles) == 1:
        return ensembles[0]
    else:
        return ", ".join(ensembles[:-1]) + ", and " + ensembles[-1]


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
    incomplete_ensembles = []
    previous_prefix = None
    for row in df.itertuples():
        if (next_prefix := row.ensemble_name[:4]) != previous_prefix:
            previous_prefix = next_prefix
            content.append("\\hline\n")

        if np.isnan(row.w0.nominal_value) or np.isnan(row.w0.std_dev):
            w0 = r"\cdots"
            Ncfg_GF = r"$\cdots$"
            incomplete_ensembles.append(row.ensemble_name)
        else:
            w0 = f"{row.w0:.02uSL}"
            Ncfg_GF = row.Ncfg_GF

        content.append(
            "{} & {} & ${}$ & ${:.02uSL}$ & ${:.02uSL}$ \\\\\n".format(
                row.ensemble_name,
                Ncfg_GF,
                w0,
                row.Q0,
                row.sigma_Q,
            )
        )

    formatted_ensembles = format_list(incomplete_ensembles)
    definitions = format_definitions({"WZeroIncompleteEnsembles": formatted_ensembles})
    return header + "".join(content) + footer, definitions


if __name__ == "__main__":
    common_table_main(format_table, definitions=True)
