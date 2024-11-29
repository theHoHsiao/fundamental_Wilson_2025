#!/usr/bin/env python3


import numpy as np

from ..tables_common import common_table_main, get_header, get_footer


def format_table(df):
    header = get_header(
        [
            r"$\beta$",
            r"fit range ($am_0$)",
            r"$\chi^2/N_{\rm d.o.f.}$",
            r"$A$",
            r"$B$",
        ],
        column_separation="1em",
    )
    footer = get_footer()
    content = []

    for row in df.itertuples():
        if np.isinf(row.chisquare):
            formatted_chisquare = r"$\infty$    "
        else:
            formatted_chisquare = f"{row.chisquare:.02f}"
        content.append(
            ("{} & ${}$ & {} & ${:.02uSL}$ & ${:.02uSL}$ \\\\\n").format(
                row.beta,
                row.bare_mass_range,
                formatted_chisquare,
                row.A,
                row.B,
            )
        )

    return header + "".join(content) + footer


if __name__ == "__main__":
    common_table_main(format_table, index_name="beta")
