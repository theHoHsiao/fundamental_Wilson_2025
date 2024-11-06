#!/usr/bin/env python3


from ..tables_common import channel_table_main
from ..plots_common import ch_tag


def format_table(df):
    header = (
        "\\begin{tabular}{|c|c|c|c|c|} \n \\hline\\hline \n"
        "M & $\\hat{m}_{M,\\,\\chi}^2$ & $L_{\\rm M}^m$ & $W_M^m$ & $\\chi^2/{\\rm N_{d.o.f}}$  \\\\ \n"
        "\\hline \n"
    )
    footer = "\\hline \n \\end{tabular}"
    content = []

    for row in df.itertuples():
        content.append(
            (
                "${}$ & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02f}$ \\\\\n"
            ).format(
                "\\rm " + ch_tag(row.channel.split("_")[1]),
                row.M,
                row.L,
                row.W,
                row.chi_sqr_dof,
            )
        )

    return header + "".join(content) + footer


if __name__ == "__main__":
    channel_table_main(format_table)
