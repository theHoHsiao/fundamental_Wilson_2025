#!/usr/bin/env python3


from ..tables_common import common_table_main, get_header, get_footer
from ..plots_common import ch_tag


def format_table(df):
    header = get_header(
        [
            "M",
            r"$\hat{f}_{M,\,\chi}^2$",
            r"$L_{\rm M}^f$",
            r"$W_M^f$",
            r"$\chi^2/{\rm N_{d.o.f}}$",
        ],
        column_separation="2em",
    )
    footer = get_footer()
    content = []

    for row in df.itertuples():
        content.append(
            (
                "${}$ & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02f}$ \\\\\n"
            ).format(
                "\\rm " + ch_tag(row.channel),
                row.F,
                row.L,
                row.W,
                row.chisquare,
            )
        )

    return header + "".join(content) + footer


if __name__ == "__main__":
    common_table_main(format_table, index_name="channel")
