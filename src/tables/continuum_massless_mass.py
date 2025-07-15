#!/usr/bin/env python3


from ..tables_common import common_table_main, get_header, get_footer
from ..plots_common import ch_tag


def format_table(df):
    header = get_header(
        [
            r"M",
            r"$\hat{m}_{M,\,\chi}^2$",
            r"$L_{\rm M}^m$",
            r"$W_M^m$",
            r"$\chi^2/{\rm N_{d.o.f}}$",
        ]
    )
    footer = get_footer()
    content = []

    for row in df.itertuples():
        content.append(
            (
                "${}$ & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02f}$ \\\\\n"
            ).format(
                "\\rm " + ch_tag(row.channel.upper()),
                row.M,
                row.L,
                row.W,
                row.chisquare,
            )
        )

    return header + "".join(content) + footer


if __name__ == "__main__":
    common_table_main(format_table, index_name="channel")