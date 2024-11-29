#!/usr/bin/env python3

from ..tables_common import common_table_main, get_header, get_footer


def format_table(df):
    header = get_header(
        [
            r"Ensemble",
            r"$am_{\mathrm{v}}$",
            r"$af_{\mathrm{v}}$",
            r"$am_{\mathrm{av}}$",
            r"$af_{\mathrm{av}}$",
            r"$am_{\mathrm{t}}$",
            r"$am_{\mathrm{at}}$",
        ],
        hlines=1,
        column_separation="1em",
    )
    footer = get_footer()
    content = []
    previous_prefix = None
    for row in df.itertuples():
        if (next_prefix := row.ensemble_name[:4]) != previous_prefix:
            previous_prefix = next_prefix
            content.append("\\hline\n")

        content.append(
            (
                "{} & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ & "
                "${:.02uSL}$ & ${:.02uSL}$ \\\\\n"
            ).format(
                row.ensemble_name,
                row.v_mass,
                row.v_decay_constant,
                row.av_mass,
                row.av_decay_constant,
                row.t_mass,
                row.at_mass,
            )
        )

    return header + "".join(content) + footer


if __name__ == "__main__":
    common_table_main(format_table)
