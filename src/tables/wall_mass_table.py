#!/usr/bin/env python3


from ..tables_common import common_table_main, get_header, get_footer


def format_table(df):
    header = get_header(
        [
            r"Ensemble",
            r"$am_{\mathrm{PCAC}}$",
            r"$am_{\mathrm{ps}}$",
            r"$af_{\mathrm{ps}}$",
            r"$am_{\mathrm{s}}$",
            r"$m_{\mathrm{ps}}L$",
            r"$f_{\mathrm{ps}}L$",
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

        content.append(
            (
                "{} & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ & "
                "${:.02uSL}$ & ${:.02uSL}$ \\\\\n"
            ).format(
                row.ensemble_name,
                row.mPCAC,
                row.ps_mass,
                row.ps_decay_constant,
                row.s_mass,
                row.ps_mass * row.Ns,
                row.ps_decay_constant * row.Ns,
            )
        )

    return header + "".join(content) + footer


if __name__ == "__main__":
    common_table_main(format_table)
