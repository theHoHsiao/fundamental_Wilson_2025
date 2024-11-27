#!/usr/bin/env python3


from ..tables_common import beta_table_main


def format_table(df):
    header = (
        "\\begin{tabular*}{\\columnwidth}{@{\\extracolsep{\\fill}}|c|c|c|c|c|} \n \\hline\\hline \n"
        "$\\beta$ & fit range ($am_0$) & $\\chi^2/N_{\\rm d.o.f.}$ & $A$ & $B$ \\\\ \\hline \n"
    )
    footer = "\\hline \n \\end{tabular*}"
    content = []

    for row in df.itertuples():
        content.append(
            ("{} & ${}$ & ${:.02f}$ & ${:.02uSL}$ & ${:.02uSL}$ \\\\\n").format(
                row.beta,
                row.bare_mass_range,
                row.chisquare,
                row.A,
                row.B,
            )
        )

    return header + "".join(content) + footer


if __name__ == "__main__":
    beta_table_main(format_table)
