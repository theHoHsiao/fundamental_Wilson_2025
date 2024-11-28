#!/usr/bin/env python3


from ..tables_common import common_table_main


def format_table(df):
    header = (
        "\\begin{tabular}{|c|c|c|c|c|c|} \n \\hline \\hline \n"
        "~~~$\\beta$~~ & fit range ($am_0$) & $\\chi^2/N_{\\rm d.o.f.}$ & ~~~$\\tilde{C}$~~~ & ~~~$Y$~~~ & ~~~$y$~~~\\\\ \\hline \n"
    )
    footer = "\\hline \n \\end{tabular}"
    content = []

    for row in df.itertuples():
        content.append(
            (
                "{} & ${}$ & ${:.02f}$ & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ \\\\\n"
            ).format(
                row.beta,
                row.bare_mass_range,
                row.chisquare,
                row.A,
                row.B,
                row.y,
            )
        )

    return header + "".join(content) + footer


if __name__ == "__main__":
    common_table_main(format_table, index_name="beta")
