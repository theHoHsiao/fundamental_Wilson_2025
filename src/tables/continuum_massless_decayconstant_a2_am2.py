#!/usr/bin/env python3


from ..tables_common import ansatze_table_main, get_header, get_footer
from ..plots_common import ch_tag


def format_table(df, ansatz):
    header = get_header(
        [
            r"M",
            r"$\hat{f}_{M,\,\chi}^2$",
            r"$L_{\rm M}^f$",
            r"$W_{\rm M}^f$",
            r"$R_{\rm M}^f$",
            r"$C_{\rm M}^f$",
            r"$\chi^2/{\rm N_{d.o.f}}$",
        ]
    )
    footer = get_footer()
    content = []

    for row in df.itertuples():
        # Helper to get attribute or return "-"
        def get_val(obj, attr, fmt="{:.02uSL}"):
            val = getattr(obj, attr, None)
            if val is None:
                return "-"
            # Return formatted string if it exists
            return fmt.format(val)

        content.append(
            (
                "${}$ & ${}$ & ${}$ & ${}$ & ${}$ & ${}$ & ${}$ \\\\\n"
            ).format(
                "\\rm " + ch_tag(row.channel),
                get_val(row, f"F_{ansatz}"),
                get_val(row, f"Lf_{ansatz}"),
                get_val(row, f"Wf_{ansatz}"),
                get_val(row, f"Rf_{ansatz}"),
                get_val(row, f"Cf_{ansatz}"),
                get_val(row, f"chisquare", fmt="{:.02f}"), # Custom format for chisq
            )
        )

    return header + "".join(content) + footer


if __name__ == "__main__":
    ansatze_table_main(format_table, index_name="channel")