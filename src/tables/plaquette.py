#!/usr/bin/env python3


from ..definitions_common import format_definitions
from ..tables_common import common_table_main, get_header, get_footer


def format_table(results, skip_missing_names=True):
    header = get_header(
        [
            "Ensemble",
            r"$N_t \times N_s^3$",
            r"$\beta$",
            "$am_0$",
            r"$N_{\mathrm{cfg}}$",
            r"$\langle \mathcal{P} \rangle$",
            "Comment",
        ],
        hlines=1,
    )
    footer = get_footer()
    content = []
    previous_beta = None
    heavy_ps_limit = 0.45

    for result in results.sort_values(
        by=["beta", "mF"], ascending=[True, False]
    ).to_dict(orient="records"):
        if skip_missing_names and result["ensemble_name"] is None:
            continue
        if result["beta"] != previous_beta:
            previous_beta = result["beta"]
            content.append(r"\hline")
        #comment = "heavy" if result["ps_mass"].nominal_value > heavy_ps_limit else ""
        content.append(
            " & ".join(
                [
                    result["ensemble_name"],
                    f"${result['Nt']} \\times {result['Ns']}^3$",
                    f"{result['beta']}",
                    f"{result['mF']}",
                    f"{result['Ncfg']}",
                    f"{result['avg_plaquette']:.02uSL}",
                    "comment",
                ]
            )
            + r" \\"
        )
    definitions = format_definitions({"HeavyPSMassLimit": heavy_ps_limit})
    return header + "\n".join(content) + footer, definitions


if __name__ == "__main__":
    common_table_main(format_table, definitions=True)
