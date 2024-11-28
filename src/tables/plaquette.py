#!/usr/bin/env python3


from ..tables_common import common_table_main


def format_table(results, skip_missing_names=True):
    header = (
        "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|}\n\\hline\\hline\n"
        + " & ".join(
            [
                "Ensemble",
                "$N_t \\times N_s^3$",
                r"$\beta$",
                "$am_0$",
                r"$N_{\mathrm{cfg}}$",
                r"$\langle \mathcal{P} \rangle$",
                "Comment",
            ]
        )
        + " \\\\\n"
    )
    footer = "\n\\hline\\hline\n\\end{tabular}"
    content = []
    previous_beta = None

    for result in results.sort_values(
        by=["beta", "mAS"], ascending=[True, False]
    ).to_dict(orient="records"):
        if skip_missing_names and result["ensemble_name"] is None:
            continue
        if result["beta"] != previous_beta:
            previous_beta = result["beta"]
            content.append(r"\hline")
        content.append(
            " & ".join(
                [
                    result["ensemble_name"],
                    f"${result['NT']} \\times {result['NS']}^3$",
                    f"{result['beta']}",
                    f"{result['mAS']}",
                    f"{result['Ncfg_spectrum']}",
                    f"{result['avg_plaquette']:.02uSL}",
                    "tbc",
                ]
            )
            + r" \\"
        )
    return header + "\n".join(content) + footer


if __name__ == "__main__":
    common_table_main(format_table)
