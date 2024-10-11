#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

from ..dump import read_files


def get_args():
    parser = ArgumentParser(
        description="Take summary ensemble data and output a LaTeX table"
    )
    parser.add_argument(
        "data_filenames", nargs="+", help="Filenames of plaquette result files"
    )
    parser.add_argument("--output_file", type=FileType("w"), default="-")
    return parser.parse_args()


def tabulate(results, skip_missing_names=True):
    header = (
        "\\begin{tabular}{ccccccccc}\n\\hline\\hline\n"
        + " & ".join(
            [
                "Ensemble",
                "$N_t \\times N_s^3$",
                r"$\beta$",
                "$am_0$",
                r"$N_{\mathrm{cfg}}$",
                r"$\delta_{\mathrm{traj}}$",
                r"$\langle \mathcal{P} \rangle$",
                r"$\tau_{\mathrm{exp}}^{\mathcal{P}}$",
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
                    f"{result['delta_traj_spectrum']}",
                    f"{result['avg_plaquette']:.02uSL}",
                    f"{result['plaq_autocorr']:.02uSL}",
                    "tbc",
                ]
            )
            + r" \\"
        )
    return header + "\n".join(content) + footer


def main():
    args = get_args()
    data = read_files(args.data_filenames)
    print(tabulate(data), file=args.output_file)


if __name__ == "__main__":
    main()
