#!/usr/bin/env python3

# from uncertainties import ufloat, UFloat

# from .bootstrap import BootstrapSampleSet

from ..tables_common import ensemble_table_main
#import numpy as np


def format_table(df):
    header = (
        "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|}\n"
        "\\hline\\hline\n"
        r"Ensemble & $\beta$ & $m_0$ & $N_t$ & $N_s$ & $N_{\rm config.}$ & $\delta_{\rm traj.}$ & $N_{\rm bin}$"
        r"& $\tau_{\mathrm{int}}^{\langle P \rangle}$ & $\tau_{\mathrm{int}}^{w_0/a}$ & $ \langle P \rangle$ & $w_0/a$ & "
        r"$Q_0$ & $\sigma_Q$  \\"
        "\n\\hline"
    )
    footer = "\\hline\\hline\n\\end{tabular}"
    content = []
    previous_prefix = None
    for row in df.itertuples():
        if (next_prefix := row.ensemble_name[:4]) != previous_prefix:
            previous_prefix = next_prefix
            content.append("\\hline\n")

        if row.Hasenbush:
            ast = "$^\\ast$"
        else:
            ast = ""

        content.append(
            (
                "{}"+ast+" & {} & {} & {} & {} & {} & {} & {} &"
                " ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$"
                "\\\\\n"
            ).format(
                row.ensemble_name,
                row.beta,
                row.mF,
                int(row.Nt),
                int(row.Ns),
                row.Ncfg,
                row.delta_traj_spec,
                row.bin_size,
                row.tau_exp_plaq,
                row.tau_exp_w0,
                row.avg_plaquette,
                row.w0,
                row.Q0,
                row.sigma_Q,
            )
        )

    return header + "".join(content) + footer


if __name__ == "__main__":
    ensemble_table_main(format_table)
