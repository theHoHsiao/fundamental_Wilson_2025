import numpy as np

from ..tables_common import common_table_main


def format_table(df):
    header = (
        "\\begin{tabular}{|c|c|c|c|c|c|c|} \n \\hline \\hline \n"
        r"Ensemble & $\delta_{\mathrm{traj}}$ & "
        r"$\tau_{\mathrm{exp}}^{\mathcal{P}}$ & $\tau_{\mathrm{exp}}^{\mathrm{ps}}$ "
        r"& $\delta_{\mathrm{traj}}^{\mathrm{GF}}$ "
        r"& $\tau_{\mathrm{exp}}^{w_0}$ & $\tau_{\mathrm{exp}}^{\mathcal{Q}}$"
        "\\\\ \\hline \n"
    )
    footer = "\\hline \n \\end{tabular}"
    content = []
    for row in df.itertuples():
        if np.isnan(row.tau_exp_w0.nominal_value):
            formatted_tau_exp_w0 = r"$\cdots$"
            formatted_delta_traj_w0 = r"$\cdots$"
        else:
            formatted_tau_exp_w0 = "{:.0f}".format(row.tau_exp_w0.nominal_value)
            formatted_delta_traj_w0 = "{}".format(row.delta_traj_w0)
        content.append(
            "{} & {} & {:.1f} & {:.1f} & {} & {} & {:.0f} \\\\\n".format(
                row.ensemble_name,
                row.delta_traj_spectrum,
                row.tau_exp_plaq.nominal_value,
                row.tau_exp_ps_correlator.nominal_value,
                formatted_delta_traj_w0,
                formatted_tau_exp_w0,
                row.tau_exp_Q.nominal_value,
            )
        )
    return header + "".join(content) + footer


if __name__ == "__main__":
    common_table_main(format_table)
