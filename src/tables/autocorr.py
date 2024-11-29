import numpy as np

from ..tables_common import common_table_main, get_header, get_footer


def format_table(df):
    header = get_header(
        [
            r"Ensemble",
            r"$\delta_{\mathrm{traj}}$",
            r"$\tau_{\mathrm{exp}}^{\mathcal{P}}$",
            r"$\tau_{\mathrm{exp}}^{\mathrm{ps}}$ ",
            r"$\delta_{\mathrm{traj}}^{\mathrm{GF}}$ ",
            r"$\tau_{\mathrm{exp}}^{w_0}$",
            r"$\tau_{\mathrm{exp}}^{\mathcal{Q}}$",
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
