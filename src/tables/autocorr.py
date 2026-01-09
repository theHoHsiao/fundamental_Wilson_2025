import numpy as np

from ..tables_common import common_table_main, get_header, get_footer


def format_table(df):
    header = get_header(
        [
            r"Ensemble",
            r"$N_{\rm cfg}$", #
            r"$N_{\rm cfg}^{\rm spec.}$",
            r"$\delta^{\rm meas}_{\mathrm{traj}}$",
            r"$\delta^{\rm spec}_{\mathrm{traj}}$",
            #r"$\tau_{\mathrm{exp}}^{\mathcal{P}}$",
            r"$\tau_{\mathrm{int}}^{\mathrm{ps-pt}}$ ",
            r"$\tau_{\mathrm{int}}^{\mathrm{ps-sm}}$ ",
            r"$\tau_{\mathrm{int}}^{\mathcal{P}}$",
            r"$\tau_{\mathrm{int}}^{w_0}$",
            r"$\tau_{\mathrm{int}}^{\mathcal{Q}}$",
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
            formatted_tau_exp_w0 = "{:.02uSL}".format(row.tau_exp_w0)
            formatted_delta_traj_w0 = "{}".format(row.delta_traj_w0)
        content.append(
            "{} & {} & {} & {} & {} & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ & {} & ${:.02uSL}$ \\\\\n".format(
                row.ensemble_name,
                row.Ncfg, 
                row.Ncfg_spec, 
                row.delta_traj_auto, #row.tau_exp_plaq, {:.02uSL} &
                row.delta_traj_spec,
                row.tau_exp_ps_correlator_point,
                row.tau_exp_ps_correlator_smear,
                row.tau_exp_plaq,
                formatted_tau_exp_w0,
                row.tau_exp_Q,
            )
        )
    return header + "".join(content) + footer


if __name__ == "__main__":
    common_table_main(format_table)
