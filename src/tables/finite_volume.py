#!/usr/bin/env python3


from ..definitions import format_definitions
from ..tables_common import common_table_main, get_header, get_footer


def format_table(df):
    beta = set(df.beta)
    mAS = set(df.mAS)
    delta_traj = set(df.delta_traj_spectrum)

    if len(beta) != 1 or len(mAS) != 1 or len(delta_traj) != 1:
        raise ValueError("Inconsistent parameters in finite volume study.")

    definitions = format_definitions(
        {
            "FiniteVolumeStudyBeta": beta.pop(),
            "FiniteVolumeStudyMass": mAS.pop(),
            "FiniteVolumeStudyDeltaTraj": delta_traj.pop(),
        }
    )

    m_ps_inf = df[df.Ns == max(df.Ns)].ps_mass.values[0]

    header = get_header(
        [
            r"$N_t \times N_s^3$",
            r"$N_{\mathrm{cfg}}$ ",
            r"$\langle\mathcal{P}\rangle$",
            r"$am_{\mathrm{PCAC}}$ ",
            r"$am_{\mathrm{ps}}$",
            r"$am_{\mathrm{v}}$ ",
            r"$m_{\mathrm{ps}}^{\mathrm{inf}} L$",
        ],
        hlines=1,
    )
    footer = get_footer()
    content = []
    for row in df.sort_values(by="Ns").itertuples():
        content.append(
            (
                r"${} \times {}$ & {} "
                r"& {:.02uSL} & {:.02uSL} & {:.02uSL} & {:.02uSL} & {:.02uSL} \\"
            ).format(
                row.Nt,
                row.Ns,
                row.Ncfg_spectrum,
                row.avg_plaquette,
                row.mPCAC,
                row.ps_mass,
                row.v_mass,
                m_ps_inf * row.Ns,
            )
        )
    return header + "\n".join(content) + footer, definitions


if __name__ == "__main__":
    common_table_main(format_table, definitions=True)
