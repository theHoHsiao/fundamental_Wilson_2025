#!/usr/bin/env python3
from collections import defaultdict

import numpy as np
import pandas as pd
# from uncertainties import ufloat, UFloat

# from .bootstrap import BootstrapSampleSet

from argparse import ArgumentParser, FileType
from ..dump import combine_df_ufloats


def get_args():
    parser = ArgumentParser()

    parser.add_argument(
        "data_filenames", nargs="+", help="Filenames of meson extraction result files"
    )
    parser.add_argument("--output_file", type=FileType("w"), default="-")
    return parser.parse_args()


def format_table(df):
    header = (
        "\\begin{tabular}{|c|c|c|c|c|c|c|}\n"
        "\\hline\\hline\n"
        r"Ensemble & $Q_0$ & $\sigma_Q$ & $\tau_{\mathrm{exp}}$ & "
        r"$N_{\mathrm{traj}}^{\mathrm{GF}}$ & $\delta_{\mathrm{traj}}^{\mathrm{GF}}$ & "
        "$w_0 / a$ \\\\\n"
        r"\hline"
    )
    footer = "\\hline\\hline\n\\end{tabular}"
    content = []
    previous_prefix = None
    for row in df.itertuples():
        if (next_prefix := row.ensemble_name[:4]) != previous_prefix:
            previous_prefix = next_prefix
            content.append("\\hline\n")

        if np.isnan(row.w0.nominal_value) or np.isnan(row.w0.std_dev):
            w0 = r"\cdots"
            num_configs = r"$\cdots"
            trajectory_step = r"$\cdots"
        else:
            w0 = f"{row.w0:.02uSL}"
            num_configs = row.num_configs
            trajectory_step = row.trajectory_step

        content.append(
            (
                "{} & ${:.02uSL}$ & ${:.02uSL}$ & ${:.02uSL}$ & {} & {} & ${}$ \\\\\n"
            ).format(
                row.ensemble_name,
                row.Q0,
                row.sigma_Q,
                row.tau_exp_Q,
                num_configs,
                trajectory_step,
                w0,
            )
        )

    return header + "".join(content) + footer


def read_files(filenames):
    search_keys = [
        "ps_mass_value",
        "ps_matrix_element_value",
        "v_mass_value",
        "v_matrix_element_value",
        "t_mass_value",
        "t_matrix_element_value",
        "av_mass_value",
        "av_matrix_element_value",
        "at_mass_value",
        "at_matrix_element_value",
        "s_mass_value",
        "s_matrix_element_value",
    ]
    data = defaultdict(list)
    for filename in filenames:
        file_data = pd.read_csv(filename).set_index("ensemble_name")
        for key in search_keys:
            if key in file_data.columns:
                data[key].append(file_data)
                break
        else:
            raise ValueError(f"Unrecognised data in {filename}.")

    data_frames = [pd.concat(obs_data) for obs_data in data.values()]

    result = pd.concat(data_frames, axis=1).reset_index()
    return combine_df_ufloats(result)


def print_table(csv):
    select_ens = csv[csv.use_in_main_plots].ensemble_name.values

    header = "\\begin{tabular}{|" + "c|" * 7 + "} \\hline \\hline \n"

    for ens in select_ens:
        tmp_ens = csv[csv.ensemble_name == ens]
        Nc = tmp_ens.Nc.values[0]
        beta = tmp_ens.beta.values[0]
        nAS = tmp_ens.nAS.values[0]
        mAS = tmp_ens.mAS.values[0]
        Ns = tmp_ens.Ns.values[0]
        Nt = tmp_ens.Nt.values[0]
        dir_template = f"Sp{Nc}b{beta}nAS{nAS}mAS{mAS}T{Nt}L{Ns}"
        print(dir_template)

        # for ch in ["ps", "v", "t", "av", "at", "s"]:

    return header


def main():
    args = get_args()
    data = read_files(args.data_filenames)
    print(data.keys())
    # print(format_table(data.sort_values(by="ensemble_name")), file=args.output_file)


if __name__ == "__main__":
    main()
