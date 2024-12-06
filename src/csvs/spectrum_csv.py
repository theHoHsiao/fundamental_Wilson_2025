#!/usr/bin/env python3

from ..csv_common import standard_csv_main


def common_column_index(column_name):
    standard_columns = [
        "ensemble_name",
        "group_family",
        "Nc",
        "Nt",
        "Ns",
        "beta",
        "nF",
        "nAS",
        "nS",
        "nADJ",
        "mF",
        "mAS",
        "mS",
        "mADJ",
        "start",
        "init_conf",
        "final_conf",
        "delta_conf_spectrum",
        "delta_conf_w0",
    ]
    perturbation = 0
    if column_name.endswith("_value"):
        column_name = column_name[:-6]
    if column_name.endswith("_uncertainty"):
        column_name = column_name[:-12]
        perturbation = 1

    if column_name in standard_columns:
        return (0, standard_columns.index(column_name), perturbation)
    if column_name.startswith("N"):
        return (1, column_name, perturbation)
    if "tau" in column_name:
        return (2, column_name, perturbation)
    if "Q" in column_name:
        return (3, column_name, perturbation)
    if "plaquette" in column_name:
        return (4, column_name, perturbation)
    if column_name == "mPCAC":
        return (5, column_name, perturbation)
    if any(
        [
            column_name.endswith(key)
            for key in ["mass", "decay_constant", "matrix_element", "Rfps", "chisquare"]
        ]
    ) and not column_name.startswith("smear_"):
        return (6, column_name, perturbation)
    if column_name.startswith("smear_"):
        return (7, column_name, perturbation)
    return (8, column_name, perturbation)


def process_df(data):
    # Currently this workflow is only used for nAS=3,
    # so hardcode this here.
    # This should be generalised elsewhere in the workflow
    # if this is applied to other theories.
    if "nAS" not in data.columns:
        data["nAS"] = 3
    return data.sort_index(axis="columns", key=lambda c: c.map(common_column_index))


if __name__ == "__main__":
    standard_csv_main(df_processor=process_df, index_name="ensemble_name")
