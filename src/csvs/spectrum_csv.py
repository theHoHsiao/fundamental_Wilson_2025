#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

from ..dump import split_df_ufloats, read_files


def get_args():
    parser = ArgumentParser(
        description="Collate many small CSV files into one large one for sharing."
    )
    parser.add_argument(
        "data_filenames", nargs="+", help="Filenames of result files to tabulate"
    )
    parser.add_argument(
        "--output_file",
        type=FileType("w"),
        default="-",
        help="Where to place the resulting CSV",
    )
    parser.add_argument(
        "--index_name",
        required=True,
        choices=["ensemble_name", "beta", "channel"],
        help="Column to index on.",
    )
    return parser.parse_args()


def common_column_order(columns):
    return columns.map(common_column_index)


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
    if any(
        [
            column_name.endswith(key)
            for key in ["mass", "decay_constant", "matrix_element", "chisquare"]
        ]
    ):
        return (5, column_name, perturbation)
    if column_name.startswith("smear_"):
        return (6, column_name, perturbation)
    return (7, column_name, perturbation)


def main():
    args = get_args()
    data = read_files(args.data_filenames, index_name=args.index_name)
    data["group_family"] = "Sp"
    data["Nc"] = 4
    split_df_ufloats(data.sort_index(axis="columns", key=common_column_order)).to_csv(
        args.output_file, index=False
    )


if __name__ == "__main__":
    main()
