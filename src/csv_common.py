#!/usr/bin/env python3

from argparse import ArgumentParser, FileType, SUPPRESS

from .dump import read_files, split_df_ufloats


def get_args(index_name):
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
        choices=["ensemble_name", "beta", "channel"],
        default=None,
        help=("Column to index on." if index_name is None else SUPPRESS),
    )
    args = parser.parse_args()
    if args.index_name is None:
        args.index_name = index_name
    return args


def standard_csv_main(index_name=None, df_processor=lambda x: x):
    args = get_args(index_name=index_name)
    data = read_files(args.data_filenames, index_name=index_name)
    data["group_family"] = "Sp"
    data["Nc"] = 4
    split_df_ufloats(df_processor(data)).to_csv(args.output_file, index=False)
