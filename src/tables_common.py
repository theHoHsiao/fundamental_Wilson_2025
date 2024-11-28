#!/usr/bin/env python3

from argparse import ArgumentParser, FileType
from .dump import read_files


def by_ensemble_name(column):
    """
    Use as a sort key when sorting a DataFrame by column name
    to avoid putting e.g. ASB2M10 before ASB2M2.
    """
    return column.apply(lambda e: ((elems := e.split("M"))[0], int(elems[1])))


def get_standard_table_args():
    parser = ArgumentParser()

    parser.add_argument(
        "data_filenames", nargs="+", help="Filenames of result files to tabulate"
    )
    parser.add_argument(
        "--output_file",
        type=FileType("w"),
        default="-",
        help="Where to place the table",
    )
    return parser.parse_args()


def common_table_main(tabulate_function, index_name="ensemble_name"):
    args = get_standard_table_args()
    data = read_files(args.data_filenames, index_name=index_name)
    if index_name == "ensemble_name":
        data = data.sort_values(by="ensemble_name", key=by_ensemble_name)
    print(tabulate_function(data), file=args.output_file)
