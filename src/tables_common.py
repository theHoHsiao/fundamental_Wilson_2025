#!/usr/bin/env python3

from argparse import ArgumentParser, FileType
from .dump import read_files
from .provenance import text_metadata, get_basic_metadata


def split_ensemble_name(name):
    """
    Split an ensemble name, e.g. ASB0M1L3, into its component parts for sorting.
    """
    name_start, name_end = name.split("M")
    if "L" in name_end:
        mass_idx, Ns = name_end.split("L")
        return name_start, int(mass_idx), int(Ns)
    else:
        # Use an arbitrarily long size to sort after all labeled sizes
        return name_start, int(name_end), 1000


def by_ensemble_name(column):
    """
    Use as a sort key when sorting a DataFrame by column name
    to avoid putting e.g. ASB2M10 before ASB2M2.
    """
    return column.apply(split_ensemble_name)


def get_standard_table_args(definitions=False):
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
    if definitions:
        parser.add_argument(
            "--definitions_file",
            type=FileType("w"),
            default="-",
            help="Where to place the generated definitions",
        )
    return parser.parse_args()


def common_table_main(tabulate_function, index_name="ensemble_name", definitions=False):
    args = get_standard_table_args(definitions=definitions)
    data = read_files(args.data_filenames, index_name=index_name)
    if index_name == "ensemble_name":
        data = data.sort_values(by="ensemble_name", key=by_ensemble_name)
    print(text_metadata(get_basic_metadata(), comment_char="%"), file=args.output_file)

    result = tabulate_function(data)
    if definitions:
        result, definitions = result
        print(definitions, file=args.definitions_file)

    print(result, file=args.output_file)
