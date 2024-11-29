#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

from .dump import read_sample_files
from .provenance import get_basic_metadata, text_metadata


def format_definitions(definitions):
    formatted_definitions = [text_metadata(get_basic_metadata(), comment_char="%")]
    formatted_definitions.extend(
        [
            r"\newcommand{{\{}}}{{{}}}".format(name, value)
            for name, value in definitions.items()
        ]
    )
    return "\n".join(formatted_definitions)


def get_standard_definition_args():
    parser = ArgumentParser()

    parser.add_argument(
        "data_filenames", nargs="+", help="Filenames of result files to tabulate"
    )
    parser.add_argument(
        "--output_file",
        type=FileType("w"),
        default="-",
        help="Where to place the definitions",
    )
    return parser.parse_args()


def common_definitions_main(definition_function, group_key="channel"):
    args = get_standard_definition_args()
    data = read_sample_files(args.data_filenames, group_key=group_key)
    print(format_definitions(definition_function(data)), file=args.output_file)
