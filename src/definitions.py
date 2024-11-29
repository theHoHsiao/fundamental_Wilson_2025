#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

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


def get_args():
    parser = ArgumentParser(description="Output a single LaTeX definition")
    parser.add_argument("name", help="Name of the command to define")
    parser.add_argument("value", help="What the command should generate")
    parser.add_argument(
        "--definitions_file",
        type=FileType("w"),
        default="-",
        help="Where to place the generated definitions",
    )
    return parser.parse_args()


def main():
    args = get_args()
    print(format_definitions({args.name: args.value}), file=args.definitions_file)


if __name__ == "__main__":
    main()
