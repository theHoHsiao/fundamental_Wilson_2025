#!/usr/bin/env python3

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
