#!/usr/bin/env python3


def format_definitions(definitions):
    formatted_definitions = [
        r"\newcommand{{\{}}}{{{}}}".format(name, value)
        for name, value in definitions.items()
    ]
    return "\n".join(formatted_definitions)
