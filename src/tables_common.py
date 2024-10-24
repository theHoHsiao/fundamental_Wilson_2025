#!/usr/bin/env python3


def by_ensemble_name(column):
    """
    Use as a sort key when sorting a DataFrame by column name
    to avoid putting e.g. ASB2M10 before ASB2M2.
    """
    return column.apply(lambda e: ((elems := e.split("M"))[0], int(elems[1])))
