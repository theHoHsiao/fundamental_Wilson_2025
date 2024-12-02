#!/usr/bin/env python3

from ..csv_common import standard_csv_main


def fit_column_index(column):
    preferred_ordering = [
        "group_family",
        "Nc",
        "beta",
        "fit_type",
        "A",
        "B",
        "y",
        "mAS_min",
        "mAS_max",
    ]
    if column in preferred_ordering:
        return (0, preferred_ordering.index(column))
    else:
        return (1, column)


def process_df(data):
    return data.sort_index(axis="columns", key=lambda c: c.map(fit_column_index))


if __name__ == "__main__":
    standard_csv_main(index_name="beta", df_processor=process_df)
