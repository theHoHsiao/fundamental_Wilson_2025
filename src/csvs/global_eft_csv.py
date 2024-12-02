#!/usr/bin/env python3

import numpy as np

from ..csv_common import standard_csv_main


def fit_column_index(column):
    preferred_ordering = [
        "group_family",
        "Nc",
        "channel",
        "fit_type",
        "continuum",
        "L",
        "W",
        "mAS_min",
        "mAS_max",
    ]
    if column in preferred_ordering:
        return (0, preferred_ordering.index(column))
    else:
        return (1, column)


def merge_channels(data):
    result = []
    for datum in data.to_dict(orient="records"):
        value = [value for value in datum.values() if not np.isnan(value.nominal_value)]
        if len(value) > 1:
            raise ValueError("Multiple observables found in one fit result!")
        elif len(value) == 0:
            result.append(np.nan)
        else:
            result.append(value.pop())
    return result


def process_df(data):
    combined_columns = ["M", "F", "R_vdfps"]
    data["continuum"] = merge_channels(data[combined_columns])
    return data.drop(columns=combined_columns).sort_index(
        axis="columns", key=lambda c: c.map(fit_column_index)
    )


if __name__ == "__main__":
    standard_csv_main(index_name=None, df_processor=process_df)
