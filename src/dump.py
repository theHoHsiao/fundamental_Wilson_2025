#!/usr/bin/env python3

import pandas as pd
from uncertainties import ufloat, UFloat


def dump_dict(data, filename):
    to_write = {}
    for k, v in data.items():
        if isinstance(v, UFloat):
            if f"{k}_value" in data or "{k}_uncertainty" in data:
                raise ValueError("Clashing keys detected.")

            to_write[f"{k}_value"] = v.nominal_value
            to_write[f"{k}_uncertainty"] = v.std_dev
        else:
            to_write[k] = v

    pd.DataFrame([to_write]).to_csv(filename, index=False)


def combine_df_ufloats(df):
    result = pd.DataFrame()
    for column_name in df.columns:
        if column_name.endswith("_uncertainty"):
            continue
        elif not column_name.endswith("_value"):
            result[column_name] = df[column_name]
        elif not f"{column_name[:-5]}uncertainty" in df.columns:
            result[column_name] = df[column_name]
        else:
            result[column_name[:-6]] = df.apply(
                lambda row: ufloat(row[column_name], row[f"{column_name[:-5]}uncertainty"]),
                axis=1,
            )
    return result
