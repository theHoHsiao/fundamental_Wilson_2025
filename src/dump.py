#!/usr/bin/env python3

from collections import defaultdict
import json

import numpy as np
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


def dump_samples(data, fp):
    to_write = {}
    for k, v in data.items():
        if isinstance(v, np.ndarray):
            to_write[k] = list(v)
        else:
            to_write[k] = v
    return json.dump(to_write, fp)


def combine_df_ufloats(df):
    result = pd.DataFrame()
    for column_name in df.columns:
        if column_name.endswith("_uncertainty"):
            continue
        elif not column_name.endswith("_value"):
            result[column_name] = df[column_name]
        elif f"{column_name[:-5]}uncertainty" not in df.columns:
            result[column_name] = df[column_name]
        else:
            result[column_name[:-6]] = df.apply(
                lambda row: ufloat(
                    row[column_name], row[f"{column_name[:-5]}uncertainty"]
                ),
                axis=1,
            )
    return result


def read_files(filenames):
    search_keys = ["Q0_value", "w0_value", "mPCAC_value"]
    data = defaultdict(list)
    for filename in filenames:
        file_data = pd.read_csv(filename)
        for key in search_keys:
            if key in file_data.columns:
                data[key].append(file_data)
                break
        else:
            raise ValueError(f"Unrecognised data in {filename}.")

    data_frames = [pd.concat(obs_data) for obs_data in data.values()]

    result = data_frames[0]
    for df in data_frames[1:]:
        result = result.join(df, on="ensemble_name")

    return combine_df_ufloats(result)
