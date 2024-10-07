#!/usr/bin/env python3

from collections import defaultdict
import json

import numpy as np
import pandas as pd
from uncertainties import ufloat, UFloat

from .bootstrap import BootstrapSampleSet


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
        file_data = pd.read_csv(filename).set_index("ensemble_name")
        for key in search_keys:
            if key in file_data.columns:
                data[key].append(file_data)
                break
        else:
            raise ValueError(f"Unrecognised data in {filename}.")

    data_frames = [pd.concat(obs_data) for obs_data in data.values()]

    result = pd.concat(data_frames, axis=1).reset_index()
    return combine_df_ufloats(result)


def read_sample_file(filename):
    with open(filename, "r") as f:
        raw_data = json.load(f)

    data = {}
    for k, v in raw_data.items():
        if isinstance(v, list):
            if len(v) == 0:
                continue

            data[k] = BootstrapSampleSet(v)
        else:
            data[k] = v

    return data


def read_sample_files(filenames):
    results = {}
    for filename in filenames:
        file_data = read_sample_file(filename)
        if file_data["ensemble_name"] not in results:
            results[file_data["ensemble_name"]] = file_data
        else:
            target = results[file_data["ensemble_name"]]
            for k, v in file_data.items():
                if "samples" not in k and k in target:
                    if target[k] != v:
                        raise ValueError(f"Inconsistent metadata in {filename}")
                elif "samples" in k and not v:
                    continue
                else:
                    target[k] = v

    return list(results.values())
