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
        elif isinstance(v, np.int64):
            to_write[k] = int(v)
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


def drop_duplicate_columns(df):
    # Verify that duplicated columns are consistent
    for column in set(df.columns):
        column_subset = df[column]
        if hasattr(column_subset, "columns"):
            # Column is name duplicated, as we get more than one column when using it
            first_column = column_subset.iloc[:, 0]
            for column_idx in range(1, len(column_subset.columns)):
                if not (first_column == column_subset.iloc[:, column_idx]).all():
                    raise ValueError(f"Inconsistent data for column {column}.")

    # Drop the duplicated columns
    return df.loc[:, ~df.columns.duplicated()].copy()


def read_files(filenames):
    # A key on which to search; only one is needed per file type.
    # (More will create duplicates.)
    search_keys = [
        "Q0_value",
        "w0_value",
        "mPCAC_value",
        "avg_plaquette_value",
        "ps_mass_value",
    ]

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

    result = drop_duplicate_columns(pd.concat(data_frames, axis=1).reset_index())
    return combine_df_ufloats(result)


def read_sample_file(filename):
    with open(filename, "r") as f:
        raw_data = json.load(f)

    data = {}
    samples_fields = [key for key in raw_data if key.endswith("_samples")]
    for samples_field in samples_fields:
        value_field = samples_field.replace("_samples", "_value")
        if value_field not in raw_data:
            raise ValueError("Bootstrap samples with no central value")
        data[samples_field] = BootstrapSampleSet(
            raw_data.pop(value_field), raw_data.pop(samples_field)
        )

    return {**data, **raw_data}


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


def read_extp_sample_files(filenames):
    results = {}
    for filename in filenames:
        file_data = read_sample_file(filename)
        if file_data["channel"] not in results:
            results[file_data["channel"]] = file_data
        else:
            target = results[file_data["channel"]]
            for k, v in file_data.items():
                if "samples" not in k and k in target:
                    if target[k] != v:
                        raise ValueError(f"Inconsistent metadata in {filename}")
                elif "samples" in k and not v:
                    continue
                else:
                    target[k] = v

    return list(results.values())
