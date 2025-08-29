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
        if isinstance(v, BootstrapSampleSet):
            v = v.to_ufloat()

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
        elif isinstance(v, BootstrapSampleSet):
            if isinstance(v.mean, np.ndarray):
                to_write[f"{k}_value"] = list(v.mean)
            else:
                to_write[f"{k}_value"] = v.mean
            to_write[f"{k}_samples"] = v.samples.tolist()
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


def split_ufloat_column(column):
    values, uncertainties = [], []
    for datum in column:
        if isinstance(datum, UFloat):
            values.append(datum.nominal_value)
            uncertainties.append(datum.std_dev)
        else:
            values.append(datum)
            uncertainties.append(np.nan)
    return pd.Series(values, name=f"{column.name}_value"), pd.Series(
        uncertainties, name=f"{column.name}_uncertainty"
    )


def split_df_ufloats(df):
    result = [pd.DataFrame()]
    for column_name in df.columns:
        if (
            df[column_name].dtype == "O"
            and df[column_name].map(lambda datum: isinstance(datum, UFloat)).any()
        ):
            if (
                f"{column_name}_value" in df.columns
                or f"{column_name}_uncertainty" in df.columns
            ):
                raise ValueError(f"Conflicting columns: {column_name}")
            result.extend(split_ufloat_column(df[column_name]))
        else:
            result.append(df[column_name])
    return pd.concat(result, axis=1)


class NotYetFound:
    """None-like class that will not compare equal with None."""

    pass


def make_consistent_column(columns):
    result_values = []
    column_name = columns.columns[0]
    columns.columns = range(len(columns.columns))
    not_yet_found = NotYetFound()

    for row in columns.to_dict(orient="records"):
        result_value = not_yet_found
        for value in row.values():
            if np.isnan(value):
                continue
            if result_value is not_yet_found:
                result_value = value
                continue
            if value != result_value:
                raise ValueError(f"Inconsistent data in column {column_name}")
            if isinstance(value, int) or isinstance(value, np.integer):
                result_value = value
        if result_value is not_yet_found:
            result_value = np.nan
        result_values.append(result_value)

    return pd.Series(result_values, name=column_name)


def drop_duplicate_columns(df):
    result_columns = []
    # Verify that duplicated columns are consistent
    for column in set(df.columns):
        column_subset = df[column]
        if hasattr(column_subset, "columns"):
            # Column is name duplicated, as we get more than one column when using it
            first_column = column_subset.iloc[:, 0]
            for column_idx in range(1, len(column_subset.columns)):
                if not (first_column == column_subset.iloc[:, column_idx]).all():
                    result_columns.append(make_consistent_column(column_subset))
                    break
            else:
                result_columns.append(first_column)
        else:
            result_columns.append(column_subset)

    return pd.concat(result_columns, axis=1)


# A key on which to search; only one is needed per file type.
# (More will create duplicates.)
key_observables = {
    "ensemble_name": (
        ["Q0", "w0", "mPCAC", "avg_plaquette", "tau_exp_ps_correlator", "tau_exp_ps_eff_w0"]
        + [f"{state}_mass" for state in ["ps", "v", "t", "av", "at", "s"]]
        + [f"f_{state}_decay_constant" for state in ["ps", "v", "av"]]
        + [f"f_{state}_matrix_element" for state in ["ps", "v", "av"]]
        + [
            f"gevp_f_{state}_E0_mass"
            for state in ["ps", "v", "t", "av", "at", "s", "rhoE1"]
        ]
        + [
            f"smear_{state}_Rfps"
            for state in ["ps", "v", "t", "av", "at", "s", "rhoE1"]
        ]
        + [f"{state}_Rfps" for state in ["ps", "v", "t", "av", "at", "s"]]
        + ["smear_rhoE1_Rmv"]
    ),
    "beta": ["A"],
    "channel": ["M", "F"],
}


def read_files(filenames, index_name="ensemble_name"):
    data = defaultdict(list)
    for filename in filenames:
        file_data = pd.read_csv(filename)
        if index_name is None:
            data[None].append(file_data)
        else:
            indexed_data = file_data.set_index(index_name)
            for observable in key_observables[index_name]:
                key = f"{observable}_value"
                if key in indexed_data.columns:
                    data[key].append(indexed_data)
                    break
            else:
                raise ValueError(f"Unrecognised data in {filename}.")

    data_frames = [pd.concat(obs_data) for obs_data in data.values()]

    result = drop_duplicate_columns(
        pd.concat(data_frames, axis=1).reset_index(drop=(index_name is None))
    )
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


def read_sample_files(filenames, group_key="ensemble_name"):
    results = {}
    for filename in filenames:
        file_data = read_sample_file(filename)
        if file_data.get(group_key) not in results:
            results[file_data.get(group_key)] = file_data
        else:
            target = results[file_data.get(group_key)]
            for k, v in file_data.items():
                if "samples" not in k and k in target:
                    if target[k] != v:
                        raise ValueError(f"Inconsistent metadata in {filename}")
                elif "samples" in k and not v:
                    continue
                else:
                    target[k] = v

    return list(results.values())
