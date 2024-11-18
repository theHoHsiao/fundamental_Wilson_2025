#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

import numpy as np

from .dump import dump_dict, dump_samples, read_sample_files


def get_args(channels=None, beta=False):
    parser = ArgumentParser()

    parser.add_argument(
        "data_filenames",
        nargs="+",
        metavar="sample_filename",
        help="Filenames of sample files used for extrapolation",
    )
    if beta:
        parser.add_argument(
            "--beta",
            type=float,
            default=None,
            help="Beta value to which to filter the data",
        )
    if channels is not None:
        parser.add_argument(
            "--channel",
            choices=channels,
            default=None,
            help="Measuring channel",
        )

    parser.add_argument("--output_file_mean", type=FileType("w"), default="-")

    parser.add_argument(
        "--output_file_samples",
        type=FileType("w"),
        default=None,
        help="Where to output the bootstrap samples for fitting results",
    )

    return parser.parse_args()


def compute_derived_spectrum(source, target, observable):
    obs_value = target[observable]

    target[f"{observable}_squared"] = obs_value**2
    target[f"log_{observable}_squared"] = np.log(obs_value**2)
    created_keys = [f"{observable}_squared", f"log_{observable}_squared"]

    if "w0_samples" in source:
        target[f"{observable}_hat_squared"] = (source["w0_samples"] * obs_value) ** 2
        created_keys.append(f"{observable}_hat_squared")
    if "ps_decay_constant_samples" in source:
        target[f"{observable}_over_ps_decay_constant"] = (
            obs_value / source["ps_decay_constant_samples"]
        )
        created_keys.append(f"{observable}_over_ps_decay_constant")
    if "mPCAC_samples" in source:
        target[f"{observable}_squared_over_mPCAC"] = np.log(
            obs_value**2 / source["mPCAC_samples"]
        )
        created_keys.append(f"log_{observable}_squared_over_mpcac")

    return created_keys


def prepare_data(filenames, observables, beta=None):
    data = read_sample_files(filenames)
    results = []
    extra_observables = set()
    for datum in data:
        datum_result = {}
        if beta is not None and (datum.get("beta") != beta):
            continue

        for observable in observables:
            if observable in datum:
                # This is not a sample set, so just copy across as-is
                datum_result[observable] = datum[observable]
                continue

            key = f"{observable}_samples"
            if key not in datum:
                break

            obs_value = datum[key]
            if np.isnan(obs_value.samples.mean()):
                break
            datum_result[observable] = obs_value

            if observable == "w0":
                datum_result["lat_a"] = 1 / obs_value
                extra_observables.add("lat_a")

            if observable.endswith("_mass") or observable.endswith("_decay_constant"):
                extra_observables.update(
                    compute_derived_spectrum(datum, datum_result, observable)
                )
        else:
            # Only append if all required observables were found
            results.append(datum_result)

    return {
        observable: np.asarray([datum[observable] for datum in results])
        for observable in (observables + list(extra_observables))
    }


def check_name_value_lengths(names, values):
    if len(values) != len(names):
        message = (
            "Wrong number of names ({len(names)}) "
            "for fit result ({len(fit_values)} values)"
        )
        raise ValueError(message)


def dump_fit_result(args, fit_result, names, **extra_columns):
    fit_values, chisquare = fit_result
    check_name_value_lengths(names, fit_values)

    keys = {
        name: getattr(args, name) for name in ["beta", "channel"] if hasattr(args, name)
    }

    results = {**keys, **{name: value for name, value in zip(names, fit_values)}}

    dump_dict(
        {"chisquare": chisquare, **results, **extra_columns},
        args.output_file_mean,
    )
    if args.output_file_samples:
        dump_samples(results, args.output_file_samples)
