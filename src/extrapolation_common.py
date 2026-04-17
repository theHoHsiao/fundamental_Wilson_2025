#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

import numpy as np
import inspect
from functools import partial

from .dump import dump_dict, dump_samples, read_sample_files
from .fitting import mass_square_fit_form, mass_square_a2_fit_form, mass_square_m4_fit_form, mass_square_a2_m4_fit_form, mass_square_a2_am2_fit_form, mass_square_full_fit_form
from .fitting import global_meson_fit


def get_args(channels=None, beta=False, ansatz=None):
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
    if ansatz is not None:
        parser.add_argument(
            "--ansatz",
            choices=["a", "a2", "m4", "a2_m4", "a2_am2", "full"],
            default=None,
            help="Ansatz to use for extrapolation",
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
        target[f"log_{observable}_squared_over_mPCAC"] = np.log(
            obs_value**2 / source["mPCAC_samples"]
        )
        created_keys.append(f"log_{observable}_squared_over_mPCAC")

    return created_keys


def get_data(filenames, observables, beta=None, cutoff=10, cutoff_observable=None):
    data = read_sample_files(filenames)
    results = []
    extra_observables = set()
    for datum in data:
        datum_result = {}
        if beta is not None and (datum.get("beta") != beta):
            continue
        

        if cutoff_observable is not None and datum[f"{cutoff_observable}_samples"].mean > cutoff:
            #print(datum[f"{cutoff_observable}_samples"].mean)
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
            f"Wrong number of names ({len(names)}) "
            f"for fit result ({len(values)} values)"
        )
        raise ValueError(message)


def dump_fit_result(args, fit_type, fit_result, names, **extra_columns):
    fit_values, chisquare = fit_result
    check_name_value_lengths(names, fit_values)

    keys = {
        name: getattr(args, name) for name in ["beta", "channel"] if hasattr(args, name)
    }

    results = {**keys, **{name: value for name, value in zip(names, fit_values)}}
    results["fit_type"] = fit_type

    dump_dict(
        {**results, "chisquare": chisquare, **extra_columns},
        args.output_file_mean,
    )
    if args.output_file_samples:
        dump_samples(results, args.output_file_samples)


def dump_fit_result_chisquare(args, fit_type, fit_result, names, **extra_columns):
    fit_values, chisquare = fit_result
    check_name_value_lengths(names, fit_values)

    keys = {
        name: getattr(args, name) for name in ["beta", "channel"] if hasattr(args, name)
    }

    results = {**keys, **{name: value for name, value in zip(names, fit_values)}}
    results["fit_type"] = fit_type

    dump_dict(
        {**results, "chisquare": chisquare, **extra_columns},
        args.output_file_mean,
    )
    if args.output_file_samples:
        dump_samples({**results, f"chisquare_{args.ansatz}": chisquare}, args.output_file_samples)



def initialize_fit_parameters(fit_form, guess):
    num_fit = get_num_params(fit_form)
    x0 = np.zeros(num_fit)

    x0[:3] = guess

    return x0


def initialize_fit_parameters_from_a(fit_form, data, lat_a_means, channel_obs_key, post_guess=None):

    init_param_a = global_meson_fit(
        partial(mass_square_fit_form, lat_a=lat_a_means),
        data["gevp_f_ps_E0_mass_hat_squared"],
        data[f"{channel_obs_key}_hat_squared"],
        initial_guess=post_guess,
    )
    popt_simple = [init_param_a[0][0].mean, init_param_a[0][1].mean, init_param_a[0][2].mean,]

    num_fit = get_num_params(fit_form)

    x0 = np.zeros(num_fit)

    x0[:3] = popt_simple

    #print(f"auto init params: {x0} for {channel_obs_key}")

    return x0


def get_num_params(func):
    # Get the list of parameters for the function
    sig = inspect.signature(func)
    # Subtract 1 because the first argument is the data (mass/lat_a)
    return len(sig.parameters) - 2


def ansatz_form(fit_prefix):
    return {
        "a": mass_square_fit_form,
        "a2": mass_square_a2_fit_form,
        "m4": mass_square_m4_fit_form,
        "a2_m4": mass_square_a2_m4_fit_form,
        "a2_am2": mass_square_a2_am2_fit_form,
        "full": mass_square_full_fit_form,
    }[fit_prefix]