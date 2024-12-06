#!/usr/bin/env python3

from argparse import ArgumentParser, FileType
import itertools

from flow_analysis.stats.autocorrelation import exp_autocorrelation_fit

import h5py
import numpy as np
from scipy.optimize import curve_fit
import uncertainties

from .bootstrap import basic_bootstrap, get_rng
from .dump import dump_dict
from .utils import get_index_separation


def get_skip(plaquettes):
    initial_plaquette = plaquettes[0][1]
    for index, (traj_index, plaquette) in enumerate(plaquettes):
        if plaquette != initial_plaquette:
            first_movement_index = index
            break
    else:
        return None
    tau_exp = max(
        exp_autocorrelation_fit(
            [plaquette for index, plaquette in plaquettes[first_movement_index:]]
        ).nominal_value,
        1,
    )
    num_possible = (len(plaquettes) - first_movement_index - 3 * int(tau_exp)) // int(
        tau_exp
    )

    for division_exponent in range(1, int(np.log2(num_possible))):
        data_len = num_possible // 2**division_exponent * int(tau_exp)
        left_plaquettes = [
            plaquette for _, plaquette in plaquettes[-2 * data_len : -data_len]
        ]
        right_plaquettes = [plaquette for _, plaquette in plaquettes[-data_len:]]
        difference = basic_bootstrap(right_plaquettes) - basic_bootstrap(
            left_plaquettes
        )
        if abs(difference.nominal_value) < data_len**0.5 * difference.std_dev:
            return len(plaquettes) - 2 * data_len
    else:
        return None


def double_gaussian(x, A1, mu1, sigma1, A2, mu2, sigma2):
    return A1 * np.exp(-((x - mu1) ** 2) / (2 * sigma1**2)) + A2 * np.exp(
        -((x - mu2) ** 2) / (2 * sigma2**2)
    )


def fit_histogram(plaquettes):
    """
    Fit the histogram of a series with a pair of gaussians.
    Returns [A1, mu1, sigma1], [A2, mu2, sigma2]
    """
    num_bins = 40
    counts, bin_edges = np.histogram(plaquettes, num_bins)
    bin_centres = (bin_edges[1:] + bin_edges[:-1]) / 2

    plaquette_mean = plaquettes.mean()
    plaquette_std = plaquettes.std()
    initial_guess = [len(plaquettes) / num_bins, plaquette_mean, plaquette_std] * 2
    initial_guess[1] += plaquette_std
    initial_guess[4] -= plaquette_std

    fit_result, fit_covariance = curve_fit(
        double_gaussian,
        bin_centres,
        counts,
        p0=initial_guess,
        sigma=(counts + 1) ** 0.5,
        absolute_sigma=True,
    )
    return fit_result[:3], fit_result[3:]


def mean_dwell_time(plaquettes, mu1, mu2):
    """
    Return the maximum number of consecutive elements in plaquettes
    for which the choice of which of mu1 or mu2 is closer remains constant.
    """
    offset = plaquettes - (mu1 + mu2) / 2
    mode = np.astype(offset / np.abs(offset), int)
    dwell_times = np.asarray(
        [sum(1 for i in group) for _, group in itertools.groupby(mode)]
    )
    return uncertainties.ufloat(dwell_times.max(), dwell_times.std())


def bimodal_autocorrelation_time(plaquettes_list):
    """
    Estimate the autocorrelation time of a series that may have
    a slow, sharp oscillation between two modes.
    """
    plaquettes = np.asarray([plaquette for index, plaquette in plaquettes_list])
    (_, mu1, sigma1), (_, mu2, sigma2) = fit_histogram(plaquettes)
    if abs(mu2 - mu1) < max(sigma1, sigma2):
        return exp_autocorrelation_fit(plaquettes)
    else:
        return mean_dwell_time(plaquettes, mu1, mu2)


def get_args():
    parser = ArgumentParser(
        description="Compute mean plaquette and autocorrelation time from HMC history"
    )

    parser.add_argument("h5file", help="The file to read")
    parser.add_argument(
        "--ensemble_name",
        default=None,
        help="Name of the ensemble to analyse. Only used for tagging output.",
    )
    parser.add_argument(
        "--beta",
        type=float,
        default=None,
        help="The beta value of the ensemble to analyse",
    )
    parser.add_argument(
        "--mAS",
        type=float,
        default=None,
        help="The antisymmetric fermion mass of the ensemble to analyse",
    )
    parser.add_argument(
        "--Nt",
        type=int,
        default=None,
        help="The temporal extent of the ensemble to analyse",
    )
    parser.add_argument(
        "--Ns",
        type=int,
        default=None,
        help="The spatial extent of the ensemble to analyse",
    )
    parser.add_argument(
        "--start",
        choices=["unit", "random"],
        default=None,
        help="The starting state of the ensemble to analyse",
    )
    parser.add_argument(
        "--output_file",
        type=FileType("w"),
        default="-",
        help=(
            "Where to output the mean and uncertainty of the average plaquette."
            "(Defaults to stdout.)"
        ),
    )
    return parser.parse_args()


def get_plaquette(filename, Nt, Ns, beta, mass, start):
    data = h5py.File(filename, "r")
    ensemble_name = f"hmc_{Nt}x{Ns}x{Ns}x{Ns}b{beta}m{mass}_{start}"
    ens = data[ensemble_name]
    plaquettes = [
        (trajectory, plaquette)
        for trajectory, plaquette in zip(ens["trajectory"], ens["plaquette"])
        if not np.isnan(plaquette)
    ]
    result = {
        "Nt": Nt,
        "Ns": Ns,
        "beta": beta,
        "mAS": mass,
        "start": start,
    }

    separation = get_index_separation([index for index, plaquette in plaquettes])
    skip = get_skip(plaquettes)
    if skip is None:
        result["tau_exp_plaq"] = np.nan
        result["plaq_therm"] = np.nan
        result["avg_plaquette"] = uncertainties.ufloat(np.nan, np.nan)
        return result

    result["tau_exp_plaq"] = (
        bimodal_autocorrelation_time(plaquettes[skip:]) * separation
    )
    result["plaq_therm"] = min(
        int(10 * result["tau_exp_plaq"].nominal_value), len(plaquettes) // 2
    )
    plaquettes_subset = [
        plaquette
        for index, plaquette in plaquettes
        if index >= skip + result["plaq_therm"]
        and index % max(1, int(result["tau_exp_plaq"].nominal_value)) == 0
    ]
    result["avg_plaquette"] = basic_bootstrap(plaquettes_subset, get_rng(filename))
    return result


def main():
    args = get_args()
    data = get_plaquette(args.h5file, args.Nt, args.Ns, args.beta, args.mAS, args.start)
    data["ensemble_name"] = args.ensemble_name
    dump_dict(data, args.output_file)


if __name__ == "__main__":
    main()
