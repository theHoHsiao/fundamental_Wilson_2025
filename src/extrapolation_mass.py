#!/usr/bin/env python3

from functools import partial

from .fitting import global_meson_fit, mass_square_fit_form, split_means_samples
from .extrapolation_common import get_args, get_data, dump_fit_result


def main():
    args = get_args(channels=["ps", "v", "t", "av", "at", "s", "rhoE1"])
    channel_obs_key = f"smear_{args.channel}_mass"
    data = get_data(args.data_filenames, ["w0", "smear_ps_mass", channel_obs_key])

    lat_a_means, _ = split_means_samples(data["lat_a"])
    fit_result = global_meson_fit(
        partial(mass_square_fit_form, lat_a=lat_a_means),
        data["smear_ps_mass_hat_squared"],
        data[f"{channel_obs_key}_hat_squared"],
    )

    dump_fit_result(args, "continuum_mass", fit_result, ["M", "L", "W"])


if __name__ == "__main__":
    main()
