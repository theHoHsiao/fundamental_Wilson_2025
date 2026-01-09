#!/usr/bin/env python3

from functools import partial

from .fitting import global_meson_fit, mass_square_fit_form, split_means_samples
from .extrapolation_common import get_args, get_data, dump_fit_result


def main():
    args = get_args(channels=["ps", "v", "t", "av", "at", "s", "rhoE1"])
    channel_obs_key = f"gevp_f_{args.channel}_E0_mass"
    data = get_data(args.data_filenames, ["w0", "gevp_f_ps_E0_mass", channel_obs_key])
    #print(data.keys())

    lat_a_means, _ = split_means_samples(data["lat_a"])
    fit_result = global_meson_fit(
        partial(mass_square_fit_form, lat_a=lat_a_means),
        data["gevp_f_ps_E0_mass_hat_squared"],
        data[f"{channel_obs_key}_hat_squared"],
    )

    dump_fit_result(args, "continuum_Wilson_ChiPT", fit_result, ["M", "Lm", "Wm"])


if __name__ == "__main__":
    main()
