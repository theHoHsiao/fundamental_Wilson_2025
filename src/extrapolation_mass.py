#!/usr/bin/env python3


from .fitting import meson_M2
from .extrapolation_common import get_args, get_data, dump_fit_result


def main():
    args = get_args(channels=["ps", "v", "t", "av", "at", "s", "rhoE1"])
    channel_obs_key = f"smear_{args.channel}_mass"
    data = get_data(args.data_filenames, ["w0", "ps_mass", channel_obs_key])

    fit_result = meson_M2(
        data["smear_ps_mass_hat_squared"],
        data["lat_a"],
        data[f"{channel_obs_key}_hat_squared"],
    )

    dump_fit_result(args, fit_result, ["M", "L", "W"])


if __name__ == "__main__":
    main()
