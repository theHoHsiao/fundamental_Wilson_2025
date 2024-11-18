#!/usr/bin/env python3


from .fitting import meson_M2
from .extrapolation_common import get_args, get_data, dump_fit_result


def main():
    args = get_args(channels=["ps", "v", "av"])
    channel_obs_key = f"{args.channel}_mass"
    data = get_data(
        args.data_filenames, ["w0", "ps_mass", "ps_decay_constant", channel_obs_key]
    )

    fit_result = meson_M2(
        data["ps_mass_hat_squared"],
        data["lat_a"],
        data[f"{channel_obs_key}_over_ps_decay_const"],
    )

    dump_fit_result(
        args,
        fit_result,
        [f"{var}_{args.channel}dfps" for var in ["R", "L", "W"]],
    )


if __name__ == "__main__":
    main()
