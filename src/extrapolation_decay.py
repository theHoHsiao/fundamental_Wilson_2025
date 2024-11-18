#!/usr/bin/env python3


from .fitting import meson_M2
from .extrapolation_common import get_args, get_data, dump_fit_result


def main():
    args = get_args(channels=["ps", "v", "av"])
    channel_obs_key = f"{args.channel}_decay_constant"
    data = get_data(args.data_filenames, ["w0", "ps_mass", channel_obs_key])

    fit_result = meson_M2(
        data["ps_mass_hat_squared"], data["lat_a"], data[channel_obs_key]
    )

    dump_fit_result(
        args,
        fit_result,
        [f"{var}_{args.channel}" for var in ["F", "L", "W"]],
    )


if __name__ == "__main__":
    main()
