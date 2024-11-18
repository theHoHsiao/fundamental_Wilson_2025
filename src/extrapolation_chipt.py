#!/usr/bin/env python3

from .fitting import meson_beta
from .extrapolation_common import get_args, get_data, dump_fit_result


def main():
    args = get_args(beta=True)
    data = get_data(
        args.data_filenames,
        ["ps_mass", "ps_decay_constant", "mAS"],
        beta=args.beta,
    )

    fit_result = meson_beta(data["ps_mass"], data["ps_decay_constant"])
    bare_mass_range = f"[{data['mAS'].min()}, {data['mAS'].max()}]"

    dump_fit_result(args, fit_result, ["A", "B"], bare_mass_range=bare_mass_range)


if __name__ == "__main__":
    main()
