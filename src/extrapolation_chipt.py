#!/usr/bin/env python3

from .fitting import global_meson_fit, linear_fit_form
from .extrapolation_common import get_args, get_data, dump_fit_result


def main():
    args = get_args(beta=True)
    data = get_data(
        args.data_filenames,
        ["ps_mass", "ps_decay_constant", "mAS"],
        beta=args.beta,
    )

    fit_result = global_meson_fit(
        linear_fit_form,
        data["ps_mass_squared"],
        data["ps_decay_constant_squared"],
    )

    dump_fit_result(
        args,
        "chipt",
        fit_result,
        ["A", "B"],
        mAS_min=data["mAS"].min(),
        mAS_max=data["mAS"].max(),
    )


if __name__ == "__main__":
    main()
