#!/usr/bin/env python3


from .fitting import global_meson_fit, linear_fit_form
from .extrapolation_common import get_args, get_data, dump_fit_result


def main():
    args = get_args(beta=True)
    data = get_data(
        args.data_filenames,
        ["ps_mass", "ps_decay_constant", "mPCAC", "mAS"],
        beta=args.beta,
    )

    fit_result = global_meson_fit(
        linear_fit_form,
        data["log_ps_decay_constant_squared"],
        data["log_ps_mass_squared_over_mPCAC"],
    )

    dump_fit_result(
        args,
        "deft",
        fit_result,
        ["A", "B"],
        y=fit_result[0][1] * 2 + 2,
        mAS_min=data["mAS"].min(),
        mAS_max=data["mAS"].max(),
    )


if __name__ == "__main__":
    main()
