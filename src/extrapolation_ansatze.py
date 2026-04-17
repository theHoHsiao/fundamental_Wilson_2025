#!/usr/bin/env python3

from functools import partial

from .fitting import global_meson_fit, split_means_samples
from .extrapolation_common import dump_fit_result_chisquare, get_args, get_data, ansatz_form
from .extrapolation_common import initialize_fit_parameters_from_a


def ansatz_parameter_names(fit_prefix):
    return {
        "a" : ["M_a", "Lm_a", "Wm_a"],
        "a2": ["M_a2", "Lm_a2", "Wm_a2", "Rm_a2"],
        "m4": ["M_m4", "Lm_m4", "Wm_m4", "Pm_m4"],
        "a2_m4": ["M_a2_m4", "Lm_a2_m4", "Wm_a2_m4", "Pm_a2_m4", "Rm_a2_m4"],
        "a2_am2": ["M_a2_am2", "Lm_a2_am2", "Wm_a2_am2", "Rm_a2_am2", "Cm_a2_am2"],
        "full": ["M_full", "Lm_full", "Wm_full", "Pm_full", "Rm_full", "Cm_full"],
    }[fit_prefix]


def main():
    args = get_args(channels=["ps", "v", "t", "av", "at", "s", "rhoE1"], ansatz=["a", "a2", "m4", "a2_m4", "a2_am2", "full"])
    channel_obs_key = f"gevp_f_{args.channel}_E0_mass"
    data = get_data(args.data_filenames, ["w0", "gevp_f_ps_E0_mass", channel_obs_key], cutoff=1, cutoff_observable=channel_obs_key)
    #print(data.keys())

    lat_a_means, _ = split_means_samples(data["lat_a"])

    fit_form = ansatz_form(args.ansatz)
    parameter_names = ansatz_parameter_names(args.ansatz)
    print(f"Using ansatz {args.ansatz} with fit form == {fit_form.__name__} == and parameters {parameter_names}")

    
    init_params = initialize_fit_parameters_from_a(fit_form, data, lat_a_means, channel_obs_key, [0.5, 2.0, -0.2])

    #init_params = initialize_fit_parameters(fit_form, [0.5, 2.0, -0.2])

    fit_result = global_meson_fit(
        partial(fit_form, lat_a=lat_a_means),
        data["gevp_f_ps_E0_mass_hat_squared"],
        data[f"{channel_obs_key}_hat_squared"],
        initial_guess=init_params,
    )

    dump_fit_result_chisquare(args, "continuum_Wilson_ChiPT", fit_result, parameter_names)


if __name__ == "__main__":
    main()
