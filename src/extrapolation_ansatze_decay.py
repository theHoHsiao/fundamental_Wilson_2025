#!/usr/bin/env python3

from functools import partial

from .fitting import global_meson_fit, split_means_samples
from .extrapolation_common import dump_fit_result_chisquare, get_args, get_data, ansatz_form
from .extrapolation_common import initialize_fit_parameters_from_a

def ansatz_parameter_names(fit_prefix):
    return {
        "a" : ["F_a", "Lf_a", "Wf_a"],
        "a2": ["F_a2", "Lf_a2", "Wf_a2", "Rf_a2"],
        "m4": ["F_m4", "Lf_m4", "Wf_m4", "Pf_m4"],
        "a2_m4": ["F_a2_m4", "Lf_a2_m4", "Wf_a2_m4", "Pf_a2_m4", "Rf_a2_m4"],
        "a2_am2": ["F_a2_am2", "Lf_a2_am2", "Wf_a2_am2", "Rf_a2_am2", "Cf_a2_am2"],
        "full": ["F_full", "Lf_full", "Wf_full", "Pf_full", "Rf_full", "Cf_full"],
    }[fit_prefix]


def main():
    args = get_args(channels=["ps", "v", "av"], ansatz=["a", "a2", "m4", "a2_m4", "a2_am2", "full"])
    channel_obs_key = f"f_{args.channel}_decay_constant"
    mass_channel_obs_key = f"gevp_f_{args.channel}_E0_mass"
    data = get_data(args.data_filenames, ["w0", "gevp_f_ps_E0_mass", channel_obs_key], cutoff=1, cutoff_observable=mass_channel_obs_key)
    #print(data.keys())

    lat_a_means, _ = split_means_samples(data["lat_a"])

    fit_form = ansatz_form(args.ansatz)
    parameter_names = ansatz_parameter_names(args.ansatz)
    print(f"Using ansatz {args.ansatz} with fit form == {fit_form.__name__} == and parameters {parameter_names}")

    init_params = initialize_fit_parameters_from_a(fit_form, data, lat_a_means, channel_obs_key, [0.01, 2.0, -0.01])
    
    #init_params = initialize_fit_parameters(fit_form, [0.01, 2.0, -0.01])

    fit_result = global_meson_fit(
        partial(fit_form, lat_a=lat_a_means),
        data["gevp_f_ps_E0_mass_hat_squared"],
        data[f"{channel_obs_key}_hat_squared"],
        initial_guess=init_params,
    )

    dump_fit_result_chisquare(args, "continuum_Wilson_ChiPT", fit_result, parameter_names)


if __name__ == "__main__":
    main()
