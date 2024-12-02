def expand_ensemble(observables, condition_key, dir_slugs=[""]):
    return [
        f"intermediary_data/{dir_template}{dir_slug}/{observable}_mean.csv".format(
            **row
        )
        for observable in observables
        for row in metadata.to_dict(orient="records")
        for dir_slug in dir_slugs
        if row[condition_key]
    ]


def all_spectrum_csvs(wildcards):
    all_channels = ["at", "av", "ps", "s", "t", "v"]
    f_channels = ["av", "ps", "v"]
    main_plot_observables = (
        [f"Rfps_{channel}" for channel in all_channels]
        + [f"decay_constant_{channel}" for channel in f_channels]
        + [f"meson_{channel}" for channel in all_channels]
        + ["mpcac", "plaquette", "top_charge", "tau_ps_correlator", "w0"]
    )
    smear_observables = (
        [f"smear_Rfps_{channel}" for channel in all_channels]
        + [f"smear_meson_{channel}" for channel in all_channels]
        + [f"gevp_smear_{obs}_rhoE1" for obs in ["Rfps", "Rmv", "meson"]]
    )
    finite_volume_observables = ["meson_ps", "meson_v", "mpcac", "plaquette"]
    files_with_duplicates = (
        expand_ensemble(main_plot_observables, "use_in_main_plots")
        + expand_ensemble(smear_observables, "use_smear")
        + expand_ensemble(finite_volume_observables, "use_in_finite_volume")
        + expand_ensemble(
            ["plaquette"],
            "use_in_plaquette_phase_diagram",
            ["_unitstart", "_randomstart"],
        )
    )
    return list(dict.fromkeys(files_with_duplicates))


rule ensemble_csv:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=all_spectrum_csvs,
        script="src/csvs/spectrum_csv.py",
    output:
        csv="data_assets/ensemble_data.csv",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.csv}"


def per_beta_results(wildcards):
    return [
        f"intermediary_data/{fit_form}_extrapolation_results/{fit_form}_b{beta}_extp_mean.csv"
        for fit_form in ["deft", "chipt"]
        for beta in [6.6, 6.65, 6.7, 6.75, 6.8]
    ]


rule per_beta_fits:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=per_beta_results,
        script="src/csvs/per_beta_csv.py",
    output:
        csv="data_assets/chipt_deft_data.csv",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.csv}"


def global_eft_results(wildcards):
    channels = {
        "mass": ("at", "av", "rhoE1", "s", "t", "v"),
        "decayconstant": ("av", "ps", "v"),
    }
    return [
        f"intermediary_data/extrapolation_results/{channel}_extp_{observable}_mean.csv"
        for observable, obs_channels in channels.items()
        for channel in obs_channels
    ] + ["intermediary_data/extrapolation_results/R_mvdfps_extp_mean.csv"]


rule global_eft_fits:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=global_eft_results,
        script="src/csvs/global_eft_csv.py",
    output:
        csv="data_assets/global_eft_data.csv",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.csv}"
