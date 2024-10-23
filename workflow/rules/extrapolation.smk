from functools import partial
import pandas as pd


metadata = pd.read_csv("metadata/ensemble_metadata.csv")


def all_samples(wildcards, observables):
    return [
        f"intermediary_data/{dir_template}/{observable}_samples.json".format(**row)
        for observable in observables
        for row in metadata.to_dict(orient="records")
        if row["use_in_main_plots"]
        if row["use_in_extrapolation"]
    ]


rule Mass_continuum_massless_extrapolation:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=partial(all_samples, observables=["w0", "smear_meson_ps", "smear_meson_v", "smear_meson_t","smear_meson_av","smear_meson_at","smear_meson_s", "gevp_smear_meson_rhoE1"]),
        script="src/extrapolation.py",
    output:
        mean=f"intermediary_data/extrapolation_results/{{channel}}_extp_mass_mean.csv",
        samples=f"intermediary_data/extrapolation_results/{{channel}}_extp_mass_samples.json",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file_mean {output.mean} --output_file_samples {output.samples} --channel {wildcards.channel}"

rule Decay_continuum_massless_extrapolation:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=partial(all_samples, observables=["w0", "meson_ps", "meson_v", "meson_av", "plaquette"]),
        script="src/extrapolation_decay.py",
    output:
        mean=f"intermediary_data/extrapolation_results/{{channel}}_extp_deacy_mean.csv",
        samples=f"intermediary_data/extrapolation_results/{{channel}}_extp_decay_samples.json",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file_mean {output.mean} --output_file_samples {output.samples} --channel {wildcards.channel}"
