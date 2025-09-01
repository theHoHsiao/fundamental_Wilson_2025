from functools import partial
import pandas as pd


metadata = pd.read_csv("metadata/ensemble_metadata.csv")


def extp_samples(wildcards, observables):
    return [
        f"intermediary_data/{dir_template}/{observable}_samples.json".format(**row)
        for observable in observables
        for row in metadata.to_dict(orient="records")
        #if row["use_in_main_plots"]
        if row["use_in_extrapolation"]
    ]


rule Mass_continuum_massless_extrapolation:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=partial(
            extp_samples,
            observables=[
                "w0",
                "meson_gevp_f_ps",
                "meson_gevp_f_v",
                "meson_gevp_f_t",
                "meson_gevp_f_av",
                "meson_gevp_f_at",
                "meson_gevp_f_s",
            ],
        ),
        script="src/extrapolation_mass.py",
    output:
        mean=f"intermediary_data/extrapolation_results/f_{{channel}}_extp_mass_mean.csv",
        samples=f"intermediary_data/extrapolation_results/f_{{channel}}_extp_mass_samples.json",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file_mean {output.mean} --output_file_samples {output.samples} --channel {wildcards.channel}"


rule Decay_continuum_massless_extrapolation:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=partial(
            extp_samples,
            observables=[
                "w0",
                "meson_gevp_f_ps",
                "decay_constant_f_ps",
                "decay_constant_f_v",
                "decay_constant_f_av",
            ],
        ),
        script="src/extrapolation_decay.py",
    output:
        mean=f"intermediary_data/extrapolation_results/f_{{channel}}_extp_decayconstant_mean.csv", 
        samples=f"intermediary_data/extrapolation_results/f_{{channel}}_extp_decayconstant_samples.json",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file_mean {output.mean} --output_file_samples {output.samples} --channel {wildcards.channel}"
