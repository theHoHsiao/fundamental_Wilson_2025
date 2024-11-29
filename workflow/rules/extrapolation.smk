from functools import partial
import pandas as pd


metadata = pd.read_csv("metadata/ensemble_metadata.csv")


def extp_samples(wildcards, observables):
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
        data=partial(
            extp_samples,
            observables=[
                "w0",
                "smear_meson_ps",
                "smear_meson_v",
                "smear_meson_t",
                "smear_meson_av",
                "smear_meson_at",
                "smear_meson_s",
                "gevp_smear_meson_rhoE1",
            ],
        ),
        script="src/extrapolation_mass.py",
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
        data=partial(
            extp_samples,
            observables=[
                "w0",
                "meson_ps",
                "decay_constant_ps",
                "decay_constant_v",
                "decay_constant_av",
            ],
        ),
        script="src/extrapolation_decay.py",
    output:
        mean=f"intermediary_data/extrapolation_results/{{channel}}_extp_decay_mean.csv",
        samples=f"intermediary_data/extrapolation_results/{{channel}}_extp_decay_samples.json",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file_mean {output.mean} --output_file_samples {output.samples} --channel {wildcards.channel}"


rule Ratio_continuum_massless_extrapolation:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=partial(
            extp_samples,
            observables=["w0", "meson_ps", "meson_v", "decay_constant_ps"],
        ),
        script="src/extrapolation_ratio.py",
    output:
        mean=f"intermediary_data/extrapolation_results/R_m{{channel}}dfps_extp_mean.csv",
        samples=f"intermediary_data/extrapolation_results/R_m{{channel}}dfps_extp_samples.json",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file_mean {output.mean} --output_file_samples {output.samples} --channel {wildcards.channel}"


rule Chipt_extrapolation:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=partial(all_samples, observables=["meson_ps", "decay_constant_ps"]),
        script="src/extrapolation_chipt.py",
    output:
        mean=f"intermediary_data/chipt_extrapolation_results/chipt_b{{beta}}_extp_mean.csv",
        samples=f"intermediary_data/chipt_extrapolation_results/chipt_b{{beta}}_extp_samples.json",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file_mean {output.mean} --output_file_samples {output.samples} --beta {wildcards.beta}"


rule DEFT_extrapolation:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=partial(
            all_samples, observables=["meson_ps", "decay_constant_ps", "mpcac"]
        ),
        script="src/extrapolation_deft.py",
    output:
        mean=f"intermediary_data/deft_extrapolation_results/deft_b{{beta}}_extp_mean.csv",
        samples=f"intermediary_data/deft_extrapolation_results/deft_b{{beta}}_extp_samples.json",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file_mean {output.mean} --output_file_samples {output.samples} --beta {wildcards.beta}"


rule Chipt_ratio:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        ps_decay="intermediary_data/extrapolation_results/ps_extp_decay_samples.json",
        v_mass="intermediary_data/extrapolation_results/v_extp_mass_samples.json",
        script="src/definitions/mvhat_over_fpshat_chiral.py",
    output:
        definitions="assets/definitions/chipt_ratio.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.ps_decay} {input.v_mass} --output_file {output.definitions}"
