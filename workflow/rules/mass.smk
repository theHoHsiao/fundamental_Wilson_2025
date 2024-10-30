import pandas as pd
from functools import partial


all_metadata = pd.read_csv("metadata/ensemble_metadata.csv")

metadata_query = "Nc == {Nc} & Nt == {Nt} & Ns == {Ns} & beta == {beta} & nAS == {nAS} & mAS == {mAS}"
metadata_lookup = partial(lookup, within=all_metadata, query=metadata_query)

dir_template = "Sp{Nc}b{beta}nAS{nAS}mAS{mAS}T{Nt}L{Ns}"

channels = ["ps", "v", "t", "av", "at", "s"]

rule fit_mass_wall:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        metadata=metadata_lookup(),
        plateau_start=lambda wildcards: metadata_lookup(cols=f"{wildcards.channel}_plateau_start"),
        plateau_end=lambda wildcards: metadata_lookup(cols=f"{wildcards.channel}_plateau_end"),
    input:
        data="data_assets/correlators_wall.h5",
        script="src/mass_wall.py",
    output:
        mean=f"intermediary_data/{dir_template}/meson_{{channel}}_mean.csv",
        samples=f"intermediary_data/{dir_template}/meson_{{channel}}_samples.json",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file_mean {output.mean} --output_file_samples {output.samples} --ensemble_name {params.metadata.ensemble_name} --beta {params.metadata.beta} --mAS {params.metadata.mAS} --Nt {params.metadata.Nt} --Ns {params.metadata.Ns} --channel {wildcards.channel} --plateau_start {params.plateau_start} --plateau_end {params.plateau_end} --min_trajectory {params.metadata.init_conf} --max_trajectory {params.metadata.final_conf} --trajectory_step {params.metadata.delta_conf_spectrum}"

rule fit_mass_smear:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        metadata=metadata_lookup(),
        plateau_start=lambda wildcards: metadata_lookup(cols=f"smear_{wildcards.channel}_plateau_start"),
        plateau_end=lambda wildcards: metadata_lookup(cols=f"smear_{wildcards.channel}_plateau_end"),
        N_sink=lambda wildcards: metadata_lookup(cols=f"{wildcards.channel}_N_sink"),
    input:
        data="data_assets/correlators_smear.h5",
        script="src/mass_smear.py",
    output:
        mean=f"intermediary_data/{dir_template}/smear_meson_{{channel}}_mean.csv",
        samples=f"intermediary_data/{dir_template}/smear_meson_{{channel}}_samples.json",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file_mean {output.mean} --output_file_samples {output.samples} --ensemble_name {params.metadata.ensemble_name} --beta {params.metadata.beta} --mAS {params.metadata.mAS} --Nt {params.metadata.Nt} --Ns {params.metadata.Ns} --channel {wildcards.channel} --plateau_start {params.plateau_start} --plateau_end {params.plateau_end} --min_trajectory {params.metadata.init_conf} --max_trajectory {params.metadata.final_conf} --trajectory_step {params.metadata.delta_conf_spectrum} --N_sink {params.N_sink} --num_source {params.metadata.smear_num_source} --epsilon {params.metadata.smear_epsilon}"

rule fit_mass_GEVP: # for rhoE1 only
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        metadata=metadata_lookup(),
    input:
        data="data_assets/correlators_smear.h5",
        script="src/mass_gevp.py",
    output:
        mean=f"intermediary_data/{dir_template}/gevp_smear_meson_rhoE1_mean.csv",
        samples=f"intermediary_data/{dir_template}/gevp_smear_meson_rhoE1_samples.json",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file_mean {output.mean} --output_file_samples {output.samples} --ensemble_name {params.metadata.ensemble_name} --beta {params.metadata.beta} --mAS {params.metadata.mAS} --Nt {params.metadata.Nt} --Ns {params.metadata.Ns} --plateau_start {params.metadata.smear_rhoE1_plateau_start} --plateau_end {params.metadata.smear_rhoE1_plateau_end} --min_trajectory {params.metadata.init_conf} --max_trajectory {params.metadata.final_conf} --trajectory_step {params.metadata.delta_conf_spectrum} --N_sink {params.metadata.rhoE1_N_sink} --num_source {params.metadata.smear_num_source} --epsilon {params.metadata.smear_epsilon} --GEVP_t0 {params.metadata.GEVP_t0}"

def mass_samples(wildcards):
    return [
        f"intermediary_data/{dir_template}/meson_{wildcards.channel}_samples.json",
        f"intermediary_data/{dir_template}/plaquette_samples.json"
    ]

rule get_meson_decay:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=mass_samples,
        script="src/decay_constant.py",
    output:
        mean=f"intermediary_data/{dir_template}/decay_constant_{{channel}}_mean.csv",
        samples=f"intermediary_data/{dir_template}/decay_constant_{{channel}}_samples.json",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --channel {wildcards.channel} --output_file_mean {output.mean} --output_file_samples {output.samples}"


def ratio_fps_samples(wildcards):
    return [
        f"intermediary_data/{dir_template}/decay_constant_ps_samples.json",
        f"intermediary_data/{dir_template}/meson_{wildcards.channel}_samples.json"
    ]

rule get_fps_ratio:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=ratio_fps_samples,
        script="src/get_Rfps.py",
    output:
        mean=f"intermediary_data/{dir_template}/Rfps_{{channel}}_mean.csv",
        samples=f"intermediary_data/{dir_template}/Rfps_{{channel}}_samples.json",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --channel {wildcards.channel} --output_file_mean {output.mean} --output_file_samples {output.samples}"

def ratio_fps_smear_samples(wildcards):
    return [
        f"intermediary_data/{dir_template}/decay_constant_ps_samples.json",
        f"intermediary_data/{dir_template}/smear_meson_{wildcards.channel}_samples.json",
    ]

def ratio_fps_rhoE1_samples(wildcards):
    return [
        f"intermediary_data/{dir_template}/decay_constant_ps_samples.json",
        f"intermediary_data/{dir_template}/gevp_smear_meson_rhoE1_samples.json",
    ]

rule get_fps_smear_meson_ratio:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=ratio_fps_smear_samples,
        script="src/get_Rfps.py",
    output:
        mean=f"intermediary_data/{dir_template}/smear_Rfps_{{channel}}_mean.csv",
        samples=f"intermediary_data/{dir_template}/smear_Rfps_{{channel}}_samples.json",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --channel {wildcards.channel} --output_file_mean {output.mean} --output_file_samples {output.samples} --smear True"

rule get_fps_smear_rhoE1_ratio:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=ratio_fps_rhoE1_samples,
        script="src/get_Rfps.py",
    output:
        mean=f"intermediary_data/{dir_template}/gevp_smear_Rfps_rhoE1_mean.csv",
        samples=f"intermediary_data/{dir_template}/gevp_smear_Rfps_rhoE1_samples.json",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --channel rhoE1 --output_file_mean {output.mean} --output_file_samples {output.samples} --smear True"





def linear_fittable_pcac_data(wildcards):
    return [
        f"intermediary_data/{dir_template}/mass_mean.csv".format(**row)
        for row in metadata.to_dict(orient="records")
        if row["use_in_linear_PCAC_fits"] and row["use_in_main_plots"]
    ]
