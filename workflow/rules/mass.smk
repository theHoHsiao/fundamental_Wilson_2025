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

def all_mass_data(wildcards):
    return [
        f"intermediary_data/{dir_template}/meson_{channel}_mean.csv".format(channel=channel, **row)
        for row in metadata.to_dict(orient="records")
        for channel in channels
        if row["use_in_main_plots"]
    ]

def mpcac_data(wildcards):
    return [
        f"intermediary_data/{dir_template}/mpcac_mean.csv".format(**row)
        for row in metadata.to_dict(orient="records")
        if row["use_in_main_plots"]
    ]

def smear_mass_data(wildcards):
    return [
        f"intermediary_data/{dir_template}/smear_meson_{channel}_mean.csv".format(channel=channel, **row)
        for row in metadata.to_dict(orient="records")
        for channel in channels
        if row["use_in_main_plots"]
        if row["use_smear"]
    ]

def gevp_mass_data(wildcards):
    return [
        f"intermediary_data/{dir_template}/gevp_smear_meson_rhoE1_mean.csv".format(**row)
        for row in metadata.to_dict(orient="records")
        if row["use_in_main_plots"]
        if row["use_smear"]
    ]



def linear_fittable_pcac_data(wildcards):
    return [
        f"intermediary_data/{dir_template}/mass_mean.csv".format(**row)
        for row in metadata.to_dict(orient="records")
        if row["use_in_linear_PCAC_fits"] and row["use_in_main_plots"]
    ]

rule mass_table:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        mass_data=all_mass_data,
        mpcac_data=mpcac_data,
        metadata_csv="metadata/ensemble_metadata.csv",
        script="src/tables/mass_table.py",
    output:
        table="assets/tables/mass_table.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.mass_data} {input.mpcac_data} --output_file {output.table}"

rule smear_mass_table:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=smear_mass_data,
        data2=gevp_mass_data,
        script="src/tables/mass_table.py",
    output:
        table="assets/tables/smear_mass_table.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.table}"
