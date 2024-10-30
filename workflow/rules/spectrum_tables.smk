import pandas as pd
from functools import partial


all_metadata = pd.read_csv("metadata/ensemble_metadata.csv")

metadata_query = "Nc == {Nc} & Nt == {Nt} & Ns == {Ns} & beta == {beta} & nAS == {nAS} & mAS == {mAS}"
metadata_lookup = partial(lookup, within=all_metadata, query=metadata_query)

dir_template = "Sp{Nc}b{beta}nAS{nAS}mAS{mAS}T{Nt}L{Ns}"

channels = ["ps", "v", "t", "av", "at", "s"]

def mpcac_data(wildcards):
    return [
        f"intermediary_data/{dir_template}/mpcac_mean.csv".format(**row)
        for row in metadata.to_dict(orient="records")
        if row["use_in_main_plots"]
    ]

def wall_mass_data(wildcards):
    return [
        f"intermediary_data/{dir_template}/meson_{channel}_mean.csv".format(channel=channel, **row)
        for row in metadata.to_dict(orient="records")
        for channel in channels
        if row["use_in_main_plots"]
    ]

def smear_mass_data(wildcards):
    return [
        f"intermediary_data/{dir_template}/smear_meson_{channel}_mean.csv".format(channel=channel, **row)
        for row in metadata.to_dict(orient="records")
        for channel in channels
        if row["use_in_main_plots"]
        if row["use_smear"]
    ] + [
        f"intermediary_data/{dir_template}/gevp_smear_meson_rhoE1_mean.csv".format(**row)
        for row in metadata.to_dict(orient="records")
        if row["use_in_main_plots"]
        if row["use_smear"]
    ]

def decay_constant_data(wildcards):
    return [
        f"intermediary_data/{dir_template}/decay_constant_{channel}_mean.csv".format(channel=channel, **row)
        for row in metadata.to_dict(orient="records")
        for channel in ["ps", "v", "av"]
        if row["use_in_main_plots"]
    ]

def smear_Rfps_data(wildcards):
    return [
        f"intermediary_data/{dir_template}/smear_Rfps_{channel}_mean.csv".format(channel=channel, **row)
        for row in metadata.to_dict(orient="records")
        for channel in channels
        if row["use_in_main_plots"]
        if row["use_smear"]
    ] + [
        f"intermediary_data/{dir_template}/gevp_smear_Rfps_rhoE1_mean.csv".format( **row)
        for row in metadata.to_dict(orient="records")
        if row["use_in_main_plots"]
        if row["use_smear"]
    ]

def smear_Rmv_data(wildcards):
    return [
        f"intermediary_data/{dir_template}/gevp_smear_Rmv_rhoE1_mean.csv".format( **row)
        for row in metadata.to_dict(orient="records")
        if row["use_in_main_plots"]
        if row["use_smear"]
    ]

rule wall_mass_table:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        mass_data=wall_mass_data,
        mpcac_data=mpcac_data,
        decay_data=decay_constant_data,
        metadata_csv="metadata/ensemble_metadata.csv",
        script="src/tables/wall_mass_table.py",
    output:
        table="assets/tables/table_VI.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.mass_data} {input.mpcac_data} {input.decay_data} --output_file {output.table}"

rule smear_mass_table:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=smear_mass_data,
        script="src/tables/smear_mass_table.py",
    output:
        table="assets/tables/table_IX.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.table}"

rule smear_mass_fps_table:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=smear_Rfps_data,
        data_Rmv=smear_Rmv_data,
        script="src/tables/smear_mass_fps_table.py",
    output:
        table="assets/tables/table_X.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} {input.data_Rmv} --output_file {output.table}"
