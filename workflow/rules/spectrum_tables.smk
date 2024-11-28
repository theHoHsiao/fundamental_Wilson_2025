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

def tau_ps_correlator_data(wildcards):
    return [
        f"intermediary_data/{dir_template}/tau_ps_correlator_mean.csv".format(**row)
        for row in metadata.to_dict(orient="records")
        if row["use_in_main_plots"]
    ]


def wall_mass_data(wildcards):
    return [
        f"intermediary_data/{dir_template}/meson_{channel}_mean.csv".format(
            channel=channel, **row
        )
        for row in metadata.to_dict(orient="records")
        for channel in channels
        if row["use_in_main_plots"]
    ]


def smear_mass_data(wildcards):
    return [
        f"intermediary_data/{dir_template}/smear_meson_{channel}_mean.csv".format(
            channel=channel, **row
        )
        for row in metadata.to_dict(orient="records")
        for channel in channels
        if row["use_in_main_plots"]
        if row["use_smear"]
    ] + [
        f"intermediary_data/{dir_template}/gevp_smear_meson_rhoE1_mean.csv".format(
            **row
        )
        for row in metadata.to_dict(orient="records")
        if row["use_in_main_plots"]
        if row["use_smear"]
    ]


def decay_constant_data(wildcards):
    return [
        f"intermediary_data/{dir_template}/decay_constant_{channel}_mean.csv".format(
            channel=channel, **row
        )
        for row in metadata.to_dict(orient="records")
        for channel in ["ps", "v", "av"]
        if row["use_in_main_plots"]
    ]


def wall_Rfps_data(wildcards):
    return [
        f"intermediary_data/{dir_template}/Rfps_{channel}_mean.csv".format(
            channel=channel, **row
        )
        for row in metadata.to_dict(orient="records")
        for channel in channels
        if row["use_in_main_plots"]
    ]


def smear_Rfps_data(wildcards):
    return [
        f"intermediary_data/{dir_template}/smear_Rfps_{channel}_mean.csv".format(
            channel=channel, **row
        )
        for row in metadata.to_dict(orient="records")
        for channel in channels
        if row["use_in_main_plots"]
        if row["use_smear"]
    ] + [
        f"intermediary_data/{dir_template}/gevp_smear_Rfps_rhoE1_mean.csv".format(**row)
        for row in metadata.to_dict(orient="records")
        if row["use_in_main_plots"]
        if row["use_smear"]
    ]


def smear_Rmv_data(wildcards):
    return [
        f"intermediary_data/{dir_template}/gevp_smear_Rmv_rhoE1_mean.csv".format(**row)
        for row in metadata.to_dict(orient="records")
        if row["use_in_main_plots"]
        if row["use_smear"]
    ]


def continuum_massless_extrapolation_mass(wildcards):
    return [
        f"intermediary_data/extrapolation_results/{channel}_extp_mass_mean.csv".format()
        for channel in ["v", "t", "av", "at", "s", "rhoE1"]
    ]


def continuum_massless_extrapolation_decay(wildcards):
    return [
        f"intermediary_data/extrapolation_results/{channel}_extp_decay_mean.csv".format()
        for channel in ["ps", "v", "av"]
    ]


def chipt_extrapolation_results(wildcards):
    return [
        f"intermediary_data/chipt_extrapolation_results/chipt_b{beta}_extp_mean.csv".format()
        for beta in [6.6, 6.65, 6.7, 6.75, 6.8]
    ]


def deft_extrapolation_results(wildcards):
    return [
        f"intermediary_data/deft_extrapolation_results/deft_b{beta}_extp_mean.csv".format()
        for beta in [6.6, 6.65, 6.7, 6.75, 6.8]
    ]


def autocorr_data(wildcards):
    return [
        f"intermediary_data/{dir_template}/{basename}.csv".format(**row)
        for basename in [
            "tau_ps_correlator_mean",
            "plaquette_mean",
            "w0_mean",
            "top_charge",
        ]
        for row in metadata.to_dict(orient="records")
        if row["use_in_main_plots"]
    ]


def finite_volume_data(wildcards):
    return [
        f"intermediary_data/{dir_template}/{basename}.csv".format(**row)
        for basename in [
            "plaquette_mean",
            "meson_ps_mean",
            "meson_v_mean",
            "mpcac_mean",
        ]
        for row in metadata.to_dict(orient="records")
        if row["use_in_finite_volume"]
    ]


rule continuum_massless_decay:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=continuum_massless_extrapolation_decay,
        script="src/tables/continuum_massless_decay.py",
    output:
        table="assets/tables/nlo_coefficients_decay.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.table}"


rule continuum_massless_mass:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=continuum_massless_extrapolation_mass,
        script="src/tables/continuum_massless_mass.py",
    output:
        table="assets/tables/nlo_coefficients_mass.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.table}"


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
        table="assets/tables/wall_results1.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.mass_data} {input.mpcac_data} {input.decay_data} --output_file {output.table}" #{input.decay_data}


rule wall_mass_table2:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        mass_data=wall_mass_data,
        mpcac_data=mpcac_data,
        decay_data=decay_constant_data,
        metadata_csv="metadata/ensemble_metadata.csv",
        script="src/tables/wall_mass_table2.py",
    output:
        table="assets/tables/wall_results2.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.mass_data} {input.mpcac_data} {input.decay_data} --output_file {output.table}"


rule wall_mass_fps_table:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=wall_Rfps_data,
        script="src/tables/wall_mass_fps_table.py",
    output:
        table="assets/tables/wall_results_ratio.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.table}"


rule smear_mass_table:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=smear_mass_data,
        script="src/tables/smear_mass_table.py",
    output:
        table="assets/tables/smear_mass.tex",
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
        table="assets/tables/smear_decay.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} {input.data_Rmv} --output_file {output.table}"


rule chipt_table:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=chipt_extrapolation_results,
        script="src/tables/chipt_table.py",
    output:
        table="assets/tables/chipt_table.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.table}"


rule deft_table:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=deft_extrapolation_results,
        script="src/tables/deft_table.py",
    output:
        table="assets/tables/deft_table.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.table}"


rule autocorr_table:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=autocorr_data,
        script="src/tables/autocorr.py",
    output:
        table="assets/tables/autocorr_table.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.table}"


rule finite_volume_table:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=finite_volume_data,
        script="src/tables/finite_volume.py",
    output:
        table="assets/tables/finite_volume_table.tex",
        definitions="assets/definitions/finite_volume_definitions.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.table} --definitions_file {output.definitions}"
