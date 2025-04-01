import pandas as pd
from functools import partial


all_metadata = pd.read_csv("metadata/ensemble_metadata.csv")

metadata_query = "Nc == {Nc} & Nt == {Nt} & Ns == {Ns} & beta == {beta} & nF == {nF} & mF == {mF}"
metadata_lookup = partial(lookup, within=all_metadata, query=metadata_query)

dir_template = "Sp{Nc}b{beta}nF{nF}mF{mF}T{Nt}L{Ns}"

channels = ["ps", "v"]


rule gevp_meson_mass:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        metadata=metadata_lookup(),
        E0_plateau_start=lambda wildcards: metadata_lookup(cols=f"{wildcards.channel}_E0_plateau_start"),
        E0_plateau_end=lambda wildcards: metadata_lookup(cols=f"{wildcards.channel}_E0_plateau_end"),
        E1_plateau_start=lambda wildcards: metadata_lookup(cols=f"{wildcards.channel}_E1_plateau_start"),
        E1_plateau_end=lambda wildcards: metadata_lookup(cols=f"{wildcards.channel}_E1_plateau_end"),
        E2_plateau_start=lambda wildcards: metadata_lookup(cols=f"{wildcards.channel}_E2_plateau_start"),
        E2_plateau_end=lambda wildcards: metadata_lookup(cols=f"{wildcards.channel}_E2_plateau_end"),
        E3_plateau_start=lambda wildcards: metadata_lookup(cols=f"{wildcards.channel}_E3_plateau_start"),
        E3_plateau_end=lambda wildcards: metadata_lookup(cols=f"{wildcards.channel}_E3_plateau_end"),
    input:
        data="data_assets/corr_sp4_FUN.h5",
        script="src/mass_gevp_meson.py",
    output:
        samples=f"intermediary_data/{dir_template}/meson_gevp_{{channel}}_samples.json",
        plot=f"intermediary_data/{dir_template}/meson_gevp_{{channel}}_eff_mass.pdf",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file_samples {output.samples} --ensemble_name {params.metadata.ensemble_name} --effmass_plot_file {output.plot} --plot_styles {plot_styles}"
        " --beta {params.metadata.beta} --mF {params.metadata.mF} --Nt {params.metadata.Nt} --Ns {params.metadata.Ns}"
        " --min_trajectory {params.metadata.init_conf} --max_trajectory {params.metadata.final_conf} --trajectory_step {params.metadata.delta_conf}"
        " --channel {wildcards.channel} --gevp_t0 {params.metadata.gevp_t0}"
        " --n_smear_min {params.metadata.n_smear_min} --n_smear_max {params.metadata.n_smear_max} --n_smear_diff {params.metadata.n_smear_diff}"
        " --E0_plateau_start {params.E0_plateau_start} --E0_plateau_end {params.E0_plateau_end}"
        " --E1_plateau_start {params.E1_plateau_start} --E1_plateau_end {params.E1_plateau_end}"
        " --E2_plateau_start {params.E2_plateau_start} --E2_plateau_end {params.E2_plateau_end}"
        " --E3_plateau_start {params.E3_plateau_start} --E3_plateau_end {params.E3_plateau_end}"


rule meson_matrix_element:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        metadata=metadata_lookup(),
        plateau_start=lambda wildcards: metadata_lookup(cols=f"{wildcards.channel}_matrix_element_plateau_start"),
        plateau_end=lambda wildcards: metadata_lookup(cols=f"{wildcards.channel}_matrix_element_plateau_end"),
    input:
        data="data_assets/corr_sp4_FUN.h5",
        script="src/matrix_element_meson.py",
    output:
        mean=f"intermediary_data/{dir_template}/meson_extraction_{{channel}}_mean.csv",
        samples=f"intermediary_data/{dir_template}/meson_extraction_{{channel}}_samples.json",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file_mean {output.mean} --output_file_samples {output.samples} --ensemble_name {params.metadata.ensemble_name}"
        " --beta {params.metadata.beta} --mF {params.metadata.mF} --Nt {params.metadata.Nt} --Ns {params.metadata.Ns}"
        " --min_trajectory {params.metadata.init_conf} --max_trajectory {params.metadata.final_conf} --trajectory_step {params.metadata.delta_conf}"
        " --n_smear_max {params.metadata.n_smear_max} --channel {wildcards.channel} --E0_plateau_start {params.plateau_start} --E0_plateau_end {params.plateau_end}"
 

def extraction_samples(wildcards):
    return [
        f"intermediary_data/{dir_template}/meson_extraction_{rep}_{channel}_samples.json".format(**row)
        for row in metadata.to_dict(orient="records")
        for channel in ["ps", "v"]
        for rep in ["f"]
    ] 


rule get_meson_decay:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        mass=f"intermediary_data/{dir_template}/meson_extraction_{{channel}}_samples.json",
        plaq=f"intermediary_data/{dir_template}/plaquette_samples.json",
        script="src/decay_constant.py",
    output:
        mean=f"intermediary_data/{dir_template}/decay_constant_{{channel}}_mean.csv",
        samples=f"intermediary_data/{dir_template}/decay_constant_{{channel}}_samples.json",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.mass} {input.plaq} --channel {wildcards.channel} --output_file_mean {output.mean} --output_file_samples {output.samples}"



rule ps_correlator_autocorrelation:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        metadata=metadata_lookup(),
        plateau_start=lambda wildcards: metadata_lookup(cols="ps_plateau_start"),
        plateau_end=lambda wildcards: metadata_lookup(cols=f"ps_plateau_end"),
    input:
        data="data_assets/corr_sp4_FUN.h5",
        script="src/ps_correlators_autocorrelation.py",
    output:
        mean=f"intermediary_data/{dir_template}/tau_ps_correlator_mean.csv",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file_mean {output.mean} --ensemble_name {params.metadata.ensemble_name} --beta {params.metadata.beta} --mAS {params.metadata.mAS} --Nt {params.metadata.Nt} --Ns {params.metadata.Ns} --plateau_start {params.plateau_start} --plateau_end {params.plateau_end} --min_trajectory {params.metadata.init_conf} --max_trajectory {params.metadata.final_conf} --trajectory_step {params.metadata.delta_conf_spectrum}"


def mass_samples(wildcards):
    return [
        f"intermediary_data/{dir_template}/meson_{wildcards.channel}_samples.json",
        f"intermediary_data/{dir_template}/plaquette_samples.json",
    ]


def ratio_fps_samples(wildcards):
    return [
        f"intermediary_data/{dir_template}/decay_constant_ps_samples.json",
        f"intermediary_data/{dir_template}/meson_{wildcards.channel}_samples.json",
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


def ratio_mv_rhoE1_samples(wildcards):
    return [
        f"intermediary_data/{dir_template}/smear_meson_v_samples.json",
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


rule get_Rmv_smear_rhoE1_ratio:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=ratio_mv_rhoE1_samples,
        script="src/get_Rmv.py",
    output:
        mean=f"intermediary_data/{dir_template}/gevp_smear_Rmv_rhoE1_mean.csv",
        samples=f"intermediary_data/{dir_template}/gevp_smear_Rmv_rhoE1_samples.json",
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
