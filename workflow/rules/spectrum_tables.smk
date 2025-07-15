import pandas as pd
from functools import partial


all_metadata = pd.read_csv("metadata/ensemble_metadata.csv")

metadata_query = "Nc == {Nc} & Nt == {Nt} & Ns == {Ns} & beta == {beta} & nF == {nF} & mF == {mF}"
metadata_lookup = partial(lookup, within=all_metadata, query=metadata_query)

dir_template = "Sp{Nc}b{beta}nF{nF}mF{mF}T{Nt}L{Ns}"

def extraction_means(wildcards):
    return [
        f"intermediary_data/{dir_template}/meson_extraction_{rep}_{channel}_mean.csv".format(**row)
        for row in metadata.to_dict(orient="records")
        for channel in ["ps", "v"]
        for rep in ["f"]
        if row["use_in_table"]
    ]


def gevp_E0_means(wildcards):
    return [
        f"intermediary_data/{dir_template}/meson_gevp_E0_{rep}_{channel}_mean.csv".format(**row)
        for row in metadata.to_dict(orient="records")
        for channel in ["ps", "v"]
        for rep in ["f"]
        if row["use_in_table"]
    ] 


def decay_constant_means(wildcards):
    return [
        f"intermediary_data/{dir_template}/decay_constant_{rep}_{channel}_mean.csv".format(**row)
        for row in metadata.to_dict(orient="records")
        for channel in ["ps", "v", "av"]
        for rep in ["f"]
        if row["use_in_table"]
    ] 

def w0_means(wildcards):
    return [
        f"intermediary_data/{dir_template}/w0_mean.csv".format(**row)
        for row in metadata.to_dict(orient="records")
        if row["use_in_w0"]
    ] 

def plaq_means(wildcards):
    return [
        f"intermediary_data/{dir_template}/plaquette_mean.csv".format(**row)
        for row in metadata.to_dict(orient="records")
        if row["use_in_table"]
    ] 


def continuum_massless_extrapolation_mass(wildcards):
    return [
        f"intermediary_data/extrapolation_results/f_{channel}_extp_mass_mean.csv".format()
        for channel in ["v", "t", "av", "at", "s"]
    ]


def continuum_massless_extrapolation_decayconstant(wildcards):
    return [
        f"intermediary_data/extrapolation_results/f_{channel}_extp_decayconstant_mean.csv".format()
        for channel in ["ps", "v", "av"]
    ]


rule f_meson_mass_table:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        decay=decay_constant_means,
        mass=extraction_means,
        E0=gevp_E0_means,
        w0=w0_means,
        script="src/tables/print_spectrum_meson_f.py",
    output:
        table="assets/tables/spectrum_meson_f.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.decay} {input.mass} {input.E0} {input.w0} --output_file {output.table}"


rule ens_table:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        plaq=plaq_means,
        w0=w0_means,
        script="src/tables/print_ens.py",
    output:
        table="assets/tables/ensemble.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.plaq} {input.w0} --output_file {output.table}"


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


rule continuum_massless_decayconstant:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=continuum_massless_extrapolation_decayconstant,
        script="src/tables/continuum_massless_decayconstant.py",
    output:
        table="assets/tables/nlo_coefficients_decayconstant.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.table}"
