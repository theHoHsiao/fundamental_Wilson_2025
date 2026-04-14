import pandas as pd
from functools import partial


all_metadata = pd.read_csv("metadata/ensemble_metadata.csv")

metadata_query = "Nc == {Nc} & Nt == {Nt} & Ns == {Ns} & beta == {beta} & nF == {nF} & mF == {mF}"
metadata_lookup = partial(lookup, within=all_metadata, query=metadata_query)

dir_template = "Sp{Nc}b{beta}nF{nF}mF{mF}T{Nt}L{Ns}"

W0_threshold = 0.28125


def autocorr_data(wildcards):
    return [
        f"intermediary_data/{dir_template}/{basename}.csv".format(**row)
        for basename in [
            "tau_ps_correlator_mean",
            "plaquette_mean",
            "w0_mean",
            "top_charge_mean",
            f"meson_gevp_E0_f_ps_mean",      
        ]
        for row in metadata.to_dict(orient="records")
        if row["use_in_extrapolation"]
    ]


def extraction_means(wildcards):
    return [
        f"intermediary_data/{dir_template}/meson_extraction_{rep}_{channel}_mean.csv".format(**row)
        for row in metadata.to_dict(orient="records")
        for channel in ["ps", "v", "av"]
        for rep in ["f"]
        if row["use_in_table"]
    ]


def gevp_E0_means(wildcards):
    return [
        f"intermediary_data/{dir_template}/meson_gevp_E0_{rep}_{channel}_mean.csv".format(**row)
        for row in metadata.to_dict(orient="records")
        for channel in ["ps", "v", "t", "av", "at", "s"]
        for rep in ["f"]
        if row["use_in_table"]
    ]


def mPCAC_means(wildcards):
    return [
        f"intermediary_data/{dir_template}/mpcac_mean.csv".format(**row)
        for row in metadata.to_dict(orient="records")
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
        if row["use_in_table"]
    ] 

def top_charge_means(wildcards):
    return [
        f"intermediary_data/{dir_template}/top_charge_mean.csv".format(**row)
        for row in metadata.to_dict(orient="records")
        if row["use_in_table"]
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

def continuum_massless_extrapolation_mass_ansatze(wildcards ):
    return [
        f"intermediary_data/extrapolation_results/f_{channel}_extp_{wildcards.ansatz}_mass_mean.csv".format()
        for channel in ["v", "t", "av", "at", "s"]
    ]


def continuum_massless_extrapolation_decay_ansatze(wildcards ):
    return [
        f"intermediary_data/extrapolation_results/f_{channel}_extp_{wildcards.ansatz}_decayconstant_mean.csv".format()
        for channel in ["ps", "v", "av"]
    ]


def continuum_massless_extrapolation_decayconstant(wildcards):
    return [
        f"intermediary_data/extrapolation_results/f_{channel}_extp_decayconstant_mean.csv".format()
        for channel in ["ps", "v", "av"]
    ]


def continuum_massless_extrapolation_mass_a2(wildcards):
    return [
        f"intermediary_data/extrapolation_results/f_{channel}_extp_a2_mass_mean.csv".format()
        for channel in ["v", "t", "av", "at", "s"]
    ]


def continuum_massless_extrapolation_decayconstant_a2(wildcards):
    return [
        f"intermediary_data/extrapolation_results/f_{channel}_extp_a2_decayconstant_mean.csv".format()
        for channel in ["ps", "v", "av"]
    ]


def continuum_massless_extrapolation_mass_a2_am2(wildcards):
    return [
        f"intermediary_data/extrapolation_results/f_{channel}_extp_a2_am2_mass_mean.csv".format()
        for channel in ["v", "t", "av", "at", "s"]
    ]


def continuum_massless_extrapolation_decayconstant_a2_am2(wildcards):
    return [
        f"intermediary_data/extrapolation_results/f_{channel}_extp_a2_am2_decayconstant_mean.csv".format()
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
        mPCAC=mPCAC_means,
        script="src/tables/print_spectrum_meson_1.py",
    output:
        table="assets/tables/spectrum_meson_1.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.decay} {input.mass} {input.E0} {input.w0} {input.mPCAC} --output_file {output.table}"


rule f_meson_mass_table2:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        decay=decay_constant_means,
        mass=extraction_means,
        E0=gevp_E0_means,
        script="src/tables/print_spectrum_meson_2.py",
    output:
        table="assets/tables/spectrum_meson_2.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.decay} {input.mass} {input.E0}  --output_file {output.table}"


rule ens_table:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        plaq=plaq_means,
        w0=w0_means,
        mass=gevp_E0_means,
        top_charge=top_charge_means,
        script="src/tables/print_ens.py",
    output:
        table="assets/tables/ensemble.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.plaq} {input.w0} {input.mass} {input.top_charge} --output_file {output.table}"


rule continuum_massless_mass_ansatze:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=continuum_massless_extrapolation_mass_ansatze,
        script="src/tables/continuum_massless_mass_ansatze.py",
    output:
        table="assets/tables/le_coefficients_mass_ansatz_{ansatz}.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.table} --ansatz {wildcards.ansatz}"

rule continuum_massless_decay_ansatze:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=continuum_massless_extrapolation_decay_ansatze,
        script="src/tables/continuum_massless_decayconstant_ansatze.py",
    output:
        table="assets/tables/le_coefficients_decayconstant_ansatz_{ansatz}.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.table} --ansatz {wildcards.ansatz}"



rule continuum_massless_mass_a2_am2:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=continuum_massless_extrapolation_mass_a2_am2,
        script="src/tables/continuum_massless_mass_a2_am2.py",
    output:
        table="assets/tables/clean_le_coefficients_mass_a2_am2.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.table} --ansatz a2_am2"


rule continuum_massless_decay_a2_am2:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=continuum_massless_extrapolation_decayconstant_a2_am2,
        script="src/tables/continuum_massless_decayconstant_a2_am2.py",
    output:
        table="assets/tables/clean_le_coefficients_decayconstant_a2_am2.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.table} --ansatz a2_am2"


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


rule continuum_massless_mass_a2:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=continuum_massless_extrapolation_mass_a2,
        script="src/tables/continuum_massless_mass_a2.py",
    output:
        table="assets/tables/nlo_coefficients_mass_a2.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.table}"


rule continuum_massless_decayconstant_a2:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=continuum_massless_extrapolation_decayconstant_a2,
        script="src/tables/continuum_massless_decayconstant_a2.py",
    output:
        table="assets/tables/nlo_coefficients_decayconstant_a2.tex",
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


rule g_VPP_KSRF:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        mv="intermediary_data/extrapolation_results/f_v_extp_a2_mass_samples.json",
        fps="intermediary_data/extrapolation_results/f_ps_extp_a2_decayconstant_samples.json",
        script="src/tables/g_VPP_KSRF.py",
    output:
        table="assets/tables/g_VPP_KSRF.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.mv} {input.fps} --output_file {output.table}"
