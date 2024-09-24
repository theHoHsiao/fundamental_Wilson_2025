import pandas as pd

flow_metadata = pd.read_csv("metadata/flow_meta.csv")
metadata_query = "Nc == {Nc} & Nt == {Nt} & Ns == {Ns} & beta == {beta} & nAS == {nAS} & mAS == {mAS}"

W0_threshold = 0.28125

dir_template = "Sp{Nc}b{beta}nAS{nAS}mAS{mAS}T{Nt}L{Ns}"


rule w0:
    params:
        metadata=lookup(within=flow_metadata, query=metadata_query),
    input:
        data=f"raw_data/flows/{dir_template}/out_wflow",
        script="src/flow.py",
    output:
        mean=f"intermediary_data/{dir_template}/w0_mean.csv",
        samples=f"intermediary_data/{dir_template}/w0_samples.csv",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python {input.script} {input.data} {W0_threshold} --output_file_mean {output.mean} --output_file_samples {output.samples} --min_trajectory {params.metadata.init_conf} --max_trajectory {params.metadata.final_conf} --trajectory_step {params.metadata.delta_conf} --ensemble_name {params.metadata.ensemble_name}"


rule topological_charge:
    params:
        metadata=lookup(within=flow_metadata, query=metadata_query),
    input:
        data=f"raw_data/flows/{dir_template}/out_wflow",
        script="src/top_charge.py",
    output:
        data=f"intermediary_data/{dir_template}/top_charge.csv",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python {input.script} {input.data} --output_file {output.data} --min_trajectory {params.metadata.init_conf} --max_trajectory {params.metadata.final_conf} --ensemble_name {params.metadata.ensemble_name}"


def all_flow_data(wildcards):
    return [
        f"intermediary_data/{dir_template}/{basename}.csv".format(**row)
        for basename in ["w0_mean", "top_charge"]
        for row in flow_metadata.to_dict(orient="records")
    ]


rule gflow_table:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=all_flow_data,
        script="src/tables/Q_table.py",
    output:
        table="assets/tables/gflow_table.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.table}"
