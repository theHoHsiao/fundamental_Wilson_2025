import pandas as pd

metadata = pd.read_csv("metadata/ensemble_metadata.csv")

W0_threshold = 0.28125


rule w0:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        metadata=lookup(within=metadata, query=metadata_query),
    input:
        data=f"data_assets/flows.h5",
        script="src/flow.py",
    output:
        mean=f"intermediary_data/{dir_template}/w0_mean.csv",
        samples=f"intermediary_data/{dir_template}/w0_samples.json",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} {W0_threshold} --output_file_mean {output.mean} --output_file_samples {output.samples} --min_trajectory {params.metadata.init_conf} --max_trajectory {params.metadata.final_conf} --trajectory_step {params.metadata.delta_conf_w0} --ensemble_name {params.metadata.ensemble_name} --beta {params.metadata.beta} --mAS {params.metadata.mAS} --Nt {params.metadata.Nt} --Ns {params.metadata.Ns}"


rule w0_flow_plot:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        metadata=lookup(within=metadata, query=metadata_query),
    input:
        data=f"data_assets/flows.h5",
        script="src/plot_w_flow.py",
    output:
        plot=f"assets/plots/w0_flow_{dir_template}.{{plot_filetype}}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} {W0_threshold} --plot_file {output.plot} --plot_styles {plot_styles} --min_trajectory {params.metadata.init_conf} --max_trajectory {params.metadata.final_conf} --trajectory_step {params.metadata.delta_conf_w0} --beta {params.metadata.beta} --mAS {params.metadata.mAS} --Nt {params.metadata.Nt} --Ns {params.metadata.Ns}"


rule topological_charge:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        metadata=lookup(within=metadata, query=metadata_query),
    input:
        data=f"data_assets/flows.h5",
        script="src/top_charge.py",
    output:
        data=f"intermediary_data/{dir_template}/top_charge_mean.csv",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.data} --min_trajectory {params.metadata.init_conf} --max_trajectory {params.metadata.final_conf} --ensemble_name {params.metadata.ensemble_name} --beta {params.metadata.beta} --mAS {params.metadata.mAS} --Nt {params.metadata.Nt} --Ns {params.metadata.Ns}"


rule topological_charge_history_plot:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        metadata=lookup(within=metadata, query=metadata_query),
    input:
        data=f"data_assets/flows.h5",
        script="src/top_charge.py",
    output:
        plot=f"assets/plots/top_charge_history_{dir_template}.pdf",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file /dev/null --min_trajectory {params.metadata.init_conf} --max_trajectory {params.metadata.final_conf} --plot_file {output.plot} --plot_styles {plot_styles} --beta {params.metadata.beta} --mAS {params.metadata.mAS} --Nt {params.metadata.Nt} --Ns {params.metadata.Ns}"


def all_flow_data(wildcards):
    return [
        f"intermediary_data/{dir_template}/{basename}.csv".format(**row)
        for basename in ["w0_mean", "top_charge"]
        for row in metadata.to_dict(orient="records")
        if row["use_in_main_plots"]
    ]


rule gflow_table:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=all_flow_data,
        script="src/tables/Q_table.py",
    output:
        table="assets/tables/gflow_table.tex",
        definitions="assets/definitions/gflow_incomplete_ensembles.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.table} --definitions_file {output.definitions}"
