from glob import glob


rule avg_plaquette:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        metadata=lookup(within=metadata, query=metadata_query),
    input:
        data="data_assets/corr_sp4_FUN.h5",
        script="src/plaquette_skip.py",
    wildcard_constraints:
        Ns=r"\d+",
    output:
        mean=f"intermediary_data/{dir_template}/plaquette_mean.csv",
        samples=f"intermediary_data/{dir_template}/plaquette_samples.json",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file_mean {output.mean} --output_file_samples {output.samples} --ensemble_name {params.metadata.ensemble_name} --beta {params.metadata.beta} --mF {params.metadata.mF} --Nt {params.metadata.Nt} --Ns {params.metadata.Ns}"
        " --min_trajectory {params.metadata.init_conf} --max_trajectory {params.metadata.final_conf}"
        " --trajectory_step {params.metadata.delta_conf_obs} --trajectory_step_auto {params.metadata.delta_conf}"


rule avg_plaquette_hmc:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        metadata=lookup(
            within=metadata, query=metadata_query + " & start_type == '{start}'"
        ),
    input:
        data="data_assets/hmc.h5",
        script="src/plaquette_hmc.py",
    output:
        mean=f"intermediary_data/{dir_template}_{{start}}start/plaquette_mean.csv",
    wildcard_constraints:
        Ns=r"\d+",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --ensemble_name {params.metadata.ensemble_name} --beta {params.metadata.beta} --mAS {params.metadata.mAS} --Nt {params.metadata.Nt} --Ns {params.metadata.Ns} --start {wildcards.start} --output_file {output.mean}"


def main_plaquette_data(wildcards):
    return [
        f"intermediary_data/{dir_template}/{basename}.csv".format(**row)
        for basename in ["plaquette_mean", "meson_gevp_E0_f_ps_mean"]
        for row in metadata.to_dict(orient="records")
        if row["use_in_main_plots"]
    ]


rule tabulate_largevolume_plaquettes:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=main_plaquette_data,
        script="src/tables/plaquette.py",
    output:
        table="assets/tables/plaquette_table.tex",
        definitions="assets/definitions/heavy_ps_limit.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.table} --definitions_file={output.definitions}"
