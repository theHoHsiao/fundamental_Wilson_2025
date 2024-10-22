from glob import glob


rule avg_plaquette:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        metadata=lookup(within=metadata, query=metadata_query),
    input:
        data="data_assets/correlators_wall.h5",
        script="src/plaquette.py",
    output:
        mean=f"intermediary_data/{dir_template}/plaquette_mean.csv",
        samples=f"intermediary_data/{dir_template}/plaquette_samples.json",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file_mean {output.mean} --output_file_samples {output.samples} --ensemble_name {params.metadata.ensemble_name} --beta {params.metadata.beta} --mAS {params.metadata.mAS} --Nt {params.metadata.Nt} --Ns {params.metadata.Ns} --min_trajectory {params.metadata.init_conf} --max_trajectory {params.metadata.final_conf} --trajectory_step {params.metadata.delta_conf_spectrum}"


def main_plaquette_data(wildcards):
    return [
        f"intermediary_data/{dir_template}/plaquette_mean.csv".format(**row)
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
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file {output.table}"


rule plot_smallvolume_plaquettes:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=glob("raw_data/hmc/out_hmc_8x8x8x8_*"),
        script="src/plaquette_hmc.py",
    output:
        plot="assets/plots/plaquette_phasediagram.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --plot_filename {output.plot} --plot_styles {plot_styles}"
