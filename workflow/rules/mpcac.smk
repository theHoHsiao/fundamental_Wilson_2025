import pandas as pd


metadata = pd.read_csv("metadata/ensemble_metadata.csv")
metadata_query = "Nc == {Nc} & Nt == {Nt} & Ns == {Ns} & beta == {beta} & nAS == {nAS} & mAS == {mAS}"

dir_template = "Sp{Nc}b{beta}nAS{nAS}mAS{mAS}T{Nt}L{Ns}"


rule fit_pcac:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        metadata=lookup(within=metadata, query=metadata_query),
    input:
        data="data_assets/correlators_wall.h5",
        script="src/mpcac.py",
    output:
        mean=f"intermediary_data/{dir_template}/mpcac_mean.csv",
        samples=f"intermediary_data/{dir_template}/mpcac_samples.json",
        plot=f"intermediary_data/{dir_template}/pcac_eff_mass.pdf",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file_mean {output.mean} --output_file_samples {output.samples} --effmass_plot_file {output.plot} --plot_styles {plot_styles} --ensemble_name {params.metadata.ensemble_name} --beta {params.metadata.beta} --mAS {params.metadata.mAS} --Nt {params.metadata.Nt} --Ns {params.metadata.Ns} --plateau_start {params.metadata.mpcac_plateau_start} --plateau_end {params.metadata.mpcac_plateau_end} --min_trajectory {params.metadata.init_conf} --max_trajectory {params.metadata.final_conf} --trajectory_step {params.metadata.delta_conf_spectrum}"


def all_pcac_data(wildcards):
    return [
        f"intermediary_data/{dir_template}/mpcac_mean.csv".format(**row)
        for row in metadata.to_dict(orient="records")
        if row["use_in_main_plots"]
    ]


def fittable_pcac_data(wildcards):
    return [
        f"intermediary_data/{dir_template}/mpcac_mean.csv".format(**row)
        for row in metadata.to_dict(orient="records")
        if row["use_in_PCAC_fits"] and row["use_in_main_plots"]
    ]


rule plot_pcac:
    input:
        plot_data=all_pcac_data,
        fit_data=fittable_pcac_data,
        script="src/plots/pcac_fits.py",
    output:
        plot="assets/plots/mpcac_vs_m0.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m src.plots.pcac_fits {input.plot_data} --fit_filenames {input.fit_data} --plot_filename {output.plot} --plot_styles {plot_styles}"
