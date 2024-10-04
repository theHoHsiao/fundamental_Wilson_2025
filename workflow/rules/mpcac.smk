import pandas as pd


pcac_metadata = pd.read_csv("metadata/pcac_meta.csv")
metadata_query = "Nc == {Nc} & Nt == {Nt} & Ns == {Ns} & beta == {beta} & nAS == {nAS} & mAS == {mAS}"

dir_template = "Sp{Nc}b{beta}nAS{nAS}mAS{mAS}T{Nt}L{Ns}"


rule fit_pcac:
    params:
        metadata=lookup(within=pcac_metadata, query=metadata_query),
    input:
        data="data_assets/correlators_wall.h5",
        script="src/mpcac.py",
    output:
        mean=f"intermediary_data/{dir_template}/mpcac.csv",
        plot=f"intermediary_data/{dir_template}/pcac_eff_mass.pdf",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python {input.script} {input.data} --output_file {output.mean} --effmass_plot_file {output.plot} --plot_styles {plot_styles} --beta {params.metadata.beta} --mAS {params.metadata.mAS} --Nt {params.metadata.Nt} --Ns {params.metadata.Ns} --plateau_start {params.metadata.mpcac_plateau_start} --plateau_end {params.metadata.mpcac_plateau_end} --min_trajectory {params.metadata.init_conf} --max_trajectory {params.metadata.final_conf} --trajectory_step {params.metadata.delta_conf}"


def all_pcac_data(wildcards):
    return [
        f"intermediary_data/{dir_template}/mpcac.csv".format(**row)
        for row in pcac_metadata.to_dict(orient="records")
    ]


rule plot_pcac:
    input:
        data=all_pcac_data,
        script="src/plots/pcac_fits.py",
    output:
        plot="assets/plots/mpcac_vs_m0.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m src.plots.pcac_fits {input.data} --plot_filename {output.plot} --plot_styles {plot_styles}"
