from glob import glob


rule tabulate_largevolume_plaquettes:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data="data_assets/correlators_wall.h5",
        script="src/plaquette.py",
        metadata="metadata/flow_meta.csv",
    output:
        table="assets/tables/plaquette_table.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --metadata {input.metadata} --output_table {output.table}"


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
