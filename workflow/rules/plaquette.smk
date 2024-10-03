from glob import glob

rule plot_smallvolume_plaquettes:
    input:
        data=glob("raw_data/hmc/out_hmc_8x8x8x8_*"),
        script="src/plaquette_hmc.py",
    output:
        plot="assets/plots/plaquette_phasediagram.pdf",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python {input.script} {input.data} --plot_filename {output.plot} --plot_styles {plot_styles}"
