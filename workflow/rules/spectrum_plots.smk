from functools import partial

metadata = pd.read_csv("metadata/ensemble_metadata.csv")


def extraction_samples(wildcards):
    return [
        f"intermediary_data/{dir_template}/meson_extraction_{rep}_{channel}_samples.json".format(**row)
        for row in metadata.to_dict(orient="records")
        for channel in ["ps", "v"]
        for rep in ["f"]
    ]


def mass_gevp_samples(wildcards):
    return [
        f"intermediary_data/{dir_template}/meson_gevp_{rep}_{channel}_samples.json".format(**row)
        for row in metadata.to_dict(orient="records")
        for channel in ["ps", "v"]
        for rep in ["f"]
    ]


def decay_samples(wildcards):
    return [
        f"intermediary_data/{dir_template}/decay_constant_{rep}_{channel}_samples.json".format(**row)
        for row in metadata.to_dict(orient="records")
        for channel in ["ps", "v"]
        for rep in ["f"]
    ]


rule plot_extrapolations_meson_mass:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=mass_gevp_samples,
        script="src/plots/w0mps_vs_meson.py",
    output:
        plot_data="assets/plots/m2_all_con_sp4as.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --plot_styles {plot_styles} --plot_file_data {output.plot_data}"


rule plot_extrapolations_meson_decay:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=mass_gevp_samples,
        decay=decay_samples,
        script="src/plots/amps_vs_decay.py",
    output:
        plot_data="assets/plots/dec2_all_con_sp4as.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} {input.decay} --plot_styles {plot_styles} --plot_file {output.plot_data}"
