from functools import partial

metadata = pd.read_csv("metadata/ensemble_metadata.csv")


def extraction_samples(wildcards):
    return [
        f"intermediary_data/{dir_template}/meson_extraction_{rep}_{channel}_samples.json".format(**row)
        for row in metadata.to_dict(orient="records")
        for channel in ["ps", "v"]
        for rep in ["f"]
        if row["use_in_extrapolation"]
    ]


def mass_gevp_samples(wildcards):
    return [
        f"intermediary_data/{dir_template}/meson_gevp_{rep}_{channel}_samples.json".format(**row)
        for row in metadata.to_dict(orient="records")
        for channel in ["ps", "v"]
        for rep in ["f"]
        if row["use_in_extrapolation"]
    ]



def w0_samples(wildcards):
    return [
        f"intermediary_data/{dir_template}/w0_samples.json".format(**row)
        for row in metadata.to_dict(orient="records")
        if row["use_in_extrapolation"]
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
        w0=w0_samples,
        script="src/plots/w0mps_vs_meson.py",
    output:
        plot_data="assets/plots/m2_all_con_sp4fund.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} {input.w0} --plot_styles {plot_styles} --plot_file_data {output.plot_data}"


rule plot_extrapolations_meson_decay:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        mass=mass_gevp_samples,
        decay=decay_samples,
        w0=w0_samples,
        script="src/plots/w0mps_vs_decay.py",
    output:
        plot_data="assets/plots/dec2_all_con_sp4fund.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.mass} {input.decay} {input.w0} --plot_styles {plot_styles} --plot_file {output.plot_data}"


rule plot_extrapolations_Rvps:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=mass_gevp_samples,
        w0=w0_samples,
        script="src/plots/w0mps_vs_Rvps.py",
    output:
        plot_data="assets/plots/Rvps_all_con_sp4fund.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} {input.w0} --plot_styles {plot_styles} --plot_file_data {output.plot_data}"

