from functools import partial

metadata = pd.read_csv("metadata/ensemble_metadata.csv")


def extraction_samples(wildcards):
    return [
        f"intermediary_data/{dir_template}/meson_extraction_{rep}_{channel}_samples.json".format(**row)
        for row in metadata.to_dict(orient="records")
        for channel in ["ps", "v", "av"]
        for rep in ["f"]
        if row["use_in_extrapolation"]
    ]


def mass_gevp_samples(wildcards):
    return [
        f"intermediary_data/{dir_template}/meson_gevp_{rep}_{channel}_samples.json".format(**row)
        for row in metadata.to_dict(orient="records")
        for channel in ["ps", "v", "t", "av", "at", "s"]
        for rep in ["f"]
        if row["use_in_extrapolation"]
    ]


def mass_smear_samples(wildcards):
    return [
        f"intermediary_data/{dir_template}/meson_smear_mass_{rep}_{channel}_samples.json".format(**row)
        for row in metadata.to_dict(orient="records")
        for channel in ["ps", "v", "t"]
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
        for channel in ["ps", "v", "av"]
        for rep in ["f"]
        if row["use_in_extrapolation"]

    ]


def mass_extp(wildcards, observables):
    return [
        f"intermediary_data/extrapolation_results/{observable}_samples.json".format(
            **row
        )
        for observable in observables
        for row in metadata.to_dict(orient="records")
        #if row["use_in_main_plots"]
        if row["use_in_extrapolation"]
    ]

rule plot_extrapolations_meson_mass:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=mass_gevp_samples,
        w0=w0_samples,
        fit_results=partial(
            mass_extp,
            observables=[
                "f_v_extp_mass",
                "f_t_extp_mass",
                "f_av_extp_mass",
                "f_at_extp_mass",
                "f_s_extp_mass",
            ],
        ),
        script="src/plots/w0mps_vs_meson.py",
    output:
        plot_data="assets/plots/m2_all_con_sp4fund.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} {input.w0} --plot_styles {plot_styles} --plot_file_data {output.plot_data} --fit_parameters {input.fit_results}"


rule plot_mpsL:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=mass_gevp_samples,
        w0=w0_samples,
        script="src/plots/mps_vs_mpsL.py",
    output:
        plot_data="assets/plots/mpsL.{plot_filetype}",
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
        fit_results=partial(
            mass_extp,
            observables=[
                "f_ps_extp_decayconstant",
                "f_v_extp_decayconstant",
                "f_av_extp_decayconstant",
            ],
        ),
        script="src/plots/w0mps_vs_decay.py",
    output:
        plot_data="assets/plots/dec2_all_con_sp4fund.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.mass} {input.decay} {input.w0} --plot_styles {plot_styles} --plot_file {output.plot_data} --fit_results {input.fit_results}"


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


rule plot_test:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        #data=mass_gevp_samples,
        w0=w0_samples,
        script="src/plots/test.py",
    output:
        plot_data="assets/plots/test.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.w0} --plot_styles {plot_styles} --plot_file_data {output.plot_data}"


rule plot_extrapolations_meson_mass_smear:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=mass_smear_samples,
        w0=w0_samples,
        script="src/plots/w0mps_vs_meson_smear.py",
    output:
        plot_data="assets/plots/m2_all_con_sp4fund_smear.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} {input.w0} --plot_styles {plot_styles} --plot_file_data {output.plot_data}"


rule compare_mass_smear_gevp:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=mass_smear_samples,
        data2=mass_gevp_samples,
        w0=w0_samples,
        script="src/plots/compare_smear_gevp.py",
    output:
        plot_data="assets/plots/compare_mass_smear_gevp.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} {input.data2} {input.w0} --plot_styles {plot_styles} --plot_file_data {output.plot_data}"


rule check_lat_a_extrapolations_meson_mass:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=mass_gevp_samples,
        w0=w0_samples,
        fit_results=partial(
            mass_extp,
            observables=[
                "f_v_extp_mass",
                "f_t_extp_mass",
                "f_av_extp_mass",
                "f_at_extp_mass",
                "f_s_extp_mass",
            ],
        ),
        script="src/plots/lat_a_vs_meson.py",
    output:
        plot_data="assets/plots/m2_all_lat_a_sp4fund.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} {input.w0} --plot_styles {plot_styles} --plot_file_data {output.plot_data} --fit_parameters {input.fit_results}"


rule check_lat_a_extrapolations_meson_decay:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        mass=mass_gevp_samples,
        decay=decay_samples,
        w0=w0_samples,
        fit_results=partial(
            mass_extp,
            observables=[
                "f_ps_extp_decayconstant",
                "f_v_extp_decayconstant",
                "f_av_extp_decayconstant",
            ],
        ),
        script="src/plots/lat_a_vs_decay.py",
    output:
        plot_data="assets/plots/dec2_all_lat_a_sp4fund.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.mass} {input.decay} {input.w0} --plot_styles {plot_styles} --plot_file {output.plot_data} --fit_results {input.fit_results}"
