from functools import partial


metadata = pd.read_csv("metadata/ensemble_metadata.csv")


def all_samples(wildcards, observables):
    return [
        f"intermediary_data/{dir_template}/{observable}_samples.json".format(**row)
        for observable in observables
        for row in metadata.to_dict(orient="records")
        if row["use_in_main_plots"]
    ]

def extp_samples(wildcards, observables):
    return [
        f"intermediary_data/{dir_template}/{observable}_samples.json".format(**row)
        for observable in observables
        for row in metadata.to_dict(orient="records")
        if row["use_in_extrapolation"]
    ]

def mass_extp(wildcards, observables):
    return [
        f"intermediary_data/extrapolation_results/{observable}_samples.json".format(**row)
        for observable in observables
        for row in metadata.to_dict(orient="records")
        if row["use_in_main_plots"]
        if row["use_in_extrapolation"]
    ]

def volume_samples(wildcards, observables):
    return [
        f"intermediary_data/{dir_template}/{observable}_samples.json".format(**row)
        for observable in observables
        for row in metadata.to_dict(orient="records")
        if row["use_in_finite_volume"]
    ]

def ASB2s_samples(wildcards, observables):
    return [
        f"intermediary_data/{dir_template}/{observable}_samples.json".format(**row)
        for observable in observables
        for row in metadata.to_dict(orient="records")
        if row["ensembles_B6p7"]
    ]

rule plot_finite_volume:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=partial(volume_samples, observables=["meson_ps","meson_v"]),
        script="src/plots/finite_volume.py"
    output:
        plot="assets/plots/mps_mv_vs_mpsL.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --plot_styles {plot_styles} --plot_file {output.plot}"

rule plot_mpsmv_vs_mpcac:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=partial(ASB2s_samples, observables=["meson_ps","meson_v", "mpcac"]),
    output:
        plot="assets/plots/mpsmv_vs_mpcac_b6p7.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --plot_styles {plot_styles} --plot_file {output.plot}"

rule plot_mpsfps_vs_mpcac:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=partial(ASB2s_samples, observables=["meson_ps", "plaquette","mpcac"]),
        script="src/plots/mpsfps_vs_mpcac.py"
    output:
        plot="assets/plots/mpsfps_vs_mpcac_b6p7.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --plot_styles {plot_styles} --plot_file {output.plot}"

rule plot_meson_mass_vs_fermion_mass:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=partial(ASB2s_samples, observables=["meson_ps","meson_v", "meson_av","meson_at", "meson_s", "mpcac"]),
        script="src/plots/meson_vs_fermion.py"
    output:
        plot="assets/plots/meson_masses_b6p7_m0.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --plot_styles {plot_styles} --plot_file {output.plot}"

rule plot_decay_constant_vs_fermion_mass:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=partial(ASB2s_samples, observables=["meson_ps","meson_v", "meson_av", "mpcac", "plaquette"]),
        script="src/plots/decay_vs_fermion.py"
    output:
        plot="assets/plots/meson_decay_b6p7_m0.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --plot_styles {plot_styles} --plot_file {output.plot}"

rule plot_meson_mass_fps_vs_mpcac:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=partial(ASB2s_samples, observables=["meson_ps","meson_v", "meson_av","meson_at", "meson_s", "mpcac", "plaquette"]),
        script="src/plots/meson_fps_vs_mpcac.py"
    output:
        plot="assets/plots/mmfps_vs_mpcac_b6p7.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --plot_styles {plot_styles} --plot_file {output.plot}"

rule plot_mv_vs_mps:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=partial(ASB2s_samples, observables=["meson_ps","meson_v", "w0"]),
        fund_data="data_assets/meson_meta_fund.csv",
        script="src/plots/mv_vs_mps.py"
    output:
        plot="assets/plots/m2v_vs_m2ps_GF_b6p7.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --plot_styles {plot_styles} --plot_file {output.plot}"

rule plot_extrapolations_meson_mass:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=partial(extp_samples, observables=["w0", "smear_meson_ps", "smear_meson_v", "smear_meson_t","smear_meson_av","smear_meson_at","smear_meson_s", "gevp_smear_meson_rhoE1"]),
        fit_results=partial(mass_extp, observables=["v_extp_mass","t_extp_mass","av_extp_mass","at_extp_mass","s_extp_mass","rhoE1_extp_mass"]),
        script="src/plots/w0mps_vs_meson.py",
    output:
        plot="assets/plots/m2_all_con_sp4as.{plot_filetype}",
        plot2="assets/plots/meson_spectrum_con.{plot_filetype}"
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --plot_styles {plot_styles} --plot_file {output.plot} --plot_file2 {output.plot2} --fit_parameters {input.fit_results}"

rule plot_extrapolations_meson_decay:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=partial(extp_samples, observables=["w0", "meson_ps", "meson_v","meson_av", "plaquette"]),
        fit_results=partial(mass_extp, observables=["ps_extp_decay", "v_extp_decay", "av_extp_decay"]),
        script="src/plots/w0mps_vs_decay.py",
    output:
        plot="assets/plots/f2_up_con_sp4as.{plot_filetype}",
        plot2="assets/plots/f2_low_con_sp4as.{plot_filetype}"
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --plot_styles {plot_styles} --plot_file {output.plot} --plot_file2 {output.plot2} --fit_parameters {input.fit_results}"


rule plot_R_mvfps_vs_mps:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=partial(extp_samples, observables=["meson_ps","meson_v", "w0", "plaquette"]),
        fund_data="data_assets/mv_fps_fund.csv",
        extp_data="intermediary_data/extrapolation_results/R_mvdfps_extp_samples.json",
        script="src/plots/R_mvfps_vs_mps.py"
    output:
        plot="assets/plots/mvfps_vs_m2ps_GF_F_vs_AS.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --plot_styles {plot_styles} --plot_file {output.plot}"

rule plot_R_mvprime_dmv_vs_mps:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=partial(extp_samples, observables=["smear_meson_ps","smear_meson_v", "w0", "gevp_smear_meson_rhoE1"]),
        script="src/plots/R_mvprime_dmv_vs_mps.py"
    output:
        plot="assets/plots/excited_vector_ratio.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --plot_styles {plot_styles} --plot_file {output.plot}"

rule plot_w0_vs_mps:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=partial(all_samples, observables=["meson_ps","w0"]),
        script="src/plots/w0_vs_mps.py"
    output:
        plot="assets/plots/w0vsmps2.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --plot_styles {plot_styles} --plot_file {output.plot}"

rule plot_R_mps_dfps_vs_mps:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=partial(all_samples, observables=["meson_ps","w0", "plaquette"]),
        script="src/plots/R_mps_dfps_vs_mps.py"
    output:
        plot="assets/plots/mfpsvsmps2.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --plot_styles {plot_styles} --plot_file {output.plot}"

rule plot_gmor:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        data=partial(ASB2s_samples, observables=["meson_ps","w0", "plaquette", "mpcac"]),
        script="src/plots/gmor.py"
    output:
        plot="assets/plots/gmor_b6p7.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --plot_styles {plot_styles} --plot_file {output.plot}"
