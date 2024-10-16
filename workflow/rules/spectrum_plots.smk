from functools import partial


metadata = pd.read_csv("metadata/ensemble_metadata.csv")


def all_samples(wildcards, observables):
    return [
        f"intermediary_data/{dir_template}/{observable}_samples.json".format(**row)
        for observable in observables
        for row in metadata.to_dict(orient="records")
        if row["use_in_main_plots"]
        if row["use_smear"]
    ]


rule plot_w0mps_vs_w0mv:
    input:
        data=partial(all_samples, observables=["w0", "smear_meson_ps", "smear_meson_v", "smear_meson_t","smear_meson_av","smear_meson_at","smear_meson_s", "GEVP_meson_rhoE1"]),
        script="src/plots/w0mps_vs_meson.py",
    output:
        plot="assets/plots/m2_all_con_sp4as.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m src.plots.w0mps_vs_meson {input.data} --plot_styles {plot_styles} --plot_file {output.plot}"
