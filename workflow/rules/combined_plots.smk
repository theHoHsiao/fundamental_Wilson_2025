from functools import partial


metadata = pd.read_csv("metadata/ensemble_metadata.csv")


def all_samples(wildcards, observables):
    return [
        f"intermediary_data/{dir_template}/{observable}_samples.json".format(**row)
        for observable in observables
        for row in metadata.to_dict(orient="records")
        if row["use_in_main_plots"]
    ]


rule plot_w0_vs_mpcac:
    input:
        data=partial(all_samples, observables=["w0", "mpcac"]),
        script="src/plots/w0_vs_mpcac.py",
    output:
        plot="assets/plots/w0_vs_mpcac.{plot_filetype}",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m src.plots.w0_vs_mpcac {input.data} --plot_styles {plot_styles} --plot_file {output.plot}"
