import pandas as pd
from functools import partial


all_metadata = pd.read_csv("metadata/ensemble_metadata.csv")

metadata_query = "Nc == {Nc} & Nt == {Nt} & Ns == {Ns} & beta == {beta} & nF == {nF} & mF == {mF}"
metadata_lookup = partial(lookup, within=all_metadata, query=metadata_query)

dir_template = "Sp{Nc}b{beta}nF{nF}mF{mF}T{Nt}L{Ns}"

def extraction_means(wildcards):
    return [
        f"intermediary_data/{dir_template}/meson_extraction_{rep}_{channel}_mean.csv".format(**row)
        for row in metadata.to_dict(orient="records")
        for channel in ["ps", "v"]
        for rep in ["f"]
    ] 



def decay_constant_means(wildcards):
    return [
        f"intermediary_data/{dir_template}/decay_constant_{rep}_{channel}_mean.csv".format(**row)
        for row in metadata.to_dict(orient="records")
        for channel in ["ps", "v"]
        for rep in ["f"]
    ] 


rule f_meson_mass_table:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        decay=decay_constant_means,
        mass=extraction_means, #passing more data than the code needs to generate all the results. To be changed...
        script="src/tables/print_matrix_element_meson_f.py",
    output:
        table="assets/tables/matrix_element_meson_f.tex",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.decay} {input.mass} --output_file {output.table}"
