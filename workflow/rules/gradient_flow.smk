import pandas as pd

flow_metadata = pd.read_csv("metadata/flow_meta.csv")
metadata_query = "Nc == {Nc} & Nt == {Nt} & Ns == {Ns} & beta == {beta} & nAS == {nAS} & mAS == {mAS}"

W0_threshold = 0.28125

dir_template = "Sp{Nc}b{beta}nAS{nAS}mAS{mAS}T{Nt}L{Ns}"


rule w0:
    params:
        metadata=lookup(within=flow_metadata, query=metadata_query),
    input:
        data=f"raw_data/flows/{dir_template}/out_wflow",
        script="src/flow.py",
    output:
        mean=f"intermediary_data/{dir_template}/w0_mean.csv",
        samples=f"intermediary_data/{dir_template}/w0_samples.csv",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python {input.script} {input.data} {W0_threshold} --output_file_mean {output.mean} --output_file_samples {output.samples} --min_trajectory {params.metadata.init_conf} --max_trajectory {params.metadata.final_conf} --trajectory_step {params.metadata.delta_conf}"
