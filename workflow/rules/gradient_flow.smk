W0_threshold = 0.28125

rule w0:
    input:
        data="raw_data/flows/{dir}/out_wflow",
        script="src/flow.py",
    output:
        mean="intermediary_data/{dir}/w0_mean.csv",
        samples="intermediary_data/{dir}/w0_samples.csv",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python {input.script} {input.data} {W0_threshold} --output_file_mean {output.mean} --output_file_samples {output.samples}"
