from glob import glob

parsing_base = "src/HiRep_parsing"


rule package_smeared:
    params:
        file_dir="raw_data/corr/{smearing}",
        script_file_name="scripts/write_meson_{smearing}.jl",
    input:
        files=glob("raw_data/corr/{smearing}/*"),
        script=f"{parsing_base}/scripts/write_meson_{{smearing}}.jl",
    output:
        h5=protected("data_assets/correlators_{smearing}.h5"),
    conda:
        "../envs/hirep_parsing.yml"
    # Start packaging early,
    # since it is time consuming and many other processes depend on it
    priority: 10
    shell:
        "cd {parsing_base} && julia {params.script_file_name} ../../{params.file_dir} ../../{output.h5}"


rule package_gflow:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        files=glob("raw_data/*/topology/out/out_flow"),
        script="src/package_flows.py",
    output:
        h5=protected("data_assets/nf2_gflow.h5"),
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.files} --h5_filename {output.h5}"


rule package_hmc:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        files=glob("raw_data/hmc/out_hmc_*"),
        script="src/package_hmc.py",
    output:
        h5=protected("data_assets/hmc.h5"),
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.files} --h5_filename {output.h5}"
