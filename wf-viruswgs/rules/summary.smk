rule summary_table:
    input:
        depth_stat=rules.depth_stat.output,
        vtype=rules.virus_type.output,
    output:
        "summary/{sample}.summary.tsv",
    log:
        "logs/summary/{sample}.summary_table.log",
    benchmark:
        "logs/summary/{sample}.summary_table.bm"
    params:
        name="{sample}",
    script:
        "../scripts/summary_table.py"


rule concat_summary:
    input:
        expand("summary/{sample}.summary.tsv", sample=config["samples"]),
    output:
        "summary/all.summary.tsv",
    log:
        "logs/summary/all_samples.summary_table.log",
    benchmark:
        "logs/summary/all_samples.summary_table.bm"
    script:
        "../scripts/concat_summary.py"
