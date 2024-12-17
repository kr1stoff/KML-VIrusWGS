ruleorder: fastp_pe > fastp_se


rule fastp_se:
    input:
        sample=[".rawdata/{sample}.1.fastq.gz"],
    output:
        trimmed="qc/fastp/{sample}.1.fastq",
        failed="qc/fastp/{sample}.1.failed.fastq",
        html="qc/fastp/{sample}.html",
        json="qc/fastp/{sample}.json",
    log:
        "logs/fastp/{sample}.log",
    benchmark:
        "logs/fastqc/{sample}.bm"
    params:
        extra="-q 15 -u 40 -l 15 --cut_right --cut_window_size 4 --cut_mean_quality 20 --correction",
    threads: config["threads"]["low"]
    wrapper:
        f"file:{workflow.basedir}/wrappers/bio/fastp"


use rule fastp_se as fastp_pe with:
    input:
        sample=[".rawdata/{sample}.1.fastq.gz", ".rawdata/{sample}.2.fastq.gz"],
    output:
        trimmed=["qc/fastp/{sample}.1.fastq", "qc/fastp/{sample}.2.fastq"],
        unpaired1="qc/fastp/{sample}.u1.fastq",
        unpaired2="qc/fastp/{sample}.u2.fastq",
        merged="qc/fastp/{sample}.merged.fastq",
        failed="qc/fastp/{sample}.failed.fastq",
        html="qc/fastp/{sample}.html",
        json="qc/fastp/{sample}.json",
    params:
        extra="--merge -q 15 -u 40 -l 15 --cut_right --cut_window_size 4 --cut_mean_quality 20 --correction",


rule qc_stat:
    input:
        expand("qc/fastp/{sample}.json", sample=config["samples"]),
    output:
        "qc/fastp/fastp.stats.tsv",
    benchmark:
        "logs/fastp/qc_stat.bm"
    log:
        "logs/fastp/qc_stat.log",
    script:
        "../scripts/fastp_all_samples_qc.py"
