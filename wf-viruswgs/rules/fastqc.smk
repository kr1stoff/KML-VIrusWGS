rule fastqc:
    input:
        ".rawdata/{sample}.1.fastq.gz",
    output:
        html="qc/fastqc/{sample}.1.html",
        zip="qc/fastqc/{sample}_1_fastqc.zip",
    params:
        extra="--quiet",
    log:
        "logs/fastqc/{sample}.1.log",
    benchmark:
        "logs/fastqc/{sample}.1.bm"
    threads: config["threads"]["low"]
    conda:
        config["conda"]["basic"]
    resources:
        mem_mb=1024,
    wrapper:
        f"file:{workflow.basedir}/wrappers/bio/fastqc"


use rule fastqc as fastqc2 with:
    input:
        ".rawdata/{sample}.2.fastq.gz",
    output:
        html="qc/fastqc/{sample}.2.html",
        zip="qc/fastqc/{sample}_2_fastqc.zip",
    log:
        "logs/fastqc/{sample}.2.log",
    benchmark:
        "logs/fastqc/{sample}.2.bm"
