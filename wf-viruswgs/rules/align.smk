use rule bwa_mem_datasets_se as bwa_mem_ref_se with:
    input:
        reads=["qc/fastp/{sample}.1.fastq"],
        # idx=multiext(config["database"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
        idx=rules.ref_index.output,
    output:
        "align/{sample}.bam",
    log:
        "logs/align/{sample}_bwa_mem_ref.log",
    benchmark:
        "logs/align/{sample}_bwa_mem_ref.bm"


use rule bwa_mem_ref_se as bwa_mem_ref_pe with:
    input:
        ["qc/fastp/{sample}.1.fastq", "qc/fastp/{sample}.2.fastq"],


use rule samtools_coverage_datasets as samtools_coverage_ref with:
    input:
        "align/{sample}.bam",
    output:
        "align/{sample}.coverage",
    log:
        "logs/align/{sample}_samtools_coverage_ref.log",
    benchmark:
        "logs/align/{sample}_samtools_coverage_ref.bm"


rule samtools_index:
    input:
        "align/{sample}.bam",
    output:
        "align/{sample}.bam.bai",
    log:
        "logs/align/{sample}.samtools_index.log",
    benchmark:
        "logs/align/{sample}.samtools_index.bm"
    params:
        extra="",  # optional params string
    threads: config["threads"]["low"]  # This value - 1 will be sent to -@
    conda:
        config["conda"]["basic"]
    wrapper:
        f"file:{workflow.basedir}/wrappers/bio/samtools/index"
