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


rule samtools_depth:
    input:
        "align/{sample}.bam",