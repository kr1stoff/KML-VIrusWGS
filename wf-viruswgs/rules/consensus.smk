rule bgzip:
    input:
        rules.freebayes.output.vcf,
    output:
        "consensus/{sample}.vcf.gz",
    params:
        extra="",  # optional
    threads: 1
    log:
        "logs/consensus/{sample}.bgzip.log",
    benchmark:
        "logs/consensus/{sample}.bgzip.bm"
    conda:
        config["conda"]["basic"]
    wrapper:
        f"file:{workflow.basedir}/wrappers/bio/bgzip"


rule bcftools_index:
    input:
        rules.bgzip.output,
    output:
        "consensus/{sample}.vcf.gz.csi",
    log:
        "logs/consensus/{sample}.index.log",
    benchmark:
        "logs/consensus/{sample}.index.bm"
    conda:
        config["conda"]["basic"]
    params:
        extra="",  # optional parameters for bcftools index
    wrapper:
        f"file:{workflow.basedir}/wrappers/bio/bcftools/index"


rule bcftools_consensus:
    input:
        vcf=rules.bgzip.output,
        ref=rules.bedtools_maskfasta.output,
    output:
        "consensus/{sample}.consensus.fa",
    log:
        "logs/consensus/{sample}.consensus.log",
    benchmark:
        "logs/consensus/{sample}.consensus.bm"
    conda:
        config["conda"]["basic"]
    params:
        extra="--haplotype A",  # optional parameters for bcftools consensus
    shell:
        """
        bcftools consensus {params.extra} -p {wildcards.sample}_ -f {input.ref} -o {output} {input.vcf}
        """
