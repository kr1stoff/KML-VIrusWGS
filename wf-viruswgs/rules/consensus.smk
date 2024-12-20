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
        # * 不用这个参数, 但是要建索引
        csi=rules.bcftools_index.output,
    output:
        cons1="consensus/{sample}.rawcons.fa",
        cons2="consensus/{sample}.consensus.fa",
        noIUPAC="consensus/{sample}.no_IUPAC_masked.fa",
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
        awk '/>/ {{print; next}} {{gsub(/[rywskmhbvd]/,"N"); print}}' {input.ref} > {output.noIUPAC} 2> {log}
        bcftools consensus {params.extra} -p {wildcards.sample}_ -f {output.noIUPAC} -o {output.cons1} {input.vcf} 2>> {log}
        sed 's/_.*//' {output.cons1} > {output.cons2} 2>> {log}
        """
