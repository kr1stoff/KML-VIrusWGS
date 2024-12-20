rule bedtools_maskfasta:
    input:
        low_depth=rules.low_depth_region.output,
        raw_ref=rules.find_ref.output.ref,
    output:
        "variants/{sample}.masked.fa",
    log:
        "logs/variants/{sample}.masked.log",
    benchmark:
        "logs/variants/{sample}.masked.bm"
    conda:
        config["conda"]["basic"]
    shell:
        "bedtools maskfasta -fi {input.raw_ref} -bed {input.low_depth} -fo {output} 2> {log}"


rule samtools_faidx:
    input:
        rules.bedtools_maskfasta.output,
    output:
        "variants/{sample}.masked.fa.fai",
    log:
        "logs/variants/{sample}.samtools_faidx.log",
    benchmark:
        "logs/variants/{sample}.samtools_faidx.bm"
    params:
        extra="",
    conda:
        config["conda"]["basic"]
    wrapper:
        f"file:{workflow.basedir}/wrappers/bio/samtools/faidx"


rule freebayes:
    input:
        alns="align/{sample}.bam",
        idxs=rules.samtools_index.output,
        ref=rules.bedtools_maskfasta.output,
        # * 不用这个参数, 但是要建索引
        fai=rules.samtools_faidx.output,
    output:
        vcf="variants/{sample}.vcf",
    log:
        "logs/variants/{sample}.log",
    benchmark:
        "logs/variants/{sample}.bm"
    params:
        extra="-F 0.3 -C 10",
    threads: config["threads"]["low"]
    resources:
        mem_mb=1024,
    conda:
        config["conda"]["basic"]
    wrapper:
        f"file:{workflow.basedir}/wrappers/bio/freebayes"


rule vcf2table:
    input:
        rules.freebayes.output.vcf,
    output:
        "variants/{sample}.tsv",
    log:
        "logs/variants/{sample}.tsv.log",
    benchmark:
        "logs/variants/{sample}.tsv.bm"
    script:
        "../scripts/vcf2table.py"
