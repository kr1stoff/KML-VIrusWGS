rule samtools_depth:
    input:
        bams=["align/{sample}.bam"],
    output:
        "depth/{sample}.depth",
    log:
        "logs/depth/{sample}.depth.log",
    params:
        # optional bed file passed to -b
        extra="-a",  # optional additional parameters as string
    conda:
        config["conda"]["basic"]
    wrapper:
        f"file:{workflow.basedir}/wrappers/bio/samtools/depth"


rule depth_stat:
    input:
        rules.samtools_depth.output,
    output:
        "depth/{sample}.depth.tsv",
    log:
        "logs/depth/{sample}.stat.log",
    benchmark:
        "logs/depth/{sample}.stat.bm"
    script:
        "../scripts/depth_stats.py"


rule depth_plot:
    input:
        rules.samtools_depth.output,
    output:
        "depth/{sample}.depth.png",
    log:
        "logs/depth/{sample}.plot.log",
    benchmark:
        "logs/depth/{sample}.plot.bm"
    script:
        "../scripts/depth_plot.py"


rule low_depth_region:
    input:
        "align/{sample}.bam",
    output:
        "depth/{sample}.low_depth.bed",
    log:
        "logs/depth/{sample}.low_depth.log",
    benchmark:
        "logs/depth/{sample}.low_depth.bm"
    params:
        min_depth="10",  # minimum depth
    conda:
        config["conda"]["basic"]
    shell:
        """
        bedtools genomecov -bga -ibam {input} \
            | awk '$4<{params.min_depth}' \
            | bedtools merge -i - \
            | awk '$3-$2>{params.min_depth}' \
            1> {output} 2> {log}
        """
