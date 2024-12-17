ruleorder: multiqc_pe > multiqc_se


rule multiqc_pe:
    input:
        expand(
            "qc/fastqc/{sample}.{pe}_fastqc.zip",
            sample=config["samples"],
            pe=["1", "2"],
        ),
    output:
        "qc/multiqc.html",
        directory("qc/multiqc_data"),
    params:
        extra="--verbose",
    log:
        "logs/multiqc.log",
    benchmark:
        "logs/multiqc.bm"
    conda:
        config["conda"]["basic"]
    wrapper:
        f"file:{workflow.basedir}/wrappers/bio/multiqc"


use rule multiqc_pe as multiqc_se with:
    input:
        expand(
            "qc/fastqc/{sample}_1_fastqc.zip",
            sample=config["samples"],
        ),
