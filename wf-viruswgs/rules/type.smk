rule virus_type:
    input:
        most_similar=rules.find_ref.output.most_similar,
        acc_type=config["database"]["acc_type"],
    output:
        "type/{sample}.type",
    log:
        "logs/type/{sample}.virus_type.log",
    benchmark:
        "logs/type/{sample}.virus_type.bm"
    script:
        "../scripts/virus_type.py"
