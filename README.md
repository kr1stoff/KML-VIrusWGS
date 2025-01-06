# KML-VirusWGS

病毒全基因组分析

## 命令行

```bash
poetry -C /data/mengxf/GitHub/KML-VirusWGS run python /data/mengxf/GitHub/KML-VirusWGS/main.py -s input.tsv -r /data/mengxf/Database/genome/HBV/BVBRC/BVBRC_genome_sequence.fasta -t /data/mengxf/Database/genome/HBV/BVBRC/accession_type_comparison.json -w result/241218
```

## 解释

1. 输入数据库文件

    例如 `HBV` 基因组, 从 `BVBRC` 勾选相关分型, 选择完整基因组下载基因组集. 测试流程中 `HBV` 共包含完整基因组 *283* 个

2. 输入分型对照表

    `json` 字典格式, 以 `accn` 作为键, 以分型作为值. 例如:

    ```json
    {
        "accn|LC753644": "C1",
        "accn|LC753643": "C1",
        "accn|LC753648": "C1",
        ....
    }
    ```
