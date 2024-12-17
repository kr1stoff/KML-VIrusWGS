# bwa 方式
bwa mem -t 16 -M -Y -R '@RG\tID:test\tSM:test' \
    /data/mengxf/Database/genome/HBV/BVBRC/BVBRC_genome_sequence.fasta \
    ../../rawdata/BD315.fastq.gz |
    samtools view -@ 4 -hbS - |
    samtools sort -@ 4 -o test.bam -
samtools coverage test.bam >test.coverage
csvtk -t -C '' sort -k coverage:nr -k meanmapq:nr test.coverage | sed '1d' | head -n1 | cut -f1 >most_similar.id
seqtk subseq /data/mengxf/Database/genome/HBV/BVBRC/BVBRC_genome_sequence.fasta most_similar.id >ref.fa
bwa index ref.fa
bwa mem -t 16 -M -Y -R '@RG\tID:test\tSM:test' ref.fa ../../rawdata/BD315.fastq.gz |
    samtools view -@ 4 -hbS - |
    samtools sort -@ 4 -o test2.bam -
samtools coverage test2.bam >test2.coverage
