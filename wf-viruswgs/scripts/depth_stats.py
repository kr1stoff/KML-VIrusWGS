import pandas as pd
import sys

# logging
sys.stderr = open(snakemake.log[0], "w")


def get_bam_stat(in_depth, out_stat) -> None:
    # 分类数据
    df = pd.read_table(in_depth, sep='\t', header=None)
    df.columns = ['Ref', 'Pos', 'Depth']
    # 深度统计表格
    meandepth = '{:.1f}'.format(df['Depth'].mean())
    coverage0 = '{:.2%}'.format(len(df[df['Depth'] > 0])/df.shape[0])
    coverage10 = '{:.2%}'.format(len(df[df['Depth'] > 10])/df.shape[0])
    coverage30 = '{:.2%}'.format(len(df[df['Depth'] > 30])/df.shape[0])
    coverage100 = '{:.2%}'.format(len(df[df['Depth'] > 100])/df.shape[0])
    uniformity = '{:.2%}'.format(len(df[df['Depth'] > df['Depth'].mean()*0.2])/df.shape[0])
    with open(out_stat, 'w') as g:
        g.write('MeanDepth\tCoverage\tDepth ≥10x\tDepth ≥30x\tDepth ≥100x\tUniformity\n')
        g.write('\t'.join([meandepth, coverage0, coverage10,
                coverage30, coverage100, uniformity]) + '\n')


in_depth = snakemake.input[0]
out_stat = snakemake.output[0]
get_bam_stat(in_depth, out_stat)
