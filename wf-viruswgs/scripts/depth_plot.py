import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import sys
matplotlib.use('Agg')


sys.stderr = open(snakemake.log[0], "w")


def plot_coverage_chart(in_depth, out_chart) -> None:
    # 输入数据
    df = pd.read_table(in_depth, sep='\t', header=None)
    df.columns = ['Ref', 'Pos', 'Depth']
    # 50个区间
    lbls = [f'bin{i}' for i in range(1, 51)]
    df['Bin'] = pd.cut(df['Pos'], 50, labels=lbls)
    # X轴刻度标签
    posrange = df['Pos'].max()-df['Pos'].min()
    if posrange < 1000000:
        xtkslbls = [f"{int(posrange/10/1000*i)}K" for i in range(10)]
    elif posrange < 1000000000:
        xtkslbls = [f"{int(posrange/10/1000000*i)}M" for i in range(10)]
    else:
        xtkslbls = [f"{int(posrange/10/1000000000*i)}G" for i in range(10)]
    dfgrp = df.groupby('Bin')
    dfgrpsum = dfgrp.sum(numeric_only=True)
    dfgrpmean = dfgrp.mean(numeric_only=True)
    dfgrpmedian = dfgrp.median(numeric_only=True)
    # 绘图
    fig = plt.figure(figsize=(10, 4))
    ax = plt.axes()
    plt.subplots_adjust(left=0.1, bottom=0.12, right=0.9, top=0.88)
    plt.bar(dfgrpsum.index, dfgrpsum['Depth'], color='skyblue')
    plt.xticks(np.arange(0, 50, step=5), labels=xtkslbls)
    # 建立次坐标轴
    ax2 = ax.twinx()
    ax2.set_ylim(0, max(dfgrpmean['Depth']))
    ax2.plot(dfgrpmean.index, dfgrpmean['Depth'], color='r', lw=0.8, label='Average Depth')
    ax2.plot(dfgrpmedian.index, dfgrpmedian['Depth'],
             color='yellow', lw=1, label='Median Depth', ls='--')
    # 将次坐标轴设置为红色
    ax2.tick_params(axis='y', colors='red')
    ax2.tick_params(axis='x', colors='w')
    ax2.spines['right'].set_color('red')
    ax2.yaxis.label.set_color('red')
    ax.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    # 其他属性
    # ax.ticklabel_format(style='plain', axis='y') #不要科学计数法
    ax.set_xlabel('Genomic position', fontsize=14)
    ax.set_ylabel('Mapped reads', fontsize=12)
    ax2.set_ylabel('Depth', fontsize=12)
    legend = plt.legend(loc=1)  # 图例的位置，1为右上角
    # 输出
    plt.savefig(out_chart)
    plt.cla()
    plt.close(fig)


in_tab = snakemake.input[0]
out_chart = snakemake.output[0]
plot_coverage_chart(in_tab, out_chart)
