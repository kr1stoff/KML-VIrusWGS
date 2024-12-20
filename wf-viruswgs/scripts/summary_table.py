import sys

# logging
sys.stderr = open(snakemake.log[0], "w")


print(snakemake.input)
depth_file = snakemake.input[0]
vtype_file = snakemake.input[1]
name = snakemake.params.name
output_file = snakemake.output[0]

# 深度
with open(depth_file, 'r') as f:
    depth_header = f.readline().strip().split('\t')
    depth_data = f.readline().strip().split('\t')

# 分型
with open(vtype_file, 'r') as f:
    vtype = f.read().strip()

# 输出
with open(output_file, 'w') as f:
    f.write('\t'.join(['Sample', 'Type'] + depth_header) + '\n')
    f.write('\t'.join([name, vtype] + depth_data) + '\n')
