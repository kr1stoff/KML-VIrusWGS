import sys
import json

sys.stderr = open(snakemake.log[0], "w")

input_most_similar = snakemake.input.most_similar
acc_type = snakemake.input.acc_type
output = snakemake.output[0]

# 对照字典
acc_type_dict = json.load(open(acc_type))

# 读 accssion
with open(input_most_similar) as f:
    most_similar_acc = f.read().strip()

# 输出分型
with open(output, 'w') as f:
    f.write(acc_type_dict[most_similar_acc])
