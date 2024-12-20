import sys

# logging
sys.stderr = open(snakemake.log[0], "w")

summary_files = snakemake.input
output_file = snakemake.output[0]

with open(output_file, 'w') as f:
    with open(summary_files[0], 'r') as g:
        f.write(g.read())

    for summary_file in summary_files[1:]:
        with open(summary_file, 'r') as g:
            f.write(g.readlines()[1])
