import sys

sys.stderr = open(snakemake.log[0], "w")


def parse_vcfinfo(info):
    """return INFO dict. DP,AC,AF..."""
    info_dict = dict()
    info_list = info.strip().split(";")
    for il in info_list:
        il_list = il.strip().split("=")
        info_dict[il_list[0]] = il_list[1]
    return info_dict


def vcf2table_normal(vfile, outfile):
    """相较于vcf2table函数, 处理未注释的VCF文件"""
    with open(vfile) as fh, open(outfile, "wt", encoding="utf-8", newline="") as gh:
        gh.write("参考基因组\t变异位置\t参考碱基\t替换碱基\t变异深度\n")
        for line in fh:
            if line.startswith("##"):
                continue
            elif line.startswith("#CHROM"):
                header_list = line.strip().split("\t")
                header_dict = {hie[1]: hie[0] for hie in enumerate(header_list)}
            else:
                linelist = line.strip().split("\t")
                chrom = linelist[header_dict["#CHROM"]]
                pos = linelist[header_dict["POS"]]
                ref = linelist[header_dict["REF"]]
                alt = linelist[header_dict["ALT"]]
                qual = float(linelist[header_dict["QUAL"]])
                info = linelist[header_dict["INFO"]]
                info_dict = parse_vcfinfo(info)
                if "DP" not in info_dict:
                    sys.stderr.write("WARNING - {}".format(info))
                    continue
                elif qual < 1:  # 过滤掉质量小于1
                    continue
                else:
                    depth = format(int(info_dict["DP"]), ",")
                    outlist = [chrom, pos, ref, alt] + [depth]
                    gh.write("\t".join(outlist) + "\n")


vcf_file = snakemake.input[0]
out_tsv = snakemake.output[0]
vcf2table_normal(vcf_file, out_tsv)
