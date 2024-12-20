import logging
import click
from pathlib import Path
from kml_viruswgs import prepare_fastq_by_samptab
from kml_viruswgs import get_sample_names_by_samptab
from kml_viruswgs import create_snakemake_configfile
from kml_viruswgs import run_snakemake


logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s - %(levelname)s - %(filename)s - %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")


@click.command()
@click.option('--sample_table', '-s', type=click.Path(exists=True), required=True, help='样本信息表.')
@click.option('--database', '-r', type=click.Path(exists=True), required=True, help='参考基因组集 FASTA, 比如 HBV 的几百个基因组.')
@click.option('--acc_type', '-t', type=click.Path(exists=True), required=True,
              help='参考基因组集的分型对照表, 里面记载 database 参数中所有基因组登录号对应的分型.')
@click.option('--work_dir', '-w', type=str, default='viruswgs_result', help='结果生成目录. [default: viruswgs_result]')
@click.help_option('-h', '--help')
def main(work_dir, sample_table, database, acc_type):
    """通用病毒全基因组分析流程"""
    logging.info(f'开始分析!')
    sample_table = Path(sample_table).resolve()
    work_dir = Path(work_dir).resolve()

    # fastq
    sample_names = get_sample_names_by_samptab(sample_table)
    prepare_fastq_by_samptab(work_dir, sample_table)

    # todo snakemake
    create_snakemake_configfile(sample_names, work_dir, database, acc_type)
    run_snakemake(work_dir)

    logging.info(f'分析完成!')


if __name__ == '__main__':
    main()
