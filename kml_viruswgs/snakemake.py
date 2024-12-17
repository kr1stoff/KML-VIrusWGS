import yaml
from pathlib import Path
from subprocess import run
import logging
from kml_viruswgs import get_conda_env_dict
from kml_viruswgs import get_threads_dict
from kml_viruswgs import get_my_scripts_path
from kml_viruswgs import get_database_dict
from kml_viruswgs import get_software_dict


def create_snakemake_configfile(sample_names, workdir):
    """
    创建 snakemake 配置文件
    :param sample_names:    样本名列表
    :param workdir:         分析结果目录
    :return:
    """
    logging.info('创建 snakemake 配置文件')
    workdir = str(Path(workdir).resolve())
    dir_temp = Path(workdir).joinpath('.temp')
    dir_temp.mkdir(exist_ok=True, parents=True)

    dict_smk = {
        'workdir': workdir,
        'samples': sample_names,
        'threads': get_threads_dict(),
        'conda': get_conda_env_dict(),
        'my_scripts': get_my_scripts_path(),
        'database': get_database_dict()
    }

    with open(f'{workdir}/.temp/snakemake.yaml', 'w') as f:
        yaml.dump(dict_smk, f)


def run_snakemake(workdir):
    """
    运行 snakemake 工作流
    :param workdir:
    :return:
    """
    logging.info('运行 snakemake')
    activate = get_software_dict()['activate']
    cores = get_threads_dict()['high']
    snakefile = Path(__file__).resolve().parents[1].joinpath('wf-lvisa/Snakefile')
    configfile = f'{workdir}/.temp/snakemake.yaml'
    logfile = f'{workdir}/.log/snakemake.log'

    cml = f"""
    source {activate} snakemake
    # use-conda
    snakemake -c {cores} --use-conda -s {snakefile} --configfile {configfile}
    """

    proc = run(cml, shell=True, executable='/bin/bash', capture_output=True, encoding='utf-8')

    # 输出出来这段日志
    with open(logfile, 'w') as f:
        f.write(f'[STDOUT]\n{proc.stdout}\n[STDERR]\n{proc.stderr}')
