from kml_viruswgs import prepare_fastq_by_samptab
from kml_viruswgs import get_sample_names_by_samptab

work_dir = '/data/mengxf/Project/KML241212_HBV_WGS/result/241217'
sample_table = '/data/mengxf/GitHub/KML-VirusWGS/templates/input.tsv'


def test_prepare():
    prepare_fastq_by_samptab(work_dir, sample_table)


def test_get():
    names = get_sample_names_by_samptab(sample_table)
    print(names)
