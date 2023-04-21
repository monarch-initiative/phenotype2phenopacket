from pathlib import Path

import click

from phenotype2phenopacket.convert.convert import convert_to_phenopackets
from phenotype2phenopacket.utils.utils import read_phenotype_annotation_file


@click.command("convert")
@click.option(
    "--phenotype-annotation",
    "-p",
    required=True,
    help="Path to phenotype.hpoa.",
    type=Path,
)
@click.option(
    "--output-dir",
    "-o",
    required=True,
    help="Path to output directory.",
    type=Path,
    default="phenopackets",
    show_default=True,
)
def convert_to_phenopackets_command(phenotype_annotation: Path, output_dir: Path):
    output_dir.mkdir(exist_ok=True)
    phenotype_annotation_df = read_phenotype_annotation_file(phenotype_annotation)
    convert_to_phenopackets(phenotype_annotation_df, output_dir)
