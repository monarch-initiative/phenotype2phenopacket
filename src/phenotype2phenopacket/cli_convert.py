from pathlib import Path

import click

from phenotype2phenopacket.convert.convert import convert_to_phenopackets


@click.command("convert")
@click.option(
    "--phenotype-annotation",
    "-p",
    required=True,
    help="Path to phenotype.hpoa.",
    type=Path,
)
@click.option(
    "--num-disease",
    "-n",
    required=False,
    help="Number of diseases to create synthetic patient phenopackets for.",
    type=int,
    default=0,
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
def convert_to_phenopackets_command(phenotype_annotation: Path, num_disease: int, output_dir: Path):
    """Convert a phenotype annotation file (phenotype.hpoa) to a set of disease phenopackets."""
    output_dir.mkdir(exist_ok=True)
    convert_to_phenopackets(phenotype_annotation, num_disease, output_dir)
