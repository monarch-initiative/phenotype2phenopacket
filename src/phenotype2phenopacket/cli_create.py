from pathlib import Path

import click

from phenotype2phenopacket.create.create import create_synthetic_patients


@click.command("create")
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
def create_synthetic_patient_command(
    phenotype_annotation: Path, num_disease: int, output_dir: Path
):
    """
    Create a set of synthetic patient phenopackets from a phenotype annotation file.

    Args:
        phenotype_annotation (Path): Path to the phenotype annotation file.
        num_disease (int): Number of diseases to create synthetic patient phenopackets (use 0 for all).
        output_dir (Path): Directory to store the generated synthetic patient phenopackets.
    """
    output_dir.mkdir(exist_ok=True)
    create_synthetic_patients(phenotype_annotation, num_disease, output_dir)
