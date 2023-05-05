from pathlib import Path

import click

from phenotype2phenopacket.add.add_genes import add_genes_to_directory
from phenotype2phenopacket.utils.utils import read_disease_pg


@click.command("add-genes")
@click.option(
    "--phenopacket-dir",
    "-p",
    required=True,
    help="Path to phenopacket directory.",
    type=Path,
)
@click.option(
    "--disease-pg",
    "-d",
    required=True,
    help="Path to disease.pg data file.",
    type=Path,
)
@click.option(
    "--hgnc-data",
    "-h",
    required=True,
    help="Path to hgnc_full_set data file.",
    type=Path,
)
@click.option(
    "--output-dir",
    "-o",
    required=True,
    help="Path to output directory.",
    type=Path,
)
def add_genes_command(
    phenopacket_dir: Path,
    disease_pg: Path,
    hgnc_data: Path,
    output_dir: Path,
):
    """Add known gene-to-phenotype relationships to a phenopacket."""
    output_dir.mkdir(exist_ok=True)
    disease_pg_df = read_disease_pg(disease_pg)
    add_genes_to_directory(
        phenopacket_dir,
        disease_pg_df,
        hgnc_data,
        output_dir,
    )
