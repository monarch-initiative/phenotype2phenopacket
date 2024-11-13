from pathlib import Path

import click

from phenotype2phenopacket.add.add_genes import add_genes_to_directory
from phenotype2phenopacket.utils.utils import read_genes_to_disease


@click.command("add-genes")
@click.option(
    "--phenopacket-dir",
    "-p",
    required=True,
    help="Path to phenopacket directory.",
    type=Path,
)
@click.option(
    "--genes-to-disease",
    "-g",
    required=True,
    help="Path to genes_to_disease.txt data file.",
    type=Path,
)
@click.option(
    "--output-dir",
    "-o",
    required=True,
    help="Path to output directory.",
    type=Path,
)
@click.option(
    "--gene-identifier",
    "-i",
    required=False,
    help="Gene identifier to update in phenopacket",
    type=click.Choice(["ensembl_id", "entrez_id", "hgnc_id"]),
)
def add_genes_command(
    phenopacket_dir: Path,
    genes_to_disease: Path,
    gene_identifier: str,
    output_dir: Path,
):
    """
    Add known gene-to-phenotype relationships to a set of phenopackets in a directory.

    Args:
        phenopacket_dir (Path): Directory containing the phenopacket files.
        genes_to_disease (Path): Path to the genes_to_disease.txt file.
        gene_identifier (str): Gene identifier to add in phenopacket
        output_dir (Path): Directory to store the updated phenopackets.
    """
    output_dir.mkdir(exist_ok=True)
    genes_to_disease_df = read_genes_to_disease(genes_to_disease)
    add_genes_to_directory(
        phenopacket_dir,
        genes_to_disease_df,
        gene_identifier,
        output_dir,
    )
