from pathlib import Path

import click

from phenotype2phenopacket.add.add_genes import add_genes_to_directory
from phenotype2phenopacket.add.add_variants import add_variants_to_directory
from phenotype2phenopacket.utils.utils import read_disease_pg, read_variant_summary


@click.command("add-variants")
@click.option(
    "--phenopacket-dir",
    "-p",
    required=True,
    help="Path to phenopacket directory.",
    type=Path,
)
@click.option(
    "--variant-summary",
    "-v",
    required=True,
    help="Path to variant summary data file.",
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
    "--clinical-significance-filter",
    "-c",
    required=True,
    help="Clinical significance filter for variants. 0:NOT_PROVIDED, 1:BENIGN, 2:LIKELY_BENIGN, "
    "3:UNCERTAIN_SIGNIFICANCE, 4:LIKELY_PATHOGENIC, 5:PATHOGENIC ",
    type=click.Choice(["0", "1", "2", "3", "4", "5"], case_sensitive=False),
)
@click.option(
    "--genome-assembly-filter",
    "-g",
    required=True,
    help="Genome assembly filter",
    type=click.Choice(["GRCh37", "GRCh38"], case_sensitive=True),
)
@click.option(
    "--output-dir",
    "-o",
    required=True,
    help="Path to output directory.",
    type=Path,
)
def add_variants_command(
    phenopacket_dir: Path,
    variant_summary: Path,
    hgnc_data: Path,
    clinical_significance_filter: str,
    genome_assembly_filter: str,
    output_dir: Path,
):
    output_dir.mkdir(exist_ok=True)
    variant_summary_df = read_variant_summary(variant_summary)
    add_variants_to_directory(
        phenopacket_dir,
        variant_summary_df,
        hgnc_data,
        int(clinical_significance_filter),
        genome_assembly_filter,
        output_dir,
    )


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
    output_dir.mkdir(exist_ok=True)
    disease_pg_df = read_disease_pg(disease_pg)
    add_genes_to_directory(
        phenopacket_dir,
        disease_pg_df,
        hgnc_data,
        output_dir,
    )
