from pathlib import Path

import click
from pheval.prepare.custom_exceptions import MutuallyExclusiveOptionError

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
    cls=MutuallyExclusiveOptionError,
    mutually_exclusive=["omim_id_list"],
)
@click.option(
    "--omim-id",
    "-i",
    required=False,
    help="OMIM ID to create synthetic patient for",
    type=str,
    default=None,
    cls=MutuallyExclusiveOptionError,
    mutually_exclusive=["omim_id_list"],
)
@click.option(
    "--omim-id-list",
    "-l",
    required=False,
    help="Path to .txt file containing OMIM IDs to create synthetic patient phenopackets,"
    "with each OMIM ID separated by a new line.",
    type=Path,
    default=None,
    cls=MutuallyExclusiveOptionError,
    mutually_exclusive=["omim_id", "num_disease"],
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
@click.option(
    "--local-ontology-cache",
    "-c",
    metavar="PATH",
    required=False,
    help="Path to the local ontology cache, e.g., path to the hp.obo.",
    default=None,
    type=Path,
)
def convert_to_phenopackets_command(
    phenotype_annotation: Path,
    num_disease: int,
    omim_id: str,
    omim_id_list: Path,
    output_dir: Path,
    local_ontology_cache: Path,
):
    """
    Convert a phenotype annotation file to a set of disease phenopackets.

    Args:
        phenotype_annotation (Path): Path to the phenotype annotation file.
        num_disease (int): Number of diseases to create phenopackets (use 0 for all).
        omim_id (str): OMIM ID to create  phenopacket for
        omim_id_list (Path): Path to the text file containing OMIM IDs to create synthetic patient phenopackets.
        output_dir (Path): Directory to store the generated phenopackets.
        local_ontology_cache (Path): Path to the local ontology cache.
    """
    output_dir.mkdir(exist_ok=True)
    convert_to_phenopackets(
        phenotype_annotation, num_disease, omim_id, omim_id_list, output_dir, local_ontology_cache
    )
