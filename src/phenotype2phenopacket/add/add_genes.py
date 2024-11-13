from pathlib import Path

import polars as pl
from phenopackets import Disease
from pheval.utils.file_utils import all_files
from pheval.utils.phenopacket_utils import (
    GeneIdentifierUpdater,
    create_gene_identifier_map,
    create_hgnc_dict,
    phenopacket_reader,
)

from phenotype2phenopacket.utils.phenopacket_utils import (
    PhenopacketInterpretationExtender,
    PhenopacketUtil,
    write_phenopacket,
)


def get_phenotype_to_disease_entries(
    genes_to_disease: pl.DataFrame, disease: Disease
) -> pl.DataFrame:
    """
    Return disease.pg entries that match the provided OMIM disease ID.

    Args:
        genes_to_disease (pl.DataFrame): DataFrame containing genes_to_disease.txt entries.
        disease (Disease): Disease object containing the OMIM disease ID.

    Returns:
        pl.DataFrame: Filtered DataFrame containing entries matching the OMIM disease ID.
    """
    return genes_to_disease.filter(pl.col("disease_id") == disease.term.id)


def add_genes(
    phenopacket_path: Path,
    genes_to_disease: pl.DataFrame,
    gene_identifier_updater: GeneIdentifierUpdater,
    output_dir: Path,
):
    """
    Add known gene-to-phenotype relationships to the interpretations of a phenopacket.

    Args:
        phenopacket_path (Path): Path to the phenopacket file.
        genes_to_disease (pl.DataFrame): DataFrame containing genes_to_disease.txt entries.
        gene_identifier_updater (GeneIdentifierUpdater): Object for updating gene identifiers.
        output_dir (Path): Directory to write the updated phenopacket.

    """
    phenopacket = phenopacket_reader(phenopacket_path)
    disease = PhenopacketUtil(phenopacket).return_phenopacket_disease()
    filtered_genes_to_disease = get_phenotype_to_disease_entries(genes_to_disease, disease)
    if len(filtered_genes_to_disease) == 0:
        print(f"No gene-to-phenotype matches: {disease.term.id}, {disease.term.label}")
    else:
        phenopacket_with_genes = PhenopacketInterpretationExtender(
            phenopacket
        ).add_gene_interpretation_to_phenopacket(
            genes_to_disease_map=filtered_genes_to_disease,
            gene_identifier_updater=gene_identifier_updater,
        )
        (
            write_phenopacket(phenopacket_with_genes, output_dir.joinpath(phenopacket_path.name))
            if phenopacket_with_genes is not None
            else None
        )


def add_genes_to_directory(
    phenopacket_dir: Path,
    genes_to_disease: pl.DataFrame,
    gene_identifier: str,
    output_dir: Path,
):
    """
    Add known gene-to-phenotype relationships to the interpretations of a directory of phenopackets.

    Args:
        phenopacket_dir (Path): Directory containing the phenopacket files.
        genes_to_disease (pl.DataFrame): DataFrame containing genes_to_disease.txt entries.
        gene_identifier (str): Gene identifier for the phenopacket.
        output_dir (Path): Directory to store the updated phenopackets.
    """
    hgnc_dict = create_hgnc_dict()
    identifier_map = create_gene_identifier_map()
    gene_identifier_updater = GeneIdentifierUpdater(
        gene_identifier=gene_identifier, hgnc_data=hgnc_dict, identifier_map=identifier_map
    )
    for phenopacket_path in all_files(phenopacket_dir):
        add_genes(phenopacket_path, genes_to_disease, gene_identifier_updater, output_dir)
