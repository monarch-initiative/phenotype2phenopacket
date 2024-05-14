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
    omim_disease_pg: pl.DataFrame, disease: Disease
) -> pl.DataFrame:
    """
    Return disease.pg entries that match the provided OMIM disease ID.

    Args:
        omim_disease_pg (pl.DataFrame): DataFrame containing disease.pg entries.
        disease (Disease): Disease object containing the OMIM disease ID.

    Returns:
        pl.DataFrame: Filtered DataFrame containing entries matching the OMIM disease ID.
    """
    return omim_disease_pg.filter(pl.col("database_id") == disease.term.id)


def add_genes(
    phenopacket_path: Path,
    disease_pg: pl.DataFrame,
    gene_identifier_updater: GeneIdentifierUpdater,
    output_dir: Path,
):
    """
    Add known gene-to-phenotype relationships to the interpretations of a phenopacket.

    Args:
        phenopacket_path (Path): Path to the phenopacket file.
        disease_pg (pl.DataFrame): DataFrame containing disease.pg entries.
        gene_identifier_updater (GeneIdentifierUpdater): Object for updating gene identifiers.
        output_dir (Path): Directory to write the updated phenopacket.

    """
    phenopacket = phenopacket_reader(phenopacket_path)
    disease = PhenopacketUtil(phenopacket).return_phenopacket_disease()
    filtered_disease_pg = get_phenotype_to_disease_entries(disease_pg, disease)
    if len(filtered_disease_pg) == 0:
        print(f"No gene-to-phenotype matches: {disease.term.id}, {disease.term.label}")
    else:
        phenopacket_with_genes = PhenopacketInterpretationExtender(
            phenopacket
        ).add_gene_interpretation_to_phenopacket(
            omim_disease_phenotype_gene_map=filtered_disease_pg,
            gene_identifier_updater=gene_identifier_updater,
        )
        (
            write_phenopacket(phenopacket_with_genes, output_dir.joinpath(phenopacket_path.name))
            if phenopacket_with_genes is not None
            else None
        )


def add_genes_to_directory(phenopacket_dir: Path, disease_pg: pl.DataFrame, output_dir: Path):
    """
    Add known gene-to-phenotype relationships to the interpretations of a directory of phenopackets.

    Args:
        phenopacket_dir (Path): Directory containing the phenopacket files.
        disease_pg (pl.DataFrame): DataFrame containing disease.pg entries.
        output_dir (Path): Directory to store the updated phenopackets.
    """
    hgnc_dict = create_hgnc_dict()
    identifier_map = create_gene_identifier_map()
    gene_identifier_updater = GeneIdentifierUpdater(
        gene_identifier="ensembl_id", hgnc_data=hgnc_dict, identifier_map=identifier_map
    )
    for phenopacket_path in all_files(phenopacket_dir):
        add_genes(phenopacket_path, disease_pg, gene_identifier_updater, output_dir)
