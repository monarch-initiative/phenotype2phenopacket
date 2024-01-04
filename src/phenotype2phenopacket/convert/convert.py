import random
from pathlib import Path

import polars as pl

from phenotype2phenopacket.utils.phenopacket_utils import (
    PhenotypeAnnotationToPhenopacketConverter,
    write_phenopacket,
)
from phenotype2phenopacket.utils.utils import (
    get_phenotype_annotation_version,
    load_ontology,
    read_phenotype_annotation_file,
)


def convert_to_phenopackets(
    phenotype_annotation: Path,
    num_disease: int,
    output_dir: Path,
):
    """
    Convert a phenotype annotation file to a set of disease-specific phenopackets.

    Args:
        phenotype_annotation (Path): Path to the phenotype annotation file.
        num_disease (int): Number of diseases to convert to phenopackets.
                           If set to 0, processes all available diseases.
        output_dir (Path): Directory path to write the generated phenopackets.
    """
    human_phenotype_ontology = load_ontology()
    phenotype_annotation_df = read_phenotype_annotation_file(phenotype_annotation)
    phenotype_annotation_version = get_phenotype_annotation_version(phenotype_annotation)
    omim_diseases = phenotype_annotation_df.filter(pl.col("database_id").str.starts_with("OMIM"))
    grouped_omim_diseases = omim_diseases.partition_by(by="database_id", maintain_order=True)
    if num_disease != 0:
        grouped_omim_diseases = random.sample(grouped_omim_diseases, num_disease)
    for omim_disease in grouped_omim_diseases:
        phenopacket_file = PhenotypeAnnotationToPhenopacketConverter(
            human_phenotype_ontology
        ).create_phenopacket(omim_disease, phenotype_annotation_version)
        write_phenopacket(
            phenopacket_file.phenopacket, output_dir.joinpath(phenopacket_file.phenopacket_path)
        )
