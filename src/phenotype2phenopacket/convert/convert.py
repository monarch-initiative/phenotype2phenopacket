import random
from pathlib import Path

import polars as pl

from phenotype2phenopacket.utils.phenopacket_utils import (
    PhenotypeAnnotationToPhenopacketConverter,
    write_phenopacket,
)
from phenotype2phenopacket.utils.utils import load_ontology


def convert_to_phenopackets(
    phenotype_annotation: pl.DataFrame,
    num_disease: int,
    output_dir: Path,
):
    """Convert phenotype annotation file to a set of disease phenopackets."""
    human_phenotype_ontology = load_ontology()
    omim_diseases = phenotype_annotation.filter(pl.col("database_id").str.starts_with("OMIM"))
    grouped_omim_diseases = omim_diseases.partition_by(by="database_id", maintain_order=True)
    if num_disease != 0:
        grouped_omim_diseases = random.sample(grouped_omim_diseases, num_disease)
    for omim_disease in grouped_omim_diseases:
        phenopacket_file = PhenotypeAnnotationToPhenopacketConverter(
            human_phenotype_ontology
        ).create_phenopacket(omim_disease)
        write_phenopacket(
            phenopacket_file.phenopacket, output_dir.joinpath(phenopacket_file.phenopacket_path)
        )
