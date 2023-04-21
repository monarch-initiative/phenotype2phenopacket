import re
from enum import Enum
from pathlib import Path

import polars as pl
from phenopackets import Disease

from phenotype2phenopacket.utils.gene_map_utils import GeneIdentifierUpdater, create_hgnc_dict
from phenotype2phenopacket.utils.phenopacket_utils import (
    PhenopacketInterpretationExtender,
    PhenopacketUtil,
    phenopacket_reader,
    write_phenopacket,
)
from phenotype2phenopacket.utils.utils import all_files


class ClinicalSignificance(Enum):
    PATHOGENIC = 5
    LIKELY_PATHOGENIC = 4
    UNCERTAIN_SIGNIFICANCE = 3
    LIKELY_BENIGN = 2
    BENIGN = 1
    NOT_PROVIDED = 0


class VariantSummaryFilter:
    def __init__(
        self,
        variant_summary_df: pl.DataFrame,
        disease: Disease,
        clinical_significance_filter: int,
        genome_assembly_filter: str,
    ):
        self.variant_summary_df = variant_summary_df
        self.disease = disease
        self.clinical_significance_filter = clinical_significance_filter
        self.genome_assembly_filter = genome_assembly_filter

    def filter_for_disease(self) -> pl.DataFrame:
        return self.variant_summary_df.filter(
            pl.col("PhenotypeIDS").str.contains(r"\b" + re.escape(self.disease.term.id) + r"\b")
        )

    @staticmethod
    def map_clinical_significance(clinical_significance) -> int:
        mapping = {
            "Pathogenic": ClinicalSignificance.PATHOGENIC.value,
            "Likely pathogenic": ClinicalSignificance.LIKELY_PATHOGENIC.value,
            "Uncertain significance": ClinicalSignificance.UNCERTAIN_SIGNIFICANCE.value,
            "Likely benign": ClinicalSignificance.LIKELY_BENIGN.value,
            "Benign": ClinicalSignificance.BENIGN.value,
            "Not provided": ClinicalSignificance.NOT_PROVIDED.value,
        }
        return mapping.get(clinical_significance, -1)

    def filter_for_clinical_significance(self, filtered_disease_variant_summary):
        mapped_clinical_significance = filtered_disease_variant_summary.select(
            [
                pl.col("ClinicalSignificance").apply(lambda x: self.map_clinical_significance(x)),
                pl.exclude("ClinicalSignificance"),
            ],
            skip_nulls=True,
        )
        removed_unknown_clinical_significance = mapped_clinical_significance.filter(
            pl.col("ClinicalSignificance") is not None
        )
        return removed_unknown_clinical_significance.filter(
            pl.col("ClinicalSignificance") >= self.clinical_significance_filter
        )

    def filter_for_genome_assembly(self, filtered_clinical_significance_variant_summary):
        return filtered_clinical_significance_variant_summary.filter(
            pl.col("Assembly") == self.genome_assembly_filter
        )

    def filter_variant_summary(self):
        disease_variant_summary = self.filter_for_disease()
        if len(disease_variant_summary) == 0:
            return None
        filtered_variants = self.filter_for_clinical_significance(disease_variant_summary)
        return self.filter_for_genome_assembly(filtered_variants)


def add_variants(
    phenopacket_path: Path,
    variant_summary_df: pl.DataFrame,
    clinical_significance_filter: int,
    genome_assembly_filter: str,
    output_dir: Path,
    gene_identifier_updater,
):
    phenopacket = phenopacket_reader(phenopacket_path)
    disease = PhenopacketUtil(phenopacket).return_phenopacket_disease()
    filtered_variant_summary = VariantSummaryFilter(
        variant_summary_df=variant_summary_df,
        disease=disease,
        clinical_significance_filter=clinical_significance_filter,
        genome_assembly_filter=genome_assembly_filter,
    ).filter_variant_summary()
    if filtered_variant_summary is None:
        print("No disease matches")
    else:
        phenopacket_with_variants = PhenopacketInterpretationExtender(
            phenopacket
        ).add_variant_interpretation_to_phenopacket(
            filtered_variant_summary, gene_identifier_updater
        )
        write_phenopacket(phenopacket_with_variants, output_dir.joinpath(phenopacket_path.name))


def add_variants_to_directory(
    phenopacket_dir: Path,
    variant_summary_df: pl.DataFrame,
    hgnc_data: Path,
    clinical_significance_filter: int,
    genome_assembly_filter: str,
    output_dir: Path,
):
    gene_identifier_updater = GeneIdentifierUpdater(
        gene_identifier="ensembl_id", hgnc_data=create_hgnc_dict(hgnc_data)
    )
    for file in all_files(phenopacket_dir):
        add_variants(
            file,
            variant_summary_df,
            clinical_significance_filter,
            genome_assembly_filter,
            output_dir,
            gene_identifier_updater,
        )
