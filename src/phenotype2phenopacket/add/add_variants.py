import json
import re
from enum import Enum
from pathlib import Path

import click
import polars as pl
from google.protobuf.json_format import Parse
from phenopackets import (
    Diagnosis,
    Disease,
    Family,
    GeneDescriptor,
    GenomicInterpretation,
    Interpretation,
    OntologyClass,
    Phenopacket,
    VariantInterpretation,
    VariationDescriptor,
    VcfRecord,
)

from phenotype2phenopacket.utils.phenopacket_utils import write_phenopacket
from phenotype2phenopacket.utils.utils import all_files, read_variant_summary


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
        return mapping.get(clinical_significance, None)

    def filter_for_clinical_significance(self, filtered_disease_variant_summary):
        mapped_clinical_significance = filtered_disease_variant_summary.select(
            [
                pl.col("ClinicalSignificance").apply(lambda x: self.map_clinical_significance(x)),
                pl.exclude("ClinicalSignificance"),
            ],
            skip_nulls=True,
        )

        return mapped_clinical_significance.filter(
            pl.col("ClinicalSignificance") >= self.clinical_significance_filter
        )

    def filter_for_genome_assembly(self, filtered_clinical_significance_variant_summary):
        return filtered_clinical_significance_variant_summary.filter(
            pl.col("Assembly") == self.genome_assembly_filter
        )

    def filter_variant_summary(self):
        disease_variant_summary = self.filter_for_disease()
        if len(disease_variant_summary) == 0:
            print("No disease matches")
            return None
        filtered_variants = self.filter_for_clinical_significance(disease_variant_summary)
        return self.filter_for_genome_assembly(filtered_variants)


def return_phenopacket_disease(phenopacket: Phenopacket) -> [str]:
    return phenopacket.diseases[0]


def add_variant_genomic_interpretation(variant_entry: dict) -> GenomicInterpretation:
    return GenomicInterpretation(
        subject_or_biosample_id="patient1",
        interpretation_status=0,
        variant_interpretation=VariantInterpretation(
            acmg_pathogenicity_classification=variant_entry["ClinicalSignificance"],
            variation_descriptor=VariationDescriptor(
                id="clinvar:" + str(variant_entry["VariationID"]),
                gene_context=GeneDescriptor(
                    value_id=variant_entry["HGNC_ID"], symbol=variant_entry["GeneSymbol"]
                ),
                vcf_record=VcfRecord(
                    genome_assembly=variant_entry["Assembly"],
                    chrom=variant_entry["Chromosome"],
                    pos=variant_entry["Start"],
                    ref=variant_entry["ReferenceAllele"]
                    if variant_entry["ReferenceAllele"] != "na"
                    else variant_entry["ReferenceAlleleVCF"],
                    alt=variant_entry["AlternateAllele"]
                    if variant_entry["AlternateAllele"] != "na"
                    else variant_entry["AlternateAlleleVCF"],
                ),
            ),
        ),
    )


def add_variant_interpretations(
        disease: Disease,
        variant_summary_df: pl.DataFrame,
        clinical_significance_filter: int,
        genome_assembly_filter: str,
):
    filtered_variant_summary = VariantSummaryFilter(
        variant_summary_df=variant_summary_df,
        disease=disease,
        clinical_significance_filter=clinical_significance_filter,
        genome_assembly_filter=genome_assembly_filter,
    ).filter_variant_summary()
    if filtered_variant_summary is None:
        return None
    genomic_interpretations = []
    for row in filtered_variant_summary.rows(named=True):
        genomic_interpretations.append(add_variant_genomic_interpretation(row))
    return genomic_interpretations


def add_diagnosis(
        disease: Disease,
        variant_summary_df: pl.DataFrame,
        clinical_significance_filter: int,
        genome_assembly_filter: str,
) -> Diagnosis:
    genomic_interpretations = add_variant_interpretations(
        disease, variant_summary_df, clinical_significance_filter, genome_assembly_filter
    )
    if genomic_interpretations is None:
        return None
    return Diagnosis(
        disease=OntologyClass(
            id=disease.term.id,
            label=disease.term.label,
        ),
        genomic_interpretations=genomic_interpretations,
    )


def create_interpretation(
        phenopacket: Phenopacket,
        variant_summary_df: pl.DataFrame,
        clinical_significance_filter: int,
        genome_assembly_filter: str,
) -> Interpretation:
    disease = return_phenopacket_disease(phenopacket)
    diagnosis = add_diagnosis(
        disease, variant_summary_df, clinical_significance_filter, genome_assembly_filter
    )
    if diagnosis is None:
        return None
    return Interpretation(
        id=disease.term.label + "-interpretation",
        progress_status=0,
        diagnosis=diagnosis,
    )


def phenopacket_reader(file: Path):
    """Reads a phenopacket file, returning its contents."""
    file = open(file, "r")
    phenopacket = json.load(file)
    file.close()
    if "proband" in phenopacket:
        return Parse(json.dumps(phenopacket), Family())
    else:
        return Parse(json.dumps(phenopacket), Phenopacket())


def add_variant_interpretations_to_phenopacket(
        phenopacket: Phenopacket,
        variant_summary_df: pl.DataFrame,
        clinical_significance_filter: int,
        genome_assembly_filter: str,
):
    interpretation = create_interpretation(
        phenopacket, variant_summary_df, clinical_significance_filter, genome_assembly_filter
    )
    if interpretation is not None:
        phenopacket.interpretations.extend([interpretation])
        return phenopacket
    else:
        return None


def add_variants(
        phenopacket_path: Path,
        variant_summary_df: pl.DataFrame,
        clinical_significance_filter: int,
        genome_assembly_filter: str,
        output_dir: Path,
):
    phenopacket = phenopacket_reader(phenopacket_path)
    phenopacket_with_variants = add_variant_interpretations_to_phenopacket(
        phenopacket, variant_summary_df, clinical_significance_filter, genome_assembly_filter
    )
    if phenopacket_with_variants is not None:
        write_phenopacket(phenopacket_with_variants, output_dir.joinpath(phenopacket_path.name))


def add_variants_to_directory(
        phenopacket_dir: Path,
        variant_summary_df: pl.DataFrame,
        clinical_significance_filter: int,
        genome_assembly_filter: str,
        output_dir: Path,
):
    for file in all_files(phenopacket_dir):
        add_variants(
            file,
            variant_summary_df,
            clinical_significance_filter,
            genome_assembly_filter,
            output_dir,
        )


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
        clinical_significance_filter: str,
        genome_assembly_filter: str,
        output_dir: Path,
):
    output_dir.mkdir(exist_ok=True)
    variant_summary_df = read_variant_summary(variant_summary)
    add_variants_to_directory(
        phenopacket_dir,
        variant_summary_df,
        int(clinical_significance_filter),
        genome_assembly_filter,
        output_dir,
    )
