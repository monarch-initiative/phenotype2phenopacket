import json
import re
import secrets
from copy import copy
from dataclasses import dataclass
from pathlib import Path
from typing import Union

import polars as pl
from google.protobuf.json_format import MessageToJson, Parse
from google.protobuf.timestamp_pb2 import Timestamp
from phenopackets import (
    Age,
    Diagnosis,
    Disease,
    Family,
    GeneDescriptor,
    GenomicInterpretation,
    Individual,
    Interpretation,
    MetaData,
    OntologyClass,
    Phenopacket,
    PhenotypicFeature,
    Resource,
    TimeElement,
    VariantInterpretation,
    VariationDescriptor,
    VcfRecord,
)

from phenotype2phenopacket.utils.gene_map_utils import GeneIdentifierUpdater


@dataclass
class OnsetTerm:
    lower_age: Union[int, str]
    upper_age: Union[int, str]


onset_hpo = {
    "HP:0011462": OnsetTerm(16, 40),
    "HP:0011460": OnsetTerm(0, 0.019230769230769232),
    "HP:0011463": OnsetTerm(1, 5),
    "HP:0003584": OnsetTerm(60, 90),
    "HP:0025709": OnsetTerm(19, 25),
    "HP:0034198": OnsetTerm(0, 0.019230769230769232),
    "HP:0003596": OnsetTerm(40, 60),
    "HP:0003621": OnsetTerm(5, 15),
    "HP:0003593": OnsetTerm(0, 1),
    "HP:4000040": OnsetTerm(0, 0),
    "HP:0011461": OnsetTerm(0, 0.019230769230769232),
    "HP:0003577": OnsetTerm(0, 0),
    "HP:0034197": OnsetTerm(0, 0.019230769230769232),
    "HP:0025708": OnsetTerm(16, 19),
    "HP:0410280": OnsetTerm(1, 15),
    "HP:0030674": OnsetTerm(0, 0.019230769230769232),
    "HP:0003623": OnsetTerm(0, 0.019230769230769232),
    "HP:0034199": OnsetTerm(0, 0.019230769230769232),
    "HP:0003581": OnsetTerm(16, 80),
    "HP:0025710": OnsetTerm(25, 40),
}


@dataclass
class PhenopacketFile:
    phenopacket: Phenopacket
    phenopacket_path: Path


def phenopacket_reader(file: Path):
    """Reads a phenopacket file, returning its contents."""
    file = open(file, "r")
    phenopacket = json.load(file)
    file.close()
    if "proband" in phenopacket:
        return Parse(json.dumps(phenopacket), Family())
    else:
        return Parse(json.dumps(phenopacket), Phenopacket())


def create_phenopacket_file_name_from_disease(disease_name: str) -> Path:
    normalised_string = re.sub(r"\W+", "_", disease_name)
    return Path(normalised_string.replace(" ", "_") + ".json")


def create_json_message(phenopacket: Phenopacket) -> str:
    """Create json message for writing to file."""
    return MessageToJson(phenopacket)


def write_phenopacket(phenopacket: Phenopacket, output_file: Path) -> None:
    """Write a phenopacket."""
    phenopacket_json = create_json_message(phenopacket)
    with open(output_file, "w") as outfile:
        outfile.write(phenopacket_json)
    outfile.close()


class PhenotypeAnnotationToPhenopacketConverter:
    def __init__(self, human_phenotype_ontology):
        self.human_phenotype_ontology = human_phenotype_ontology

    @staticmethod
    def create_individual(onset_range: OnsetTerm = None) -> Individual:
        secrets_random_num = secrets.SystemRandom()
        """Create an Individual object."""
        age = (
            secrets_random_num.randint(onset_range.lower_age, onset_range.upper_age)
            if onset_range is not None
            else None
        )
        if onset_range.upper_age == 0 and onset_range.lower_age == 0:
            age = None
        return Individual(
            id="patient1",
            time_at_last_encounter=TimeElement(age=Age(iso8601duration=f"P{age}Y"))
            if age is not None
            else None,
        )

    def create_onset(self, phenotype_annotation_entry: dict) -> TimeElement:
        """Create an Onset object."""
        if phenotype_annotation_entry["onset"] is not None:
            rels = self.human_phenotype_ontology.entity_alias_map(
                phenotype_annotation_entry["onset"]
            )
            term = "".join(rels[(list(rels.keys())[0])])
            return TimeElement(
                ontology_class=OntologyClass(id=phenotype_annotation_entry["onset"], label=term)
            )
        else:
            return None

    def create_modifier(self, phenotype_annotation_entry: dict) -> [OntologyClass]:
        """Create a Modifier."""
        if phenotype_annotation_entry["modifier"] is not None:
            try:
                rels = self.human_phenotype_ontology.entity_alias_map(
                    phenotype_annotation_entry["modifier"]
                )
                term = "".join(rels[(list(rels.keys())[0])])
                return [OntologyClass(id=phenotype_annotation_entry["modifier"], label=term)]
            except IndexError:
                return [OntologyClass(id=phenotype_annotation_entry["modifier"])]
        else:
            return None

    def create_phenotypic_feature(self, phenotype_annotation_entry: dict) -> PhenotypicFeature:
        """Create a PhenotypicFeature object."""
        if phenotype_annotation_entry["aspect"] == "P":
            rels = self.human_phenotype_ontology.entity_alias_map(
                phenotype_annotation_entry["hpo_id"]
            )
            hpo_term = "".join(rels[(list(rels.keys())[0])])
            return PhenotypicFeature(
                type=OntologyClass(id=phenotype_annotation_entry["hpo_id"], label=hpo_term),
                onset=self.create_onset(phenotype_annotation_entry),
                modifiers=self.create_modifier(phenotype_annotation_entry),
                excluded=True if phenotype_annotation_entry["qualifier"] == "NOT" else False,
            )
        else:
            return None

    @staticmethod
    def create_disease(phenotype_annotation_entry: dict) -> Disease:
        """Create a Disease object."""
        return Disease(
            term=OntologyClass(
                id=phenotype_annotation_entry["database_id"],
                label=phenotype_annotation_entry["disease_name"],
            )
        )

    @staticmethod
    def create_omim_resource() -> Resource:
        """Create OMIM resource."""
        return Resource(
            id="omim",
            name="Online Mendelian Inheritance in Man",
            url="https://www.omim.org",
            version="hp/releases/2023-04-18",
            namespace_prefix="OMIM",
            iri_prefix="https://omim.org/entry/",
        )

    @staticmethod
    def create_human_phenotype_ontology_resource() -> Resource:
        """Create human phenotype ontology resource."""
        return Resource(
            id="hp",
            name="human phenotype ontology",
            url="http://purl.obolibrary.org/obo/hp.owl",
            version="hp/releases/2023-04-05",
            namespace_prefix="HP",
            iri_prefix="http://purl.obolibrary.org/obo/HP_",
        )

    def create_metadata(self) -> MetaData:
        """Create metadata"""
        timestamp = Timestamp()
        timestamp.GetCurrentTime()
        return MetaData(
            created=timestamp,
            created_by="phenotype2phenopacket",
            resources=[
                self.create_human_phenotype_ontology_resource(),
                self.create_omim_resource(),
            ],
            phenopacket_schema_version="2.0",
        )

    def create_phenopacket(
        self, omim_disease_df: pl.DataFrame, onset: OnsetTerm = None
    ) -> PhenopacketFile:
        phenotypic_features = []
        phenotype_annotation_entry = ""
        for phenotype_annotation_entry in omim_disease_df.rows(named=True):
            phenotype_annotation_entry = phenotype_annotation_entry
            phenotypic_feature = self.create_phenotypic_feature(phenotype_annotation_entry)
            if phenotypic_feature is not None:
                phenotypic_features.append(phenotypic_feature)
        return PhenopacketFile(
            phenopacket=Phenopacket(
                id=phenotype_annotation_entry["disease_name"].lower().replace(" ", "_"),
                subject=self.create_individual(onset),
                phenotypic_features=phenotypic_features,
                diseases=[self.create_disease(phenotype_annotation_entry)],
                meta_data=self.create_metadata(),
            ),
            phenopacket_path=create_phenopacket_file_name_from_disease(
                phenotype_annotation_entry["disease_name"]
            ),
        )


class PhenopacketUtil:
    def __init__(self, phenopacket: Phenopacket):
        self.phenopacket = phenopacket

    def return_phenopacket_disease(self) -> [str]:
        return self.phenopacket.diseases[0]


class PhenopacketInterpretationExtender:
    def __init__(self, phenopacket: Phenopacket):
        self.phenopacket = phenopacket

    @staticmethod
    def add_variant_genomic_interpretation(variant_entry: dict, gene_identifier_updater):
        try:
            return GenomicInterpretation(
                subject_or_biosample_id="patient1",
                interpretation_status=0,
                variant_interpretation=VariantInterpretation(
                    acmg_pathogenicity_classification=variant_entry["ClinicalSignificance"],
                    variation_descriptor=VariationDescriptor(
                        id="clinvar:" + str(variant_entry["VariationID"]),
                        gene_context=GeneDescriptor(
                            value_id=gene_identifier_updater.find_identifier(
                                variant_entry["GeneSymbol"]
                            ),
                            symbol=variant_entry["GeneSymbol"],
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
        except TypeError:
            print("N/A value", variant_entry)
            return None
        except KeyError:
            print(f"Unable to find gene_symbol for {variant_entry['entrez_id']}")
            return None

    @staticmethod
    def add_gene_genomic_interpretation(
        gene_to_phenotype_entry: dict, gene_identifier_updater: GeneIdentifierUpdater
    ):
        try:
            gene_symbol = gene_identifier_updater.obtain_gene_symbol_from_identifier(
                str(gene_to_phenotype_entry["entrez_id"])
            )
            return GenomicInterpretation(
                subject_or_biosample_id="patient1",
                interpretation_status=4
                if gene_to_phenotype_entry["disease_name"].startswith("?") is False
                else 0,
                gene=GeneDescriptor(
                    value_id=gene_identifier_updater.find_identifier(gene_symbol),
                    symbol=gene_symbol,
                ),
            )
        except KeyError:
            print(f"Unable to find gene_symbol for {gene_to_phenotype_entry['entrez_id']}")
            return None
        except TypeError:
            print("N/A value", gene_to_phenotype_entry)
            return None

    def create_variant_genomic_interpretations(
        self, filtered_variant_summary: pl.DataFrame, gene_identifier_updater
    ):
        genomic_interpretations = []
        for variant_entry in filtered_variant_summary.rows(named=True):
            genomic_interpretation = self.add_variant_genomic_interpretation(
                variant_entry, gene_identifier_updater
            )
            if genomic_interpretation is not None:
                genomic_interpretations.append(genomic_interpretation)
        return genomic_interpretations

    def create_gene_genomic_interpretations(
        self, omim_disease_phenotype_gene_map, gene_identifier_updater
    ):
        genomic_interpretations = []
        for phenotype_entry in omim_disease_phenotype_gene_map.rows(named=True):
            genomic_interpretation = self.add_gene_genomic_interpretation(
                phenotype_entry, gene_identifier_updater
            )
            if genomic_interpretation is not None:
                genomic_interpretations.append(genomic_interpretation)
        return genomic_interpretations

    def create_variant_diagnosis(
        self, filtered_variant_summary: pl.DataFrame, disease: Disease, gene_identifier_updater
    ):
        return Diagnosis(
            disease=OntologyClass(
                id=disease.term.id,
                label=disease.term.label,
            ),
            genomic_interpretations=self.create_variant_genomic_interpretations(
                filtered_variant_summary, gene_identifier_updater
            ),
        )

    def create_gene_diagnosis(
        self,
        omim_disease_phenotype_gene_map: pl.DataFrame,
        gene_identifier_updater,
        disease: Disease,
    ):
        genomic_interpretations = self.create_gene_genomic_interpretations(
            omim_disease_phenotype_gene_map, gene_identifier_updater
        )
        return (
            Diagnosis(
                disease=OntologyClass(
                    id=disease.term.id,
                    label=disease.term.label,
                ),
                genomic_interpretations=genomic_interpretations,
            )
            if genomic_interpretations is not None
            else None
        )

    def create_variant_interpretation(self, filtered_variant_summary, gene_identifier_updater):
        phenopacket_util = PhenopacketUtil(self.phenopacket)
        disease = phenopacket_util.return_phenopacket_disease()
        return Interpretation(
            id=disease.term.label + "-interpretation",
            progress_status=0,
            diagnosis=self.create_variant_diagnosis(
                filtered_variant_summary, disease, gene_identifier_updater
            ),
        )

    @staticmethod
    def add_clinvar_resource(phenopacket: Phenopacket):
        phenopacket.meta_data.resources.extend(
            [
                Resource(
                    id="clinvar",
                    name="Clinical Variation",
                    url="https://www.ncbi.nlm.nih.gov/clinvar/",
                    version="2023-04-06",
                    namespace_prefix="clinvar",
                    iri_prefix="https://www.ncbi.nlm.nih.gov/clinvar/variation/",
                )
            ]
        )
        return phenopacket

    def create_gene_interpretation(self, omim_disease_phenotype_gene_map, gene_identifier_updater):
        phenopacket_util = PhenopacketUtil(self.phenopacket)
        disease = phenopacket_util.return_phenopacket_disease()
        diagnosis = self.create_gene_diagnosis(
            omim_disease_phenotype_gene_map, gene_identifier_updater, disease
        )
        return (
            Interpretation(
                id=disease.term.label + "-interpretation",
                progress_status=0,
                diagnosis=diagnosis,
            )
            if diagnosis is not None
            else None
        )

    def add_variant_interpretation_to_phenopacket(
        self, filtered_variant_summary, gene_identifier_updater
    ):
        phenopacket_copy = copy(self.phenopacket)
        interpretation = self.create_variant_interpretation(
            filtered_variant_summary, gene_identifier_updater
        )
        if interpretation is not None:
            phenopacket_copy.interpretations.extend(
                [
                    self.create_variant_interpretation(
                        filtered_variant_summary, gene_identifier_updater
                    )
                ]
            )
            return self.add_clinvar_resource(phenopacket_copy)
        else:
            pass

    def add_gene_interpretation_to_phenopacket(
        self, omim_disease_phenotype_gene_map, gene_identifier_updater
    ):
        phenopacket_copy = copy(self.phenopacket)
        interpretations = [
            self.create_gene_interpretation(
                omim_disease_phenotype_gene_map, gene_identifier_updater
            )
        ]
        if interpretations is not None:
            phenopacket_copy.interpretations.extend(interpretations)
            return phenopacket_copy
