import re
from dataclasses import dataclass
from pathlib import Path

import polars as pl
from google.protobuf.json_format import MessageToJson
from google.protobuf.timestamp_pb2 import Timestamp
from phenopackets import (  # Diagnosis,; GeneDescriptor,; GenomicInterpretation,; Interpretation,
    Disease,
    Individual,
    MetaData,
    OntologyClass,
    Phenopacket,
    PhenotypicFeature,
    Resource,
    TimeElement,
)


@dataclass
class PhenopacketFile:
    phenopacket: Phenopacket
    phenopacket_path: Path


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
    def create_individual() -> Individual:
        """Create an Individual object."""
        return Individual(id="patient1")

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

    def create_phenopacket(self, omim_disease_df: pl.DataFrame) -> PhenopacketFile:
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
                subject=self.create_individual(),
                phenotypic_features=phenotypic_features,
                diseases=[self.create_disease(phenotype_annotation_entry)],
                meta_data=self.create_metadata(),
            ),
            phenopacket_path=create_phenopacket_file_name_from_disease(
                phenotype_annotation_entry["disease_name"]
            ),
        )
