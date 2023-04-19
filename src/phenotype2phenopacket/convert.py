import re
from pathlib import Path

import click
import polars as pl
from google.protobuf.json_format import MessageToJson
from google.protobuf.timestamp_pb2 import Timestamp
from oaklib.implementations.pronto.pronto_implementation import ProntoImplementation
from oaklib.resource import OntologyResource
from phenopackets import (
    Disease,
    Individual,
    MetaData,
    OntologyClass,
    Phenopacket,
    PhenotypicFeature,
    Resource,
    TimeElement,
)


def load_ontology():
    """Loads human phenotype ontology."""
    resource = OntologyResource(slug="hp.obo", local=False)
    return ProntoImplementation(resource)


def create_json_message(phenopacket: Phenopacket) -> str:
    """Creates json message for writing to file."""
    return MessageToJson(phenopacket)


def write_phenopacket(phenopacket: Phenopacket, output_file: Path) -> None:
    """Writes a phenopacket."""
    phenopacket_json = create_json_message(phenopacket)
    with open(output_file, "w") as outfile:
        outfile.write(phenopacket_json)
    outfile.close()


def create_phenopacket_file_name_from_disease(disease_name: str) -> Path:
    normalised_string = re.sub(r"\W+", "_", disease_name)
    return Path(normalised_string.replace(" ", "_") + ".json")


def read_phenotype_annotation_file(phenotype_annotation_file_path: Path) -> pl.DataFrame:
    """Read a phenotype annotation file and return a dataframe"""
    return pl.read_csv(phenotype_annotation_file_path, separator="\t", comment_char="#")


def create_individual() -> Individual:
    return Individual(id="patient1")


def create_onset(human_phenotype_ontology, phenotype_annotation_entry: dict) -> TimeElement:
    if phenotype_annotation_entry["onset"] is not None:
        rels = human_phenotype_ontology.entity_alias_map(phenotype_annotation_entry["onset"])
        term = "".join(rels[(list(rels.keys())[0])])
        return TimeElement(
            ontology_class=OntologyClass(id=phenotype_annotation_entry["onset"], label=term)
        )
    else:
        return None


def create_modifier(human_phenotype_ontology, phenotype_annotation_entry: dict) -> [OntologyClass]:
    if phenotype_annotation_entry["modifier"] is not None:
        rels = human_phenotype_ontology.entity_alias_map(phenotype_annotation_entry["modifier"])
        term = "".join(rels[(list(rels.keys())[0])])
        return [OntologyClass(id=phenotype_annotation_entry["modifier"], label=term)]
    else:
        return None


def create_phenotypic_feature(
    human_phenotype_ontology, phenotype_annotation_entry: dict
) -> PhenotypicFeature:
    if phenotype_annotation_entry["aspect"] == "P":
        rels = human_phenotype_ontology.entity_alias_map(phenotype_annotation_entry["hpo_id"])
        hpo_term = "".join(rels[(list(rels.keys())[0])])
        return PhenotypicFeature(
            type=OntologyClass(id=phenotype_annotation_entry["hpo_id"], label=hpo_term),
            onset=create_onset(human_phenotype_ontology, phenotype_annotation_entry),
            modifiers=create_modifier(human_phenotype_ontology, phenotype_annotation_entry),
        )
    else:
        return None


def create_disease(phenotype_annotation_entry: dict) -> Disease:
    return Disease(
        term=OntologyClass(
            id=phenotype_annotation_entry["database_id"],
            label=phenotype_annotation_entry["disease_name"],
        )
    )


def create_phenotype_ontology_resource() -> Resource:
    return Resource(
        id="hp",
        name="human phenotype ontology",
        url="http://purl.obolibrary.org/obo/hp.owl",
        version="hp/releases/2023-04-05",
        namespace_prefix="HP",
        iri_prefix="http://purl.obolibrary.org/obo/HP_",
    )


def create_metadata() -> MetaData:
    timestamp = Timestamp()
    timestamp.GetCurrentTime()
    return MetaData(
        created=timestamp,
        created_by="phenotype2phenopacket",
        resources=[create_phenotype_ontology_resource()],
        phenopacket_schema_version="2.0",
    )


def convert_to_phenopackets(phenotype_annotation: pl.DataFrame, output_dir: Path):
    human_phenotype_ontology = load_ontology()
    omim_diseases = phenotype_annotation.filter(pl.col("database_id").str.starts_with("OMIM"))
    grouped_omim_diseases = omim_diseases.partition_by(by="database_id", maintain_order=True)
    for omim_disease in grouped_omim_diseases:
        phenotypic_features = []
        phenotype_entry = ""
        for phenotype in omim_disease.rows(named=True):
            phenotypic_feature = create_phenotypic_feature(
                human_phenotype_ontology=human_phenotype_ontology,
                phenotype_annotation_entry=phenotype,
            )
            if phenotypic_feature is not None:
                phenotypic_features.append(phenotypic_feature)
            else:
                pass
            phenotype_entry = phenotype
        phenopacket = Phenopacket(
            id=phenotype_entry["disease_name"].lower().replace(" ", "_"),
            subject=create_individual(),
            phenotypic_features=phenotypic_features,
            diseases=[create_disease(phenotype_entry)],
            meta_data=create_metadata(),
        )
        write_phenopacket(
            phenopacket,
            output_dir.joinpath(
                create_phenopacket_file_name_from_disease(phenotype_entry["disease_name"])
            ),
        )


@click.command("convert")
@click.option(
    "--phenotype-annotation",
    "-p",
    required=True,
    help="Path to phenotype.hpoa.",
    type=Path,
)
@click.option(
    "--output-dir",
    "-o",
    required=True,
    help="Path to output directory.",
    type=Path,
)
def convert_to_phenopackets_command(phenotype_annotation: Path, output_dir: Path):
    phenotype_annotation_df = read_phenotype_annotation_file(phenotype_annotation)
    convert_to_phenopackets(phenotype_annotation_df, output_dir)
