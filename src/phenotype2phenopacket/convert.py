import re
from collections import defaultdict
from pathlib import Path

import click
import pandas as pd
import polars as pl
from google.protobuf.json_format import MessageToJson
from google.protobuf.timestamp_pb2 import Timestamp
from oaklib.implementations.pronto.pronto_implementation import ProntoImplementation
from oaklib.resource import OntologyResource
from phenopackets import (
    Diagnosis,
    Disease,
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
)


def read_disease_pg(disease_pg: Path) -> pl.DataFrame:
    disease = pl.read_csv(disease_pg, separator="|")
    return disease.filter(pl.col("database_id").str.starts_with("OMIM"))


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


def read_hgnc_data(hgnc_data_file: Path) -> pd.DataFrame:
    return pd.read_csv(
        hgnc_data_file,
        delimiter="\t",
        dtype=str,
    )


def create_hgnc_dict(hgnc_data_file: Path) -> defaultdict:
    """Creates reference for updating gene symbols and identifiers."""
    hgnc_df = read_hgnc_data(hgnc_data_file)
    hgnc_data = defaultdict(dict)
    for _index, row in hgnc_df.iterrows():
        previous_names = []
        hgnc_data[row["symbol"]]["ensembl_id"] = row["ensembl_gene_id"]
        hgnc_data[row["symbol"]]["hgnc_id"] = row["hgnc_id"]
        hgnc_data[row["symbol"]]["entrez_id"] = row["entrez_id"]
        hgnc_data[row["symbol"]]["refseq_accession"] = row["refseq_accession"]
        previous = str(row["prev_symbol"]).split("|")
        for p in previous:
            previous_names.append(p.strip('"'))
        hgnc_data[row["symbol"]]["previous_symbol"] = previous_names

    return hgnc_data


def create_gene_identifier_map(hgnc_data_file: Path) -> dict:
    hgnc_df = read_hgnc_data(hgnc_data_file)
    identifier_map = {}
    for _index, row in hgnc_df.iterrows():
        identifier_map[row["ensembl_gene_id"]] = row["symbol"]
        identifier_map[row["hgnc_id"]] = row["symbol"]
        identifier_map[row["entrez_id"]] = row["symbol"]
        identifier_map[row["refseq_accession"]] = row["symbol"]
    return identifier_map


class GeneIdentifierUpdater:
    def __init__(self, gene_identifier: str, hgnc_data: dict = None, identifier_map: dict = None):
        self.hgnc_data = hgnc_data
        self.gene_identifier = gene_identifier
        self.identifier_map = identifier_map

    def find_identifier(self, gene_symbol: str) -> str:
        """Finds the specified gene identifier for a gene symbol."""
        if gene_symbol in self.hgnc_data.keys():
            return self.hgnc_data[gene_symbol][self.gene_identifier]
        else:
            for _symbol, data in self.hgnc_data.items():
                for prev_symbol in data["previous_symbol"]:
                    if prev_symbol == gene_symbol:
                        return data[self.gene_identifier]

    def obtain_gene_symbol_from_identifier(self, query_gene_identifier: str) -> str:
        """
        Obtain gene symbol from a gene identifier. (e.g.)
        "
        obtain_gene_symbol_from_identifier(query_gene_identifier="HGNC:5")
        "
        """
        return self.identifier_map[query_gene_identifier]


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
        try:
            rels = human_phenotype_ontology.entity_alias_map(phenotype_annotation_entry["modifier"])
            term = "".join(rels[(list(rels.keys())[0])])
            return [OntologyClass(id=phenotype_annotation_entry["modifier"], label=term)]
        except IndexError:
            return [OntologyClass(id=phenotype_annotation_entry["modifier"])]
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


def get_phenotype_to_disease_entries(
    omim_disease_pg: pl.DataFrame, phenotype_annotation_entry: dict
) -> pl.DataFrame:
    return omim_disease_pg.filter(
        pl.col("database_id") == phenotype_annotation_entry["database_id"]
    )


def convert_phenotype_map_list_to_genomic_interpretation(
    omim_disease_phenotype_gene_map: pl.DataFrame, gene_identifier_updator: GeneIdentifierUpdater
) -> [GenomicInterpretation]:
    genomic_interpretations = []
    if len(omim_disease_phenotype_gene_map) == 0:
        return None
    else:
        for row in omim_disease_phenotype_gene_map.rows(named=True):
            gene_symbol = gene_identifier_updator.obtain_gene_symbol_from_identifier(
                str(row["entrez_id"])
            )
            genomic_interpretations.append(
                GenomicInterpretation(
                    subject_or_biosample_id="patient1",
                    interpretation_status=4 if row["disease_name"].startswith("?") is False else 0,
                    gene=GeneDescriptor(
                        value_id=gene_identifier_updator.find_identifier(gene_symbol),
                        symbol=gene_symbol,
                    ),
                )
            )
        return genomic_interpretations


def create_disease(phenotype_annotation_entry: dict) -> Disease:
    return Disease(
        term=OntologyClass(
            id=phenotype_annotation_entry["database_id"],
            label=phenotype_annotation_entry["disease_name"],
        )
    )


def create_phenotype_ontology_resource() -> Resource:
    return Resource(
        id="omim",
        name="Online Mendelian Inheritance in Man",
        url="https://www.omim.org",
        version="hp/releases/2023-04-18",
        namespace_prefix="OMIM",
        iri_prefix="https://omim.org/entry/",
    )


def create_omim_resource() -> Resource:
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
        resources=[create_phenotype_ontology_resource(), create_omim_resource()],
        phenopacket_schema_version="2.0",
    )


# def get_genes_from_omim_api(phenotype_annotation_entry: dict) -> [dict]:
#     url = (
#         f"https://api.omim.org/api/entry?mimNumber="
#         f"{phenotype_annotation_entry['database_id'].replace('OMIM:', '')}"
#         f"&format=json&include=geneMap&apiKey=Ec9viXaESlGZTviYKTPUUg"
#     )
#     response = requests.get(url, timeout=60)
#     json_data = response.json()
#     return (
#         json_data["omim"]["entryList"][0]["entry"]["phenotypeMapList"]
#         if "phenotypeMapList"
#         in json_data.get("omim", {}).get("entryList", [{}])[0].get("entry", {})
#         else None
#     )


# def convert_phenotype_map_list_to_genomic_interpretation(
#         phenotype_map_list: [dict], gene_identifier_updator: GeneIdentifierUpdater
# ) -> [GenomicInterpretation]:
#     genomic_interpretations = []
#     for phenotype_map in phenotype_map_list:
#         if phenotype_map["phenotypeMap"]["phenotype"].startswith("?"):
#             continue
#         elif "approvedGeneSymbols" in phenotype_map["phenotypeMap"]:
#             genomic_interpretations.append(
#                 GenomicInterpretation(
#                     subject_or_biosample_id="patient1",
#                     interpretation_status=4
#                     if phenotype_map["phenotypeMap"]["phenotypeMappingKey"] == 3
#                     else 0,
#                     gene=GeneDescriptor(
#                         value_id=gene_identifier_updator.find_identifier(
#                             phenotype_map["phenotypeMap"]["approvedGeneSymbols"]
#                         ),
#                         symbol=phenotype_map["phenotypeMap"]["approvedGeneSymbols"],
#                     ),
#                 )
#             )
#     return genomic_interpretations


def create_interpretations(
    phenotype_annotation_entry: dict,
    gene_identifier_updator: GeneIdentifierUpdater,
    disease_pg: pl.DataFrame,
) -> Interpretation:
    omim_disease_pg = get_phenotype_to_disease_entries(disease_pg, phenotype_annotation_entry)
    genomic_interpretations = convert_phenotype_map_list_to_genomic_interpretation(
        omim_disease_pg, gene_identifier_updator
    )
    if not genomic_interpretations:
        return None
    else:
        return Interpretation(
            id=phenotype_annotation_entry["disease_name"] + "-interpretation",
            progress_status=0,
            diagnosis=Diagnosis(
                disease=OntologyClass(
                    id=phenotype_annotation_entry["database_id"],
                    label=phenotype_annotation_entry["disease_name"],
                ),
                genomic_interpretations=genomic_interpretations,
            ),
        )


def convert_to_phenopackets(
    phenotype_annotation: pl.DataFrame,
    hgnc_data_file: Path,
    disease_pg_df: pl.DataFrame,
    output_dir: Path,
):
    human_phenotype_ontology = load_ontology()
    omim_diseases = phenotype_annotation.filter(pl.col("database_id").str.starts_with("OMIM"))
    grouped_omim_diseases = omim_diseases.partition_by(by="database_id", maintain_order=True)
    hgnc_dict = create_hgnc_dict(hgnc_data_file)
    identifier_map = create_gene_identifier_map(hgnc_data_file)
    gene_identifier_updator = GeneIdentifierUpdater(
        gene_identifier="ensembl_id", hgnc_data=hgnc_dict, identifier_map=identifier_map
    )
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
        interpretations = create_interpretations(
            phenotype_entry, gene_identifier_updator, disease_pg_df
        )
        phenopacket = Phenopacket(
            id=phenotype_entry["disease_name"].lower().replace(" ", "_"),
            subject=create_individual(),
            phenotypic_features=phenotypic_features,
            interpretations=[interpretations] if interpretations is not None else None,
            diseases=[create_disease(phenotype_entry)],
            meta_data=create_metadata(),
        )
        if phenopacket.interpretations:
            output_dir.joinpath("g2p_phenopackets").mkdir(exist_ok=True)
            write_phenopacket(
                phenopacket,
                output_dir.joinpath(
                    f"g2p_phenopackets/{create_phenopacket_file_name_from_disease(phenotype_entry['disease_name'])}"
                ),
            )
        else:
            output_dir.joinpath("phenopackets").mkdir(exist_ok=True)
            write_phenopacket(
                phenopacket,
                output_dir.joinpath(
                    f"phenopackets/{create_phenopacket_file_name_from_disease(phenotype_entry['disease_name'])}"
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
    "--hgnc-genes",
    "-h",
    required=True,
    help="Path to hgnc data file.",
    type=Path,
)
@click.option(
    "--disease-pg",
    "-d",
    required=True,
    help="Path to disease.pg.",
    type=Path,
)
@click.option(
    "--output-dir",
    "-o",
    required=True,
    help="Path to output directory.",
    type=Path,
)
def convert_to_phenopackets_command(
    phenotype_annotation: Path, hgnc_genes: Path, disease_pg: Path, output_dir: Path
):
    phenotype_annotation_df = read_phenotype_annotation_file(phenotype_annotation)
    disease_pg_df = read_disease_pg(disease_pg)
    convert_to_phenopackets(phenotype_annotation_df, hgnc_genes, disease_pg_df, output_dir)
