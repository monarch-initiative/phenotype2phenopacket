from pathlib import Path

import pandas as pd
import polars as pl
from oaklib import OntologyResource
from oaklib.implementations import ProntoImplementation
from ontobio import OntologyFactory


def is_float(element: any) -> bool:
    """Checks whether an element is a float or not."""
    if element is None:
        return False
    try:
        float(element)
        return True
    except ValueError:
        return False


def read_disease_pg(disease_pg: Path) -> pl.DataFrame:
    """Read a disease.pg file, returning a polars dataframe."""
    disease = pl.read_csv(
        disease_pg,
        separator="|",
        new_columns=[
            "database_id",
            "gene_mim_number",
            "disease_name",
            "entrez_id",
            "diagnosis_status",
            "inheritance",
        ],
        has_header=False,
    )
    return disease.filter(pl.col("database_id").str.starts_with("OMIM"))


def load_ontology():
    """Load human phenotype ontology."""
    resource = OntologyResource(slug="hp.obo", local=False)
    return ProntoImplementation(resource)


def load_ontology_factory():
    """Load human phenotype ontology factory."""
    ontology_factory = OntologyFactory()
    return ontology_factory.create("hp")


def read_hgnc_data(hgnc_data_file: Path) -> pd.DataFrame:
    """Read hgnc complete set text file, returning a polars dataframe."""
    return pd.read_csv(
        hgnc_data_file,
        delimiter="\t",
        dtype=str,
    )


def read_phenotype_annotation_file(phenotype_annotation_file_path: Path) -> pl.DataFrame:
    """Read a phenotype annotation file and return a polars dataframe"""
    return pl.read_csv(phenotype_annotation_file_path, separator="\t", comment_char="#")


def all_files(directory: Path) -> list[Path]:
    """Obtains all files from a given directory."""
    files = [path for path in directory.iterdir()]
    files.sort()
    return files
