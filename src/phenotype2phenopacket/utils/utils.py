from pathlib import Path

import pandas as pd
import polars as pl
from oaklib import OntologyResource
from oaklib.implementations import ProntoImplementation
from ontobio import Ontology, OntologyFactory


def is_float(element: any) -> bool:
    """
    Checks whether an element is a float or not.

    Args:
        element (any): The element to be checked.

    Returns:
        bool: True if the element is a float, False otherwise.
    """
    if element is None:
        return False
    try:
        float(element)
        return True
    except ValueError:
        return False


def read_disease_pg(disease_pg: Path) -> pl.DataFrame:
    """
    Read a disease.pg file and return a filtered Polars DataFrame.

    This function reads the contents of a 'disease.pg' file using Polars read_csv method
    and constructs a DataFrame. It filters the DataFrame to include only rows where the 'database_id'
    column starts with 'OMIM'.

    Args:
        disease_pg (Path): The path to the 'disease.pg' file.

    Returns:
        pl.DataFrame: A filtered Polars DataFrame containing specific columns and rows
                      where 'database_id' starts with 'OMIM'.
    """
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
    """
    Load the Human Phenotype Ontology (HPO).

    Returns:
        ProntoImplementation: An instance of ProntoImplementation containing the loaded HPO.
    """
    resource = OntologyResource(slug="hp.obo", local=False)
    return ProntoImplementation(resource)


def load_ontology_factory() -> Ontology:
    """
    Load human phenotype ontology factory.

    Returns:
        Ontology: An instance of the hp Ontology class.
    """
    ontology_factory = OntologyFactory()
    return ontology_factory.create("hp")


def read_hgnc_data(hgnc_data_file: Path) -> pd.DataFrame:
    """
    Read HGNC data from a file and return it as a Pandas DataFrame.

    Returns:
        pd.DataFrame: DataFrame containing the HGNC data.
    """
    return pd.read_csv(
        hgnc_data_file,
        delimiter="\t",
        dtype=str,
    )


def read_phenotype_annotation_file(phenotype_annotation_file_path: Path) -> pl.DataFrame:
    """
    Read a phenotype annotation file and return a Polars DataFrame.

    This function reads the contents of a phenotype annotation file using Polars read_csv method.
    The file is assumed to have tab-separated values and comments prefixed with '#'.

    Args:
        phenotype_annotation_file_path (Path): The path to the phenotype annotation file.

    Returns:
        pl.DataFrame: A Polars DataFrame containing the contents of the phenotype annotation file.
    """
    return pl.read_csv(phenotype_annotation_file_path, separator="\t", comment_char="#")


def get_phenotype_annotation_version(phenotype_annotation_file_path: Path) -> str:
    """
    Read a phenotype annotation file and return the version of the phenotype annotation.

    This function reads the contents of a phenotype annotation file and extracts the version information.
    It searches for a line starting with '#version' and extracts the version value after the colon (if found).

    Args:
        phenotype_annotation_file_path (Path): The path to the phenotype annotation file.

    Returns:
        str: The version of the phenotype annotation extracted from the file.
             Returns None if the version information is not found.
    """
    version = None
    with open(phenotype_annotation_file_path, "r") as f:
        for line in f:
            if line.startswith("#version"):
                version = line.split("#version: ")[1].strip("\n")
    f.close()
    return version


def all_files(directory: Path) -> list[Path]:
    """
    Obtains all files from a given directory.

    Args:
        directory (Path): The directory path.

    Returns:
        list[Path]: A list of Path objects representing all files in the directory.
    """
    files = [path for path in directory.iterdir()]
    files.sort()
    return files
