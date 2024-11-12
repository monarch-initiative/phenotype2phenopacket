import random
from dataclasses import dataclass
from pathlib import Path
from typing import List

import pandas as pd
import polars as pl
from oaklib import OntologyResource, get_adapter
from oaklib.implementations import ProntoImplementation


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


def read_genes_to_disease(genes_to_disease: Path) -> pl.DataFrame:
    """
    Read the genes_to_disease.txt file and return a Polars DataFrame.
    Args:
        genes_to_disease (Path): Path to the genes_to_disease.txt file.
    Returns:
        pl.DataFrame: A  Polars DataFrame containing the contents of the genes_to_disease.txt.
    """
    return pl.read_csv(genes_to_disease, separator="\t")


def load_ontology(local_cached_ontology: Path = None):
    """
    Load the Human Phenotype Ontology (HPO).
    Args:
        local_cached_ontology(Path): Path to the local cached ontology.
    Returns:
        ProntoImplementation: An instance of ProntoImplementation containing the loaded HPO.
    """
    if local_cached_ontology is None:
        return get_adapter("sqlite:obo:hp")
    else:
        resource = OntologyResource(slug=str(local_cached_ontology), local=True)
        return ProntoImplementation(resource)


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


def filter_phenotype_annotation(phenotype_annotation_file_path: Path) -> pl.DataFrame:
    """
    Filter the phenotype annotation, retaining data for OMIM diseases describing phenotype information.
    Args:
        phenotype_annotation_file_path (Path): The path to the phenotype annotation file.

    Returns:
        pl.DataFrame: A Polars DataFrame containing the contents of the filtered phenotype annotation file.
    """
    phenotype_annotation_df = read_phenotype_annotation_file(phenotype_annotation_file_path)
    return phenotype_annotation_df.filter(pl.col("database_id").str.starts_with("OMIM")).filter(
        pl.col("aspect") == "P"
    )


def group_phenotype_annotation(phenotype_annotation_df: pl.DataFrame) -> List[pl.DataFrame]:
    """
    Group the phenotype annotation by database ID.

    Args:
        phenotype_annotation_df (pl.Dataframe): The filtered phenotype annotation dataframe.

    Returns:
        List[pl.DataFrame]: A list of Polars DataFrame grouped according to the database ID.
    """
    return phenotype_annotation_df.partition_by(by="database_id", maintain_order=True)


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


@dataclass
class PhenotypeAnnotation:
    """
    Class to represent a phenotype annotation data.

    Attributes:
        df (pl.DataFrame): The phenotype annotation represented as a dataframe.
        version (str): The version of the phenotype annotation.
    """

    df: pl.DataFrame
    version: str


def return_phenotype_annotation_data(phenotype_annotation_file_path: Path) -> PhenotypeAnnotation:
    """
    Read a phenotype annotation file and return the phenotype annotation data as a Polars dataframe
    and the version of the phenotype annotation as a PhenotypeAnnotation object.

    Args:
        phenotype_annotation_file_path (Path): The path to the phenotype annotation file.

    Returns:
        PhenotypeAnnotation: The phenotype annotation data object containing the phenotype
        data for all OMIM diseases and version of the file.
    """
    return PhenotypeAnnotation(
        df=filter_phenotype_annotation(phenotype_annotation_file_path),
        version=get_phenotype_annotation_version(phenotype_annotation_file_path),
    )


def read_omim_id_list(omim_id_list_file_path: Path) -> List[str]:
    """
    Read a txt file that contains OMIM IDs separated by a new line.

    Args:
        omim_id_list_file_path (Path): path to the txt file containing OMIM IDs.

    Returns:
        List[str]: A list of OMIM ids contained in the file.
    """
    with open(omim_id_list_file_path) as f:
        lines = [line.rstrip("\n") for line in f]
    f.close()
    return lines


def filter_diseases(
    num_disease: int,
    omim_id: str,
    omim_id_list: Path,
    phenotype_annotation_data: PhenotypeAnnotation,
) -> List[pl.DataFrame]:
    """
    Filter the phenotype annotation data to either only a specific disease, a specific number of diseases,
    or a list of specific diseases.

    Args:
        num_disease (int): The number of diseases to be filtered. If set to 0, returns all available diseases.
        omim_id (str): The specific OMIM ID used to filter diseases. If provided, filters based on this ID.
        omim_id_list (Path): Path to the file containing OMIM IDs to filter for.
        phenotype_annotation_data (PhenotypeAnnotation): The PhenotypeAnnotation data to be filtered.

    Returns:
        List[pd.DataFrame]: A list of pandas DataFrames representing the filtered diseases.

    """
    if omim_id_list is not None:
        omim_ids = read_omim_id_list(omim_id_list)
        selected_diseases = []
        for omim_id in omim_ids:
            selected_diseases.append(
                phenotype_annotation_data.df.filter(pl.col("database_id") == omim_id)
            )
        return selected_diseases
    if omim_id is not None:
        grouped_omim_disease = phenotype_annotation_data.df.filter(pl.col("database_id") == omim_id)
        return (
            [grouped_omim_disease.clone() for _ in range(num_disease)]
            if num_disease != 0
            else [grouped_omim_disease]
        )
    if num_disease != 0 and omim_id is None:
        return random.choices(
            group_phenotype_annotation(phenotype_annotation_data.df), k=num_disease
        )
    if num_disease == 0 and omim_id is None and omim_id_list is None:
        return group_phenotype_annotation(phenotype_annotation_data.df)
