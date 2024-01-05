from pathlib import Path

import polars as pl
from oaklib.implementations import ProntoImplementation
from ontobio import Ontology

from phenotype2phenopacket.utils.phenopacket_utils import (
    PhenotypeAnnotationToPhenopacketConverter,
    SyntheticPatientGenerator,
    write_phenopacket,
)
from phenotype2phenopacket.utils.utils import (
    filter_diseases,
    load_ontology,
    load_ontology_factory,
    return_phenotype_annotation_data,
)


def create_synthetic_patient_phenopacket(
    human_phenotype_ontology: ProntoImplementation,
    omim_disease: pl.DataFrame,
    ontology_factory: Ontology,
    output_dir: Path,
    hpoa_version: str,
):
    """
    Create a synthetic patient phenopacket from a set of phenotype entries for a specific OMIM disease.

    Args:
        human_phenotype_ontology: An instance of ProntoImplementation containing the loaded HPO.
        omim_disease (pl.DataFrame): DataFrame containing phenotype entries for a specific OMIM disease.
        ontology_factory: Created HPO ontology from OntologyFactory
        output_dir (Path): The directory path to write the generated phenopacket.
        hpoa_version (str): The version of the Human Phenotype Ontology Annotation.

    """
    synthetic_patient_generator = SyntheticPatientGenerator(
        omim_disease, human_phenotype_ontology, ontology_factory
    )
    patient_terms = synthetic_patient_generator.patient_term_annotation_set()
    phenopacket_file = PhenotypeAnnotationToPhenopacketConverter(
        human_phenotype_ontology
    ).create_phenopacket(patient_terms, hpoa_version, synthetic_patient_generator.get_onset_range())
    write_phenopacket(
        phenopacket_file.phenopacket, output_dir.joinpath(phenopacket_file.phenopacket_path)
    )


def create_synthetic_patients(
    phenotype_annotation: Path, num_disease: int, omim_id: str, omim_id_list: Path, output_dir: Path
):
    """
    Create a set of synthetic patient phenopackets from a phenotype annotation file.

    Args:
        phenotype_annotation (Path): Path to the phenotype annotation file.
        num_disease (int): Number of diseases to generate synthetic patient phenopackets for.
                           If set to 0, processes all available diseases.
        omim_id (str) : OMIM ID to generate synthetic patient phenopackets for.
        omim_id_list (Path) : Path to the list of OMIM IDs to generate synthetic patient phenopackets for.
        output_dir (Path): Directory path to write the generated phenopackets.

    """
    phenotype_annotation_data = return_phenotype_annotation_data(phenotype_annotation)
    human_phenotype_ontology = load_ontology()
    ontology_factory = load_ontology_factory()
    grouped_omim_diseases = filter_diseases(
        num_disease, omim_id, omim_id_list, phenotype_annotation_data
    )
    for omim_disease in grouped_omim_diseases:
        create_synthetic_patient_phenopacket(
            human_phenotype_ontology,
            omim_disease,
            ontology_factory,
            output_dir,
            phenotype_annotation_data.version,
        )
