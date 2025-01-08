from pathlib import Path

import polars as pl
from oaklib.implementations import ProntoImplementation

from phenotype2phenopacket.utils.phenopacket_utils import (
    PhenotypeAnnotationToPhenopacketConverter,
    SyntheticPatientGenerator,
    write_phenopacket,
)
from phenotype2phenopacket.utils.utils import (
    filter_diseases,
    load_ontology,
    read_omim_id_list,
    return_phenotype_annotation_data,
)


def create_synthetic_patient_phenopacket(
    human_phenotype_ontology: ProntoImplementation,
    omim_disease: pl.DataFrame,
    output_dir: Path,
    pt_id: str,
    hpoa_version: str,
):
    """
    Create a synthetic patient phenopacket from a set of phenotype entries for a specific OMIM disease.

    Args:
        human_phenotype_ontology: An instance of ProntoImplementation containing the loaded HPO.
        omim_disease (pl.DataFrame): DataFrame containing phenotype entries for a specific OMIM disease.
        output_dir (Path): The directory path to write the generated phenopacket.
        pt_id (str): The patient ID.
        hpoa_version (str): The version of the Human Phenotype Ontology Annotation.

    """
    synthetic_patient_generator = SyntheticPatientGenerator(omim_disease, human_phenotype_ontology)
    patient_terms = synthetic_patient_generator.patient_term_annotation_set()
    phenopacket_file = PhenotypeAnnotationToPhenopacketConverter(
        human_phenotype_ontology
    ).create_phenopacket(
        omim_disease_df=patient_terms,
        hpoa_version=hpoa_version,
        pt_id=pt_id,
        onset=synthetic_patient_generator.get_onset_range(),
    )
    write_phenopacket(
        phenopacket_file.phenopacket, output_dir.joinpath(phenopacket_file.phenopacket_path)
    )


def create_synthetic_patients(
    phenotype_annotation: Path,
    num_disease: int,
    omim_id: str,
    omim_id_list: Path,
    output_dir: Path,
    local_cached_ontology: Path,
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
        local_cached_ontology (Path): Path to the local cached ontology.

    """
    phenotype_annotation_data = return_phenotype_annotation_data(phenotype_annotation)
    human_phenotype_ontology = load_ontology(local_cached_ontology)
    grouped_omim_diseases = filter_diseases(
        num_disease, omim_id, omim_id_list, phenotype_annotation_data
    )
    omim_ids = (
        read_omim_id_list(omim_id_list) if omim_id_list else [None] * len(grouped_omim_diseases)
    )
    for omim_id, omim_disease in zip(omim_ids, grouped_omim_diseases):
        if len(omim_disease) == 0:
            id_message = f" for {omim_id}!" if omim_id else "!"
            print(f"Skipping... Could not find any phenotype entries{id_message}")
            continue
        create_synthetic_patient_phenopacket(
            human_phenotype_ontology,
            omim_disease,
            output_dir,
            None,
            phenotype_annotation_data.version,
        )
