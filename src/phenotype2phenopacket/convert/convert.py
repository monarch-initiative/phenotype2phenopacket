from pathlib import Path

from phenotype2phenopacket.utils.phenopacket_utils import (
    PhenotypeAnnotationToPhenopacketConverter,
    write_phenopacket,
)
from phenotype2phenopacket.utils.utils import (
    filter_diseases,
    load_ontology,
    read_omim_id_list,
    return_phenotype_annotation_data,
)


def convert_to_phenopackets(
    phenotype_annotation: Path,
    num_disease: int,
    omim_id: str,
    omim_id_list: Path,
    output_dir: Path,
    local_cached_ontology: Path,
):
    """
    Convert a phenotype annotation file to a set of disease-specific phenopackets.

    Args:
        phenotype_annotation (Path): Path to the phenotype annotation file.
        num_disease (int): Number of diseases to convert to phenopackets.
                           If set to 0, processes all available diseases.
        omim_id (str) : OMIM ID to generate synthetic patient phenopackets for.
        omim_id_list (Path) : Path to the list of OMIM IDs to generate synthetic patient phenopackets for.
        output_dir (Path): Directory path to write the generated phenopackets.
        local_cached_ontology (Path): Path to the local cached ontology.
    """
    human_phenotype_ontology = load_ontology(local_cached_ontology)
    phenotype_annotation_data = return_phenotype_annotation_data(phenotype_annotation)
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
        phenopacket_file = PhenotypeAnnotationToPhenopacketConverter(
            human_phenotype_ontology
        ).create_phenopacket(omim_disease, phenotype_annotation_data.version, None)
        write_phenopacket(
            phenopacket_file.phenopacket, output_dir.joinpath(phenopacket_file.phenopacket_path)
        )
