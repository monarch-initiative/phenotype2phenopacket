import re
import secrets
import threading
import warnings
from copy import copy
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import List, Union

import polars as pl
from google.protobuf.timestamp_pb2 import Timestamp
from oaklib.implementations import ProntoImplementation
from phenopackets import (
    Age,
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
from pheval.utils.phenopacket_utils import GeneIdentifierUpdater, create_json_message

from phenotype2phenopacket.utils.utils import is_float


@dataclass
class FrequencyTerm:
    """
    Represents a frequency range with upper and lower values.

    Attributes:
        lower (Union[int, str]): The lower bound of the frequency range.
            Can be an integer or a string representing the lower frequency value.
        upper (Union[int, str]): The upper bound of the frequency range.
            Can be an integer or a string representing the upper frequency value.
    """

    lower: Union[int, str]
    upper: Union[int, str]


frequency_hpo = {
    "HP:0040281": FrequencyTerm(0, 0),
    "HP:0040282": FrequencyTerm(1, 4),
    "HP:0040283": FrequencyTerm(5, 29),
    "HP:0040284": FrequencyTerm(30, 79),
    "HP:0040285": FrequencyTerm(80, 99),
    "HP:0040286": FrequencyTerm(100, 100),
    "HP:0040280": FrequencyTerm(100, 100),
}


@dataclass
class OnsetTerm:
    """
    Represents a range of onset ages with the lowest and highest values.

    Attributes:
        lower_age (Union[int, str]): The lowest age of onset.
            Can be an integer or a string representing the lowest age value.
        upper_age (Union[int, str]): The highest age of onset.
            Can be an integer or a string representing the highest age value.
    """

    lower_age: Union[int, str]
    upper_age: Union[int, str]


onset_hpo = {
    "HP:0011462": OnsetTerm(16, 40),
    "HP:0011460": OnsetTerm(0, 0),
    "HP:0011463": OnsetTerm(1, 5),
    "HP:0003584": OnsetTerm(60, 90),
    "HP:0025709": OnsetTerm(19, 25),
    "HP:0034198": OnsetTerm(0, 0),
    "HP:0003596": OnsetTerm(40, 60),
    "HP:0003621": OnsetTerm(5, 15),
    "HP:0003593": OnsetTerm(0, 1),
    "HP:4000040": OnsetTerm(0, 0),
    "HP:0011461": OnsetTerm(0, 0),
    "HP:0003577": OnsetTerm(0, 0),
    "HP:0034197": OnsetTerm(0, 0),
    "HP:0025708": OnsetTerm(16, 19),
    "HP:0410280": OnsetTerm(1, 15),
    "HP:0030674": OnsetTerm(0, 0),
    "HP:0003623": OnsetTerm(0, 0),
    "HP:0034199": OnsetTerm(0, 0),
    "HP:0003581": OnsetTerm(16, 80),
    "HP:0025710": OnsetTerm(25, 40),
}


@dataclass
class PhenopacketFile:
    """
    Represents a proband presented as a Phenopacket object and the corresponding phenopacket file name.

    Attributes:
        phenopacket: Phenopacket object of the proband.
        phenopacket_path: Path to the phenopacket file.
    """

    phenopacket: Phenopacket
    phenopacket_path: Path


def create_phenopacket_file_name_from_disease(disease_name: str) -> Path:
    """
    Create a Phenopacket file name from the disease.

    Args:
        disease_name (str): The name of the disease.
    """
    normalised_string = re.sub(r"\W+", "_", disease_name)
    return Path(normalised_string.replace(" ", "_") + ".json")


def write_phenopacket(phenopacket: Phenopacket, output_file: Path) -> None:
    """
    Write a Phenopacket object to a file in JSON format.

    Args:
        phenopacket (Phenopacket): The Phenopacket object to be written.
        output_file (Path): The Path object representing the file to write the Phenopacket data.
    """
    phenopacket_json = create_json_message(phenopacket)
    suffix = 1
    if "_patient_" not in output_file.stem:
        while Path(
            output_file.parents[0].joinpath(f"{output_file.stem}_patient_{suffix}.json")
        ).is_file():
            suffix += 1
        output_file = output_file.parents[0].joinpath(f"{output_file.stem}_patient_{suffix}.json")
    else:
        pass
    with open(output_file, "w") as file:
        file.write(phenopacket_json)
    file.close()


class SyntheticPatientGenerator:
    """Class for generating synthetic patients."""

    def __init__(self, disease_df: pl.DataFrame, ontology: ProntoImplementation):
        """
        Initialise the SyntheticPatientGenerator class

        Args:
            disease_df (pl.DataFrame): The dataframe containing the annotation data for a specific disease.
            ontology (ProntoImplementation): An instance of ProntoImplementation containing the loaded HPO.

        """
        self.disease_df = disease_df
        self.ontology = ontology
        self.lower_age = 0
        self.upper_age = 0
        self.filtered_df = []
        self.secret_rand = secrets.SystemRandom()

    def get_number_of_terms(self) -> int:
        """
        Get the number of terms to ascribe from the full set.

        Returns:
            int: Number of terms to ascribe from the full set.
        """
        if len(self.disease_df) == 1:
            return 1
        number_of_terms = 0
        while number_of_terms == 0:
            number_of_terms = int(len(self.disease_df) * (self.secret_rand.uniform(0.2, 0.75)))
        return number_of_terms

    @staticmethod
    def shuffle_dataframe(disease_df: pl.DataFrame) -> pl.DataFrame:
        """
        Shuffle the rows of a dataframe.

        Args:
            disease_df (pl.DataFrame): The dataframe containing the annotation data for a specific disease

        Returns:
            pl.DataFrame: The shuffled dataframe
        """
        return disease_df.sample(fraction=1, shuffle=True)

    def add_frequency(self) -> pl.DataFrame:
        """
        Add random frequency to annotations without one defined.

        Returns:
            pl.DataFrame: DataFrame with random frequency added to annotations without one defined.
        """
        return self.disease_df.with_columns(
            [
                pl.when(self.disease_df["frequency"].is_null())
                .then(pl.lit(self.secret_rand.uniform(0, 1)))
                .otherwise(self.disease_df["frequency"])
                .alias("frequency")
            ]
        )

    def get_onset_range(self) -> OnsetTerm:
        """
        Get the onset range from a set of annotations for a disease.

        Returns:
            OnsetTerm: An instance of OnsetTerm representing the onset range from the annotations.
        """
        for phenotype_entry in self.disease_df.rows(named=True):
            if phenotype_entry["onset"] is not None:
                if self.lower_age < int(onset_hpo[phenotype_entry["onset"]].lower_age):
                    self.lower_age = int(onset_hpo[phenotype_entry["onset"]].lower_age)
                if self.upper_age < int(onset_hpo[phenotype_entry["onset"]].upper_age):
                    self.upper_age = int(onset_hpo[phenotype_entry["onset"]].upper_age)
        return OnsetTerm(lower_age=self.lower_age, upper_age=self.upper_age)

    def check_hpo_frequency(self, phenotype_entry: dict) -> None:
        """
        Filter annotations based on HPO-defined frequencies.

        Args:
            phenotype_entry (dict): A dictionary representing an HPO phenotype entry.

        Notes:
            This method filters phenotype entries based on their defined frequencies from the HPO database.
            If the frequency is 100-100 or if a random frequency falls within the defined limits,
            the phenotype_entry will be appended to the filtered list.
        """
        frequency_limits = frequency_hpo[phenotype_entry["frequency"]]
        if (
            frequency_limits.lower == 100
            and frequency_limits.upper == 100
            and phenotype_entry not in self.filtered_df
        ):
            self.filtered_df.append(phenotype_entry)
        else:
            random_frequency = self.secret_rand.uniform(0, 100)
            if (
                float(frequency_limits.lower) < random_frequency < float(frequency_limits.upper)
                and phenotype_entry not in self.filtered_df
            ):
                self.filtered_df.append(phenotype_entry)

    def check_frequency_threshold(
        self, frequency: float, phenotype_entry: dict, random_frequency: float
    ):
        """
        Check if patient frequency meets the filter for the disease frequency.

        Args:
            frequency (float): The disease frequency threshold.
            phenotype_entry (dict): A dictionary representing an HPO phenotype entry.
            random_frequency (float): The random frequency of the patient.

        Notes:
            This method compares the patient's random frequency against the disease frequency threshold.
            If the patient's random frequency is less than or equal to the disease frequency and
            the phenotype_entry is not already in the filtered list, it will be appended.
        """
        if random_frequency <= float(frequency) and phenotype_entry not in self.filtered_df:
            self.filtered_df.append(phenotype_entry)
        else:
            pass

    def check_percentage_frequency(self, phenotype_entry: dict):
        """
        Filter annotations with percentage-based frequency.

        Args:
            phenotype_entry (dict): A dictionary representing an HPO phenotype entry.
        """
        frequency = phenotype_entry["frequency"].strip("%")
        random_frequency = self.secret_rand.uniform(0, 100)
        self.check_frequency_threshold(frequency, phenotype_entry, random_frequency)

    def check_fraction_frequency(self, phenotype_entry: dict):
        """
        Filter annotations with fraction-based frequency.

        Args:
            phenotype_entry (dict): A dictionary representing an HPO phenotype entry.
        """
        random_frequency = self.secret_rand.uniform(0, 1)
        frequency = float(Fraction(phenotype_entry["frequency"]))
        self.check_frequency_threshold(frequency, phenotype_entry, random_frequency)

    def check_float_frequency(self, phenotype_entry: dict):
        """
        Filter annotations with float-based frequency.

        Args:
            phenotype_entry (dict): A dictionary representing an HPO phenotype entry.
        """
        random_frequency = self.secret_rand.uniform(0, 1)
        frequency = float((phenotype_entry["frequency"]))
        self.check_frequency_threshold(frequency, phenotype_entry, random_frequency)

    def check_frequency(self, phenotype_entry: dict):
        """
        Filter phenotype entries based on different types of frequency representations.

        Args:
            phenotype_entry (dict): A dictionary representing an HPO phenotype entry.

        Notes:
            This method categorises and filters phenotype entries based on different types of frequency representations:
            HPO-defined frequencies starting with "HP:",
            percentage-based frequencies ending with "%",
            frequencies represented as fractions (containing "/"),
            and frequencies represented as floating-point numbers.
        """
        if str(phenotype_entry["frequency"]).startswith("HP:"):
            self.check_hpo_frequency(phenotype_entry)
        elif str(phenotype_entry["frequency"]).endswith("%"):
            self.check_percentage_frequency(phenotype_entry)
        elif "/" in str(phenotype_entry["frequency"]):
            self.check_fraction_frequency(phenotype_entry)
        elif is_float(phenotype_entry["frequency"]) is True:
            self.check_float_frequency(phenotype_entry)

    def filter_phenotype_entries(self, frequency_df: pl.DataFrame, max_number: int):
        """
        Filter annotations based on frequency.

        Args:
            frequency_df (pl.DataFrame): DataFrame containing phenotype frequency data.
            max_number (int): Maximum number of phenotype entries to filter.

        Returns:
            pl.DataFrame: Filtered phenotype entries DataFrame.

        Notes:
            This method filters phenotype entries based on frequency constraints until the maximum number is reached.
            It sets a time limit for execution and handles timeouts by returning filtered results or sampled data.
        """
        time_limit = 15

        def worker():
            try:
                while len(self.filtered_df) < max_number:
                    for phenotype_entry in frequency_df.rows(named=True):
                        if len(self.filtered_df) >= max_number:
                            break
                        self.check_frequency(phenotype_entry)
            except Exception as e:
                print("Error in worker thread:", e)

        thread = threading.Thread(target=worker)
        thread.daemon = True
        stop_event = threading.Event()
        thread.start()
        thread.join(timeout=time_limit)

        if thread.is_alive():
            stop_event.set()
            print("Timed out!")
            if len(self.filtered_df) == 0:
                return frequency_df.sample(n=max_number)
            else:
                return pl.from_dicts(self.filtered_df, infer_schema_length=len(self.filtered_df))
        else:
            return pl.from_dicts(self.filtered_df, infer_schema_length=len(self.filtered_df))

    def get_patient_terms(self) -> pl.DataFrame:
        """
        Get patient terms filtered on frequency thresholds.

        Returns:
            pl.DataFrame: DataFrame containing patient terms filtered on frequency thresholds.

        Notes:
            This method retrieves patient terms by applying various filters based on frequency thresholds.
            It shuffles the dataframe, adds frequency information, and filters phenotype entries
            based on frequency constraints to obtain patient terms.
        """
        return self.filter_phenotype_entries(
            self.shuffle_dataframe(self.add_frequency()), self.get_number_of_terms()
        )

    def get_number_of_terms_to_randomise(self, patient_terms: pl.DataFrame) -> int:
        """
        Get the number of terms to randomise from the filtered frequency set.

        Args:
            patient_terms (pl.DataFrame): DataFrame containing the filtered frequency set of patient terms.

        Returns:
            int: Number of terms to be randomly selected from the filtered frequency set.
        """
        return self.secret_rand.randint(0, int(len(patient_terms)))

    def get_number_of_steps_for_randomisation(self) -> int:
        """
        Get the number of steps to take in the range of 1-5 for making a term more or less specific.

        Returns:
            int: Number of steps to adjust a term's specificity during randomization (1 to 5).
        """
        return self.secret_rand.randint(1, 5)

    def return_less_or_more_specific(self) -> float:
        """
        Generate a random float between 0 and 1.

        Returns:
            float: A random float value between 0 and 1.
        """
        return self.secret_rand.random()

    def subsample_patient_terms(self, patient_terms: pl.dataframe) -> pl.dataframe:
        """
        Get a subsample of patient terms to make more or less specific.


        Args:
            patient_terms (pl.DataFrame): DataFrame containing patient terms.

        Returns:
            pl.DataFrame: Subsample of patient terms for making more or less specific.
        """
        return patient_terms.sample(self.get_number_of_terms_to_randomise(patient_terms))

    def get_children_of_term(self, phenotype_entry: dict, steps: int) -> dict:
        """
        Get a child term of an HPO ID from the specified number of steps.

        Args:
            phenotype_entry (dict): A dictionary representing an HPO phenotype entry.
            steps (int): The number of steps to search for child terms.

        Returns:
            dict: A dictionary representing the updated phenotype entry with a new HPO ID.

        Notes:
            This method retrieves a child term of an HPO ID based on the specified number of steps.
            It iterates through the ontology to find descendants of the given HPO ID
            and selects a child term randomly if available within the specified steps.
        """
        term_id = phenotype_entry["hpo_id"]
        for _i in range(steps):
            descendants = list(self.ontology.incoming_relationships(term_id))
            if descendants:
                term_id = self.secret_rand.choice(descendants)[1]
            else:
                break
        phenotype_entry["hpo_id"] = term_id
        return phenotype_entry

    def get_parents_of_terms(self, phenotype_entry: dict, steps: int) -> dict:
        """
        Get a parent term of an HPO ID from the specified number of steps.

        Args:
            phenotype_entry (dict): A dictionary representing an HPO phenotype entry.
            steps (int): The number of steps to search for parent terms.

        Returns:
            dict: A dictionary representing the updated phenotype entry with a new HPO ID.

        Notes:
            This method retrieves a parent term of an HPO ID based on the specified number of steps.
            It iterates through the ontology to find hierarchical parents of the given HPO ID
            and selects a parent term randomly if available within the specified steps.
        """
        term_id = phenotype_entry["hpo_id"]
        rels = self.ontology.entity_alias_map(term_id)
        term = "".join(rels[(list(rels.keys())[0])]) if rels else ""
        if term.startswith("Abnormality of"):
            return phenotype_entry
        for _i in range(steps):
            parents = self.ontology.hierarchical_parents(term_id)
            if not parents:
                warnings.warn(f"No parents found for term {term_id}", stacklevel=2)
                return phenotype_entry
            parent = self.secret_rand.choice(parents)
            rels = self.ontology.entity_alias_map(parent)
            term = "".join(rels[(list(rels.keys())[0])]) if rels else ""
            if (
                term.startswith("Abnormality of")
                or term_id == "HP:0000118"
                or term_id == "HP:0032443"
            ):
                break
            else:
                term_id = parent
        phenotype_entry["hpo_id"] = term_id
        return phenotype_entry

    @staticmethod
    def remove_terms_to_be_randomised(
        patient_terms: pl.DataFrame, subset: pl.DataFrame
    ) -> pl.DataFrame:
        """
        Remove terms selected for randomisation from patient terms.

        Args:
            patient_terms (pl.DataFrame): DataFrame containing patient terms.
            subset (pl.DataFrame): DataFrame representing terms selected for randomisation.

        Returns:
            pl.DataFrame: DataFrame with selected terms removed for randomisation.
        """
        subset_list = subset.to_dicts()
        for phenotype_entry in subset_list:
            patient_terms = patient_terms.filter(pl.col("hpo_id") != phenotype_entry["hpo_id"])
        return patient_terms

    def alter_term_specificity(
        self, new_phenotype_terms: List[dict], phenotype_entry: dict
    ) -> List[dict]:
        """
        Alter the specificity of the HPO ID.

        Args:
            new_phenotype_terms (List[dict]): List containing updated phenotype terms.
            phenotype_entry (dict): A dictionary representing an HPO phenotype entry.

        Returns:
            List[dict]: Updated list of phenotype terms with altered specificity.

        Notes:
            This method adjusts the specificity of the HPO ID by making it less specific
            if a randomly generated float is less than 0.5; otherwise, it makes the term more specific.
        """
        if self.return_less_or_more_specific() < 0.5:
            new_phenotype_terms.append(
                self.get_parents_of_terms(
                    phenotype_entry, self.get_number_of_steps_for_randomisation()
                )
            )
        else:
            new_phenotype_terms.append(
                self.get_children_of_term(
                    phenotype_entry, self.get_number_of_steps_for_randomisation()
                )
            )
        return new_phenotype_terms

    def patient_term_annotation_set(self) -> pl.DataFrame:
        """
        Get the final patient term annotation set.

        This method generates the final set of patient term annotations by performing the following steps:
        1. Retrieves patient terms and returns them if only one term is present.
        2. Sub-samples patient terms to create a smaller subset.
        3. Alters the specificity of each phenotype entry in the subset.
        4. Removes terms selected for randomisation from the full patient terms.
        5. Combines the altered terms with the filtered terms to form the final patient term set.

        Returns:
            pl.DataFrame: DataFrame containing the final patient term annotation set.
        """
        patient_terms = self.get_patient_terms()
        if len(patient_terms) == 1:
            return patient_terms
        patient_terms_sub_sample = self.subsample_patient_terms(patient_terms)
        new_phenotype_terms = []
        for phenotype_entry in patient_terms_sub_sample.rows(named=True):
            new_phenotype_terms = self.alter_term_specificity(new_phenotype_terms, phenotype_entry)
        patient_terms_filtered = self.remove_terms_to_be_randomised(
            patient_terms, patient_terms_sub_sample
        )
        final_patient_terms = patient_terms_filtered.to_dicts() + new_phenotype_terms
        return pl.from_dicts(final_patient_terms, infer_schema_length=len(final_patient_terms))


class PhenotypeAnnotationToPhenopacketConverter:
    """Class for converting a set of phenotype annotations to a phenopacket format."""

    def __init__(self, human_phenotype_ontology: ProntoImplementation):
        """
        Initialise the PhenotypeAnnotationToPhenopacketConverter class.

        Args:
            human_phenotype_ontology (ProntoImplementation): An instance of ProntoImplementation
            containing the loaded HPO.
        """

        self.human_phenotype_ontology = human_phenotype_ontology
        self.secrets_random_num = secrets.SystemRandom()

    def create_individual(
        self, onset_range: OnsetTerm = None, pt_id: str = "patient1"
    ) -> Individual:
        """
        Create an Individual object.

        Args:
            onset_range (OnsetTerm, optional): An OnsetTerm object representing the age range of onset.
            Defaults to None.
            pt_id (str, optional): An identifier for the patient.

        Returns:
            Individual: An instance of the Individual class.
        """
        age = (
            self.secrets_random_num.randint(onset_range.lower_age, onset_range.upper_age)
            if onset_range is not None
            else None
        )
        if onset_range is None or onset_range.upper_age == 0 and onset_range.lower_age == 0:
            age = None
        return Individual(
            id="patient1" if pt_id is None else pt_id,
            time_at_last_encounter=(
                TimeElement(age=Age(iso8601duration=f"P{age}Y")) if age is not None else None
            ),
        )

    def create_onset(self, phenotype_annotation_entry: dict) -> TimeElement:
        """
        Create an Onset object.

        This method creates a TimeElement representing the onset information for a phenotype annotation entry.

        Args:
            phenotype_annotation_entry (dict): A dictionary representing a phenotype annotation entry.

        Returns:
            TimeElement or None: An instance of TimeElement if onset information is available, else None.
        """
        if phenotype_annotation_entry["onset"] is not None:
            rels = self.human_phenotype_ontology.entity_alias_map(
                phenotype_annotation_entry["onset"]
            )
            term = "".join(rels[(list(rels.keys())[0])]) if rels else None
            return TimeElement(
                ontology_class=OntologyClass(id=phenotype_annotation_entry["onset"], label=term)
            )
        else:
            return None

    def create_modifier(self, phenotype_annotation_entry: dict) -> List[OntologyClass]:
        """
        Create a Modifier.

        This method creates an OntologyClass representing the modifier information for a phenotype annotation entry.

        Args:
            phenotype_annotation_entry (dict): A dictionary representing a phenotype annotation entry.

        Returns:
            list[OntologyClass] or None: A list containing an OntologyClass if modifier information exists,
                                          otherwise returns None.
        """
        if phenotype_annotation_entry["modifier"] is not None:
            ontology_class = []
            for modifier in list(set(phenotype_annotation_entry["modifier"].split(";"))):
                try:
                    rels = self.human_phenotype_ontology.entity_alias_map(modifier)
                    term = "".join(rels[(list(rels.keys())[0])])
                    ontology_class.append(OntologyClass(id=modifier, label=term))
                except IndexError:
                    ontology_class.append(OntologyClass(id=modifier))
            return ontology_class
        else:
            return None

    def create_phenotypic_feature(self, phenotype_annotation_entry: dict) -> PhenotypicFeature:
        """
        Create a PhenotypicFeature object.

        This method generates a PhenotypicFeature object based on a phenotype annotation entry.

        Args:
            phenotype_annotation_entry (dict): A dictionary representing a phenotype annotation entry.

        Returns:
            PhenotypicFeature or None: An instance of PhenotypicFeature if the aspect is 'P' (phenotypic),
            otherwise None.
        """
        if phenotype_annotation_entry["aspect"] == "P":
            rels = self.human_phenotype_ontology.entity_alias_map(
                phenotype_annotation_entry["hpo_id"]
            )
            hpo_term = "".join(rels[(list(rels.keys())[0])]) if rels else None
            return PhenotypicFeature(
                type=OntologyClass(id=phenotype_annotation_entry["hpo_id"], label=hpo_term),
                onset=self.create_onset(phenotype_annotation_entry),
                modifiers=self.create_modifier(phenotype_annotation_entry),
                excluded=True if phenotype_annotation_entry["qualifier"] == "NOT" else False,
            )
        else:
            return None

    def create_phenotypic_features(self, omim_disease_df: pl.DataFrame) -> List[PhenotypicFeature]:
        """
        Create a list of phenotypic features.

        This method generates a list of PhenotypicFeature objects based on a DataFrame of OMIM disease entries.

        Args:
            omim_disease_df (pl.DataFrame): DataFrame containing OMIM disease entries.

        Returns:
            List[PhenotypicFeature]: A list of PhenotypicFeature objects derived from the provided DataFrame.
        """
        phenotypic_features = []
        for phenotype_annotation_entry in omim_disease_df.rows(named=True):
            phenotypic_feature = self.create_phenotypic_feature(phenotype_annotation_entry)
            if phenotypic_feature is not None:
                phenotypic_features.append(phenotypic_feature)
        return phenotypic_features

    @staticmethod
    def create_disease(phenotype_annotation_entry: dict) -> Disease:
        """
        Create a Disease object.
        Args:
            phenotype_annotation_entry (dict): A dictionary representing a phenotype annotation entry.

        Returns:
            Disease: An instance of Disease representing the disease information.
        """
        return Disease(
            term=OntologyClass(
                id=phenotype_annotation_entry["database_id"],
                label=phenotype_annotation_entry["disease_name"],
            )
        )

    @staticmethod
    def create_omim_resource() -> Resource:
        """
        Create OMIM resource.

        Returns:
            Resource: An instance of Resource representing the OMIM database resource.
        """
        return Resource(
            id="omim",
            name="Online Mendelian Inheritance in Man",
            url="https://www.omim.org",
            version="hp/releases/2023-04-18",
            namespace_prefix="OMIM",
            iri_prefix="https://omim.org/entry/",
        )

    @staticmethod
    def create_human_phenotype_ontology_resource(hpoa_version: str) -> Resource:
        """
        Create human phenotype ontology resource.

        Returns:
            Resource: An instance of Resource representing the human phenotype ontology resource.
        """
        return Resource(
            id="hp",
            name="human phenotype ontology",
            url="http://purl.obolibrary.org/obo/hp.owl",
            version="hp/releases/" + hpoa_version,
            namespace_prefix="HP",
            iri_prefix="http://purl.obolibrary.org/obo/HP_",
        )

    def create_metadata(self, hpoa_version: str) -> MetaData:
        """
        Create metadata.

        This method generates metadata for a Phenopacket including timestamp, creator information,
        associated resources, and Phenopacket schema version.

        Args:
            hpoa_version (str): Version of the Human Phenotype Ontology Annotation.

        Returns:
            MetaData: Metadata information for the Phenopacket.
        """
        timestamp = Timestamp()
        timestamp.GetCurrentTime()
        return MetaData(
            created=timestamp,
            created_by="phenotype2phenopacket",
            resources=[
                self.create_human_phenotype_ontology_resource(hpoa_version),
                self.create_omim_resource(),
            ],
            phenopacket_schema_version="2.0",
        )

    def create_phenopacket(
        self,
        omim_disease_df: pl.DataFrame,
        hpoa_version: str,
        pt_id: str,
        onset: OnsetTerm = None,
    ) -> PhenopacketFile:
        """
        Create a Phenopacket object.

        This method generates a PhenopacketFile containing a Phenopacket object with information
        from OMIM disease DataFrame, Human Phenotype Ontology Annotation version, and onset term (if available).

        Args:
            omim_disease_df (pl.DataFrame): DataFrame containing phenotype annotation disease entries.
            hpoa_version (str): Version of the Human Phenotype Ontology Annotation.
            pt_id (str): The patient ID. If not given, defaults to "patient1" in create_individual()
            onset (OnsetTerm, optional): An OnsetTerm object representing the age range of onset. Defaults to None.

        Returns:
            PhenopacketFile: A class instance containing the phenopacket file name and
            the generated Phenopacket.
        """
        phenotypic_features = self.create_phenotypic_features(omim_disease_df)
        phenotype_annotation_entry = omim_disease_df.to_dicts()[0]
        return PhenopacketFile(
            phenopacket=Phenopacket(
                id=phenotype_annotation_entry["disease_name"].lower().replace(" ", "_"),
                subject=self.create_individual(onset, pt_id),
                phenotypic_features=phenotypic_features,
                diseases=[self.create_disease(phenotype_annotation_entry)],
                meta_data=self.create_metadata(hpoa_version),
            ),
            phenopacket_path=create_phenopacket_file_name_from_disease(
                phenotype_annotation_entry["database_id"]
            ),
        )


class PhenopacketUtil:
    """Phenopacket utility class."""

    def __init__(self, phenopacket: Phenopacket):
        """
        Initialise the PhenopacketUtil class.

        Args:
            phenopacket(Phenopacket): The phenopacket object.
        """
        self.phenopacket = phenopacket

    def return_phenopacket_disease(self) -> Disease:
        """
        Return the disease object from a phenopacket.

        Returns:
            Disease: The proband disease.
        """
        return self.phenopacket.diseases[0]


class PhenopacketInterpretationExtender:
    """Class for extending the phenopacket interpretations field."""

    def __init__(self, phenopacket: Phenopacket):
        """
        Initialise the PhenopacketInterpretationExtender class.

        Args:
            phenopacket(Phenopacket): The phenopacket object.
        """
        self.phenopacket = phenopacket

    @staticmethod
    def create_gene_genomic_interpretation(
        gene_to_disease_entry: dict, gene_identifier_updater: GeneIdentifierUpdater
    ) -> GenomicInterpretation:
        """
        Create genomic interpretation for a gene-to-phenotype relationship.

        This method generates a GenomicInterpretation object based on a gene-to-phenotype relationship entry.

        Args:
            gene_to_disease_entry (dict): A dictionary representing a gene-to-disease relationship.
            gene_identifier_updater (GeneIdentifierUpdater): An instance of GeneIdentifierUpdater.

        Returns:
            GenomicInterpretation or None: An instance of GenomicInterpretation representing the interpretation
                                           of the gene-to-phenotype relationship or None if unsuccessful.
        """
        try:
            gene_symbol = gene_to_disease_entry["gene_symbol"]
            return GenomicInterpretation(
                subject_or_biosample_id="patient1",
                interpretation_status=4,
                gene=GeneDescriptor(
                    value_id=gene_identifier_updater.find_identifier(gene_symbol),
                    symbol=gene_symbol,
                ),
            )
        except KeyError:
            print(f"Unable to find gene_symbol for {gene_to_disease_entry['entrez_id']}")
            return None
        except TypeError:
            print("N/A value", gene_to_disease_entry)
            return None

    def create_gene_genomic_interpretations(
        self,
        genes_to_disease_map: pl.DataFrame,
        gene_identifier_updater: GeneIdentifierUpdater,
    ) -> List[GenomicInterpretation]:
        """
        Create a list of genomic interpretations for known gene-to-phenotype relationships.

        This method generates a list of GenomicInterpretation objects based on a DataFrame
        containing known gene-to-phenotype relationships.

        Args:
            genes_to_disease_map (pl.DataFrame): DataFrame containing gene-to-disease mappings.
            gene_identifier_updater (GeneIdentifierUpdater): An instance of GeneIdentifierUpdater.

        Returns:
            List[GenomicInterpretation]: A list of GenomicInterpretation objects representing the interpretations
                  of gene-to-phenotype relationships.
        """
        genomic_interpretations = []
        for disease_entry in genes_to_disease_map.rows(named=True):
            genomic_interpretation = self.create_gene_genomic_interpretation(
                disease_entry, gene_identifier_updater
            )
            if genomic_interpretation is not None:
                genomic_interpretations.append(genomic_interpretation)
        return genomic_interpretations

    def create_gene_diagnosis(
        self,
        genes_to_disease_map: pl.DataFrame,
        gene_identifier_updater: GeneIdentifierUpdater,
        disease: Disease,
    ) -> Diagnosis:
        """
        Create a diagnosis object for known gene-to-phenotype relationships.

        This method generates a Diagnosis object based on known gene-to-phenotype relationships
        provided in a DataFrame and a Disease object.

        Args:
            genes_to_disease_map (pl.DataFrame): DataFrame containing gene-to-disease mappings.
            gene_identifier_updater (GeneIdentifierUpdater): An instance of GeneIdentifierUpdater.
            disease (Disease): An instance of Disease representing the disease information.

        Returns:
            Diagnosis or None: A Diagnosis object representing the diagnosis based on gene-to-phenotype relationships,
                               or None if no genomic interpretations were found.
        """
        genomic_interpretations = self.create_gene_genomic_interpretations(
            genes_to_disease_map, gene_identifier_updater
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

    def create_gene_interpretation(
        self,
        genes_to_disease_map: pl.DataFrame,
        gene_identifier_updater: GeneIdentifierUpdater,
    ) -> Interpretation:
        """
        Create an interpretation object for known gene-to-phenotype relationships.

        This method generates an Interpretation object based on known gene-to-phenotype relationships
        provided in a DataFrame and a GeneIdentifierUpdater instance.

        Args:
            genes_to_disease_map (pl.DataFrame): DataFrame containing gene-to-disease mappings.
            gene_identifier_updater (GeneIdentifierUpdater): An instance of GeneIdentifierUpdater.

        Returns:
            Interpretation or None: An Interpretation object representing the interpretation based on gene-to-phenotype
                                    relationships, or None if no diagnosis was created.
        """
        phenopacket_util = PhenopacketUtil(self.phenopacket)
        disease = phenopacket_util.return_phenopacket_disease()
        diagnosis = self.create_gene_diagnosis(
            genes_to_disease_map, gene_identifier_updater, disease
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

    def add_gene_interpretation_to_phenopacket(
        self,
        genes_to_disease_map: pl.DataFrame,
        gene_identifier_updater: GeneIdentifierUpdater,
    ) -> Phenopacket:
        """
        Add interpretations to a Phenopacket.

        This method adds gene-based interpretations to a copy of the Phenopacket.

        Args:
            genes_to_disease_map (pl.DataFrame): DataFrame containing gene-to-disease mappings.
            gene_identifier_updater (GeneIdentifierUpdater): An instance of GeneIdentifierUpdater.

        Returns:
            Phenopacket or None: A copy of the Phenopacket with added gene-based interpretations,
            or None if interpretations were not created.
        """
        phenopacket_copy = copy(self.phenopacket)
        interpretations = [
            self.create_gene_interpretation(genes_to_disease_map, gene_identifier_updater)
        ]
        if interpretations is not None:
            phenopacket_copy.interpretations.extend(interpretations)
            return phenopacket_copy
