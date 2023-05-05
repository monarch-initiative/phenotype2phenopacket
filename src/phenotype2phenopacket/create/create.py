import random
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import polars as pl

from typing import Union

from phenotype2phenopacket.utils.phenopacket_utils import (
    write_phenopacket, PhenotypeAnnotationToPhenopacketConverter,
    onset_hpo, OnsetTerm
)
from phenotype2phenopacket.utils.utils import load_ontology, is_float


@dataclass
class FrequencyTerm:
    lower: Union[int, str]
    upper: Union[int, str]


frequency_hpo = {
    "HP:0040281": FrequencyTerm(0, 0),
    "HP:0040282": FrequencyTerm(1, 4),
    "HP:0040283": FrequencyTerm(5, 29),
    "HP:0040284": FrequencyTerm(30, 79),
    "HP:0040285": FrequencyTerm(80, 99),
    "HP:0040286": FrequencyTerm(100, 100),
}





class SyntheticPatientGenerator:
    def __init__(self, disease_df: pl.DataFrame, ontology):
        self.disease_df = disease_df
        self.ontology = ontology
        self.lower_age = 0
        self.upper_age = 0
        self.filtered_df = []

    def get_number_of_terms(self):
        """Get number of terms to ascribe from full set."""
        if len(self.disease_df) == 1:
            return 1
        return int(len(self.disease_df) * (random.uniform(0.2, 0.75)))

    @staticmethod
    def shuffle_dataframe(disease_df):
        """Shuffle dataframe."""
        return disease_df.sample(fraction=1, shuffle=True)

    def add_frequency(self):
        """Add random frequency to annotations without one defined."""
        updated_df = []
        for row in self.disease_df.rows(named=True):
            if row["frequency"] is None:
                synthetic_frequency = random.uniform(0, 1)
                row["frequency"] = synthetic_frequency
                updated_df.append(row)
            elif row["frequency"] is not None:
                updated_df.append(row)
        return pl.from_dicts(updated_df)

    def get_onset_range(self):
        """Get the onset range from a set of annotations for a disease."""
        for phenotype_entry in self.disease_df.rows(named=True):
            if phenotype_entry["onset"] is not None:
                if self.lower_age < float(onset_hpo[phenotype_entry["onset"]].lower_age):
                    self.lower_age = float(onset_hpo[phenotype_entry["onset"]].lower_age)
                if self.upper_age < float(onset_hpo[phenotype_entry["onset"]].upper_age):
                    self.upper_age = float(onset_hpo[phenotype_entry["onset"]].upper_age)
        return OnsetTerm(lower_age=self.lower_age, upper_age=self.upper_age)

    def check_hpo_frequency(self, phenotype_entry):
        """Filter with HPO defined frequencies."""
        frequency_limits = frequency_hpo[phenotype_entry["frequency"]]
        random_frequency = random.uniform(0, 100)
        if float(frequency_limits.lower) < random_frequency < float(
                frequency_limits.upper) and phenotype_entry not in self.filtered_df:
            self.filtered_df.append(phenotype_entry)

    def check_frequency_threshold(self, frequency, phenotype_entry, random_frequency):
        """Check if patient frequency meets the filter for the disease frequency."""
        if random_frequency <= float(frequency) and phenotype_entry not in self.filtered_df:
            self.filtered_df.append(phenotype_entry)
        else:
            pass

    def check_percentage_frequency(self, phenotype_entry):
        """Filter with percentage frequency."""
        frequency = phenotype_entry["frequency"].strip("%")
        random_frequency = random.uniform(0, 100)
        self.check_frequency_threshold(frequency, phenotype_entry, random_frequency)

    def check_fraction_frequency(self, phenotype_entry):
        """Filter with fraction frequency."""
        random_frequency = random.uniform(0, 1)
        frequency = float(Fraction(phenotype_entry["frequency"]))
        self.check_frequency_threshold(frequency, phenotype_entry, random_frequency)

    def check_float_frequency(self, phenotype_entry):
        random_frequency = random.uniform(0, 1)
        frequency = float((phenotype_entry["frequency"]))
        self.check_frequency_threshold(frequency, phenotype_entry, random_frequency)

    def check_frequency(self, phenotype_entry):
        """Filter for frequency."""
        if str(phenotype_entry["frequency"]).startswith("HP:"):
            self.check_hpo_frequency(phenotype_entry)
        elif str(phenotype_entry["frequency"]).endswith("%"):
            self.check_percentage_frequency(phenotype_entry)
        elif "/" in str(phenotype_entry["frequency"]):
            self.check_fraction_frequency(phenotype_entry)
        elif is_float(phenotype_entry["frequency"]) is True:
            self.check_float_frequency(phenotype_entry)

    def filter_phenotype_entries(self, frequency_df: pl.DataFrame, max_number: int):
        """Filter annotations based on frequency."""
        while len(self.filtered_df) < max_number:
            for phenotype_entry in frequency_df.rows(named=True):
                if len(self.filtered_df) >= max_number:
                    break
                self.check_frequency(phenotype_entry)
        return pl.from_dicts(self.filtered_df)

    def get_patient_terms(self):
        """Get patient terms, filtered on frequency thresholds."""
        return self.filter_phenotype_entries(self.shuffle_dataframe(self.add_frequency()),
                                             self.get_number_of_terms())

    @staticmethod
    def get_number_of_terms_to_randomise(patient_terms: pl.DataFrame):
        """Get number of terms to randomise from filtered frequency set."""
        return random.randint(0, int(len(patient_terms)))

    @staticmethod
    def get_number_of_steps_for_randomisation():
        """Get the number of steps to take in range 1-5 for making a term more/less specific."""
        return random.randint(1, 5)

    @staticmethod
    def return_less_or_more_specific():
        """Generate a float between 0-1."""
        return random.random()

    def subsample_patient_terms(self, patient_terms: pl.dataframe):
        """Get a subsample of patient terms to make more/less specific."""
        return patient_terms.sample(self.get_number_of_terms_to_randomise(patient_terms))

    def get_children_of_term(self, phenotype_entry: dict, steps: int):
        """Get a child term of a hpo id from the number of steps specified."""
        term_id = phenotype_entry["hpo_id"]
        for i in range(steps):
            descendants = self.ontology.descendants(phenotype_entry["hpo_id"])
            descendant = random.choice(list(descendants))
            term_id = descendant
        phenotype_entry["hpo_id"] = term_id
        return phenotype_entry

    def get_parents_of_terms(self, phenotype_entry: dict, steps: int):
        """Get a parent term of a hpo id from the number of steps specified."""
        term_id = phenotype_entry["hpo_id"]
        for i in range(steps):
            parents = self.ontology.hierarchical_parents(term_id)
            parent = random.choice(parents)
            rels = self.ontology.entity_alias_map(
                parent
            )
            term = "".join(rels[(list(rels.keys())[0])])
            if term.startswith("Abnormality of") or term_id == 'HP:0000118':
                break
            else:
                term_id = parent
        phenotype_entry["hpo_id"] = term_id
        return phenotype_entry

    @staticmethod
    def remove_terms_to_be_randomised(patient_terms: pl.DataFrame, subset: pl.DataFrame):
        """Remove terms selected to be randomised from patient terms."""
        subset_list = subset.to_dicts()
        for phenotype_entry in subset_list:
            patient_terms = patient_terms.filter(pl.col("hpo_id") != phenotype_entry["hpo_id"])
        return patient_terms

    def alter_term_specificity(self, new_phenotype_terms, phenotype_entry):
        """Alter the hpo id specificity - making less specific if the float is less than 0.5,
        otherwise the term is made more specific."""
        if self.return_less_or_more_specific() < 0.5:
            new_phenotype_terms.append(
                self.get_parents_of_terms(phenotype_entry, self.get_number_of_steps_for_randomisation()))
        else:
            new_phenotype_terms.append(
                self.get_children_of_term(phenotype_entry, self.get_number_of_steps_for_randomisation()))

    def patient_term_annotation_set(self):
        """Get the final patient term annotation set."""
        patient_terms = self.get_patient_terms()
        if len(patient_terms) == 1:
            return patient_terms
        patient_terms_sub_sample = self.subsample_patient_terms(patient_terms)
        new_phenotype_terms = []
        for phenotype_entry in patient_terms_sub_sample.rows(named=True):
            self.alter_term_specificity(new_phenotype_terms, phenotype_entry)
        patient_terms_filtered = self.remove_terms_to_be_randomised(patient_terms, patient_terms_sub_sample)
        final_patient_terms = patient_terms_filtered.to_dicts() + new_phenotype_terms
        return pl.from_dicts(final_patient_terms)


# 1. add missing random frequency,
# 2. randomly sort dataframe
# 3. then pick the entries, checking it doesn't exceed things
# 4. then do an onset from those terms/ or do that straight away

def create_synthetic_patient(phenotype_annotation: pl.DataFrame, num_disease: int, output_dir: Path):
    human_phenotype_ontology = load_ontology()
    omim_diseases = phenotype_annotation.filter(pl.col("database_id").str.starts_with("OMIM")).filter(
        pl.col("aspect") == "P")
    grouped_omim_diseases = omim_diseases.partition_by(by="database_id", maintain_order=True)
    if num_disease != 0:
        grouped_omim_diseases = random.sample(grouped_omim_diseases, num_disease)
    for omim_disease in grouped_omim_diseases:
        synthetic_patient_generator = SyntheticPatientGenerator(omim_disease, human_phenotype_ontology)
        patient_terms = synthetic_patient_generator.patient_term_annotation_set()
        phenopacket_file = PhenotypeAnnotationToPhenopacketConverter(
            human_phenotype_ontology
        ).create_phenopacket(patient_terms, synthetic_patient_generator.get_onset_range())
        write_phenopacket(
            phenopacket_file.phenopacket, output_dir.joinpath(phenopacket_file.phenopacket_path))
