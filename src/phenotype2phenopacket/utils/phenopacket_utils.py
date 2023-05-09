import json
import re
import secrets
from copy import copy
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Union

import polars as pl
from google.protobuf.json_format import MessageToJson, Parse
from google.protobuf.timestamp_pb2 import Timestamp
from phenopackets import (
    Age,
    Diagnosis,
    Disease,
    Family,
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

from phenotype2phenopacket.utils.gene_map_utils import GeneIdentifierUpdater
from phenotype2phenopacket.utils.utils import is_float


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


@dataclass
class OnsetTerm:
    lower_age: Union[int, str]
    upper_age: Union[int, str]


onset_hpo = {
    "HP:0011462": OnsetTerm(16, 40),
    "HP:0011460": OnsetTerm(0, 0.019230769230769232),
    "HP:0011463": OnsetTerm(1, 5),
    "HP:0003584": OnsetTerm(60, 90),
    "HP:0025709": OnsetTerm(19, 25),
    "HP:0034198": OnsetTerm(0, 0.019230769230769232),
    "HP:0003596": OnsetTerm(40, 60),
    "HP:0003621": OnsetTerm(5, 15),
    "HP:0003593": OnsetTerm(0, 1),
    "HP:4000040": OnsetTerm(0, 0),
    "HP:0011461": OnsetTerm(0, 0.019230769230769232),
    "HP:0003577": OnsetTerm(0, 0),
    "HP:0034197": OnsetTerm(0, 0.019230769230769232),
    "HP:0025708": OnsetTerm(16, 19),
    "HP:0410280": OnsetTerm(1, 15),
    "HP:0030674": OnsetTerm(0, 0.019230769230769232),
    "HP:0003623": OnsetTerm(0, 0.019230769230769232),
    "HP:0034199": OnsetTerm(0, 0.019230769230769232),
    "HP:0003581": OnsetTerm(16, 80),
    "HP:0025710": OnsetTerm(25, 40),
}


@dataclass
class PhenopacketFile:
    phenopacket: Phenopacket
    phenopacket_path: Path


def phenopacket_reader(file: Path):
    """Reads a phenopacket file, returning its contents."""
    file = open(file, "r")
    phenopacket = json.load(file)
    file.close()
    if "proband" in phenopacket:
        return Parse(json.dumps(phenopacket), Family())
    else:
        return Parse(json.dumps(phenopacket), Phenopacket())


def create_phenopacket_file_name_from_disease(disease_name: str) -> Path:
    normalised_string = re.sub(r"\W+", "_", disease_name)
    return Path(normalised_string.replace(" ", "_") + ".json")


def create_json_message(phenopacket: Phenopacket) -> str:
    """Create json message for writing to file."""
    return MessageToJson(phenopacket)


def write_phenopacket(phenopacket: Phenopacket, output_file: Path) -> None:
    """Write a phenopacket."""
    phenopacket_json = create_json_message(phenopacket)
    with open(output_file, "w") as outfile:
        outfile.write(phenopacket_json)
    outfile.close()


class SyntheticPatientGenerator:
    def __init__(self, disease_df: pl.DataFrame, ontology, ontology_factory):
        self.disease_df = disease_df
        self.ontology = ontology
        self.ontology_factory = ontology_factory
        self.lower_age = 0
        self.upper_age = 0
        self.filtered_df = []
        self.secret_rand = secrets.SystemRandom()

    def get_number_of_terms(self):
        """Get number of terms to ascribe from full set."""
        if len(self.disease_df) == 1:
            return 1
        return int(len(self.disease_df) * (self.secret_rand.uniform(0.2, 0.75)))

    @staticmethod
    def shuffle_dataframe(disease_df: pl.DataFrame):
        """Shuffle dataframe."""
        return disease_df.sample(fraction=1, shuffle=True)

    def add_frequency(self):
        """Add random frequency to annotations without one defined."""
        updated_df = []
        for row in self.disease_df.rows(named=True):
            if row["frequency"] is None:
                synthetic_frequency = self.secret_rand.uniform(0, 1)
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

    def check_hpo_frequency(self, phenotype_entry: dict):
        """Filter with HPO defined frequencies."""
        frequency_limits = frequency_hpo[phenotype_entry["frequency"]]
        random_frequency = self.secret_rand.uniform(0, 100)
        if (
            float(frequency_limits.lower) < random_frequency < float(frequency_limits.upper)
            and phenotype_entry not in self.filtered_df
        ):
            self.filtered_df.append(phenotype_entry)

    def check_frequency_threshold(
        self, frequency: float, phenotype_entry: dict, random_frequency: float
    ):
        """Check if patient frequency meets the filter for the disease frequency."""
        if random_frequency <= float(frequency) and phenotype_entry not in self.filtered_df:
            self.filtered_df.append(phenotype_entry)
        else:
            pass

    def check_percentage_frequency(self, phenotype_entry: dict):
        """Filter with percentage frequency."""
        frequency = phenotype_entry["frequency"].strip("%")
        random_frequency = self.secret_rand.uniform(0, 100)
        self.check_frequency_threshold(frequency, phenotype_entry, random_frequency)

    def check_fraction_frequency(self, phenotype_entry: dict):
        """Filter with fraction frequency."""
        random_frequency = self.secret_rand.uniform(0, 1)
        frequency = float(Fraction(phenotype_entry["frequency"]))
        self.check_frequency_threshold(frequency, phenotype_entry, random_frequency)

    def check_float_frequency(self, phenotype_entry: dict):
        """Filter with float threshold."""
        random_frequency = self.secret_rand.uniform(0, 1)
        frequency = float((phenotype_entry["frequency"]))
        self.check_frequency_threshold(frequency, phenotype_entry, random_frequency)

    def check_frequency(self, phenotype_entry: dict):
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
        return self.filter_phenotype_entries(
            self.shuffle_dataframe(self.add_frequency()), self.get_number_of_terms()
        )

    def get_number_of_terms_to_randomise(self, patient_terms: pl.DataFrame):
        """Get number of terms to randomise from filtered frequency set."""
        return self.secret_rand.randint(0, int(len(patient_terms)))

    def get_number_of_steps_for_randomisation(self):
        """Get the number of steps to take in range 1-5 for making a term more/less specific."""
        return self.secret_rand.randint(1, 5)

    def return_less_or_more_specific(self):
        """Generate a float between 0-1."""
        return self.secret_rand.random()

    def subsample_patient_terms(self, patient_terms: pl.dataframe):
        """Get a subsample of patient terms to make more/less specific."""
        return patient_terms.sample(self.get_number_of_terms_to_randomise(patient_terms))

    def get_children_of_term(self, phenotype_entry: dict, steps: int):
        """Get a child term of a hpo id from the number of steps specified."""
        term_id = phenotype_entry["hpo_id"]
        for _i in range(steps):
            descendants = self.ontology_factory.children(term_id)
            if descendants:
                term_id = self.secret_rand.choice(descendants)
            else:
                break
        phenotype_entry["hpo_id"] = term_id
        return phenotype_entry

    def get_parents_of_terms(self, phenotype_entry: dict, steps: int):
        """Get a parent term of a hpo id from the number of steps specified."""
        term_id = phenotype_entry["hpo_id"]
        rels = self.ontology.entity_alias_map(term_id)
        term = "".join(rels[(list(rels.keys())[0])])
        if term.startswith("Abnormality of"):
            return phenotype_entry
        for _i in range(steps):
            parents = self.ontology.hierarchical_parents(term_id)
            parent = self.secret_rand.choice(parents)
            rels = self.ontology.entity_alias_map(parent)
            term = "".join(rels[(list(rels.keys())[0])])
            if term.startswith("Abnormality of") or term_id == "HP:0000118":
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

    def alter_term_specificity(self, new_phenotype_terms: [dict], phenotype_entry: dict):
        """Alter the hpo id specificity - making less specific if the float is less than 0.5,
        otherwise the term is made more specific."""
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

    def patient_term_annotation_set(self):
        """Get the final patient term annotation set."""
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
        return pl.from_dicts(final_patient_terms)


class PhenotypeAnnotationToPhenopacketConverter:
    def __init__(self, human_phenotype_ontology):
        self.human_phenotype_ontology = human_phenotype_ontology
        self.secrets_random_num = secrets.SystemRandom()

    def create_individual(self, onset_range: OnsetTerm = None) -> Individual:
        """Create an Individual object."""
        age = (
            self.secrets_random_num.randint(onset_range.lower_age, onset_range.upper_age)
            if onset_range is not None
            else None
        )
        if onset_range is None or onset_range.upper_age == 0 and onset_range.lower_age == 0:
            age = None
        return Individual(
            id="patient1",
            time_at_last_encounter=TimeElement(age=Age(iso8601duration=f"P{age}Y"))
            if age is not None
            else None,
        )

    def create_onset(self, phenotype_annotation_entry: dict) -> TimeElement:
        """Create an Onset object."""
        if phenotype_annotation_entry["onset"] is not None:
            rels = self.human_phenotype_ontology.entity_alias_map(
                phenotype_annotation_entry["onset"]
            )
            term = "".join(rels[(list(rels.keys())[0])])
            return TimeElement(
                ontology_class=OntologyClass(id=phenotype_annotation_entry["onset"], label=term)
            )
        else:
            return None

    def create_modifier(self, phenotype_annotation_entry: dict) -> [OntologyClass]:
        """Create a Modifier."""
        if phenotype_annotation_entry["modifier"] is not None:
            try:
                rels = self.human_phenotype_ontology.entity_alias_map(
                    phenotype_annotation_entry["modifier"]
                )
                term = "".join(rels[(list(rels.keys())[0])])
                return [OntologyClass(id=phenotype_annotation_entry["modifier"], label=term)]
            except IndexError:
                return [OntologyClass(id=phenotype_annotation_entry["modifier"])]
        else:
            return None

    def create_phenotypic_feature(self, phenotype_annotation_entry: dict) -> PhenotypicFeature:
        """Create a PhenotypicFeature object."""
        if phenotype_annotation_entry["aspect"] == "P":
            rels = self.human_phenotype_ontology.entity_alias_map(
                phenotype_annotation_entry["hpo_id"]
            )
            hpo_term = "".join(rels[(list(rels.keys())[0])])
            return PhenotypicFeature(
                type=OntologyClass(id=phenotype_annotation_entry["hpo_id"], label=hpo_term),
                onset=self.create_onset(phenotype_annotation_entry),
                modifiers=self.create_modifier(phenotype_annotation_entry),
                excluded=True if phenotype_annotation_entry["qualifier"] == "NOT" else False,
            )
        else:
            return None

    def create_phenotypic_features(self, omim_disease_df: pl.DataFrame):
        """Create a list of phenotypic features."""
        phenotypic_features = []
        for phenotype_annotation_entry in omim_disease_df.rows(named=True):
            phenotypic_feature = self.create_phenotypic_feature(phenotype_annotation_entry)
            if phenotypic_feature is not None:
                phenotypic_features.append(phenotypic_feature)
        return phenotypic_features

    @staticmethod
    def create_disease(phenotype_annotation_entry: dict) -> Disease:
        """Create a Disease object."""
        return Disease(
            term=OntologyClass(
                id=phenotype_annotation_entry["database_id"],
                label=phenotype_annotation_entry["disease_name"],
            )
        )

    @staticmethod
    def create_omim_resource() -> Resource:
        """Create OMIM resource."""
        return Resource(
            id="omim",
            name="Online Mendelian Inheritance in Man",
            url="https://www.omim.org",
            version="hp/releases/2023-04-18",
            namespace_prefix="OMIM",
            iri_prefix="https://omim.org/entry/",
        )

    @staticmethod
    def create_human_phenotype_ontology_resource() -> Resource:
        """Create human phenotype ontology resource."""
        return Resource(
            id="hp",
            name="human phenotype ontology",
            url="http://purl.obolibrary.org/obo/hp.owl",
            version="hp/releases/2023-04-05",
            namespace_prefix="HP",
            iri_prefix="http://purl.obolibrary.org/obo/HP_",
        )

    def create_metadata(self) -> MetaData:
        """Create metadata"""
        timestamp = Timestamp()
        timestamp.GetCurrentTime()
        return MetaData(
            created=timestamp,
            created_by="phenotype2phenopacket",
            resources=[
                self.create_human_phenotype_ontology_resource(),
                self.create_omim_resource(),
            ],
            phenopacket_schema_version="2.0",
        )

    def create_phenopacket(
        self, omim_disease_df: pl.DataFrame, onset: OnsetTerm = None
    ) -> PhenopacketFile:
        """Create a Phenopacket object."""
        phenotypic_features = self.create_phenotypic_features(omim_disease_df)
        phenotype_annotation_entry = omim_disease_df.to_dicts()[0]
        return PhenopacketFile(
            phenopacket=Phenopacket(
                id=phenotype_annotation_entry["disease_name"].lower().replace(" ", "_"),
                subject=self.create_individual(onset),
                phenotypic_features=phenotypic_features,
                diseases=[self.create_disease(phenotype_annotation_entry)],
                meta_data=self.create_metadata(),
            ),
            phenopacket_path=create_phenopacket_file_name_from_disease(
                phenotype_annotation_entry["disease_name"]
            ),
        )


class PhenopacketUtil:
    def __init__(self, phenopacket: Phenopacket):
        self.phenopacket = phenopacket

    def return_phenopacket_disease(self) -> [str]:
        """Return the disease object from a phenopacket."""
        return self.phenopacket.diseases[0]


class PhenopacketInterpretationExtender:
    def __init__(self, phenopacket: Phenopacket):
        self.phenopacket = phenopacket

    @staticmethod
    def create_gene_genomic_interpretation(
        gene_to_phenotype_entry: dict, gene_identifier_updater: GeneIdentifierUpdater
    ):
        """Create genomic interpretation for a gene-to-phenotype relationship."""
        try:
            gene_symbol = gene_identifier_updater.obtain_gene_symbol_from_identifier(
                str(gene_to_phenotype_entry["entrez_id"])
            )
            return GenomicInterpretation(
                subject_or_biosample_id="patient1",
                interpretation_status=4
                if gene_to_phenotype_entry["disease_name"].startswith("?") is False
                else 0,
                gene=GeneDescriptor(
                    value_id=gene_identifier_updater.find_identifier(gene_symbol),
                    symbol=gene_symbol,
                ),
            )
        except KeyError:
            print(f"Unable to find gene_symbol for {gene_to_phenotype_entry['entrez_id']}")
            return None
        except TypeError:
            print("N/A value", gene_to_phenotype_entry)
            return None

    def create_gene_genomic_interpretations(
        self,
        omim_disease_phenotype_gene_map: pl.DataFrame,
        gene_identifier_updater: GeneIdentifierUpdater,
    ):
        """Create list of genomic interpretations for known gene-to-phenotype relationships."""
        genomic_interpretations = []
        for phenotype_entry in omim_disease_phenotype_gene_map.rows(named=True):
            genomic_interpretation = self.create_gene_genomic_interpretation(
                phenotype_entry, gene_identifier_updater
            )
            if genomic_interpretation is not None:
                genomic_interpretations.append(genomic_interpretation)
        return genomic_interpretations

    def create_gene_diagnosis(
        self,
        omim_disease_phenotype_gene_map: pl.DataFrame,
        gene_identifier_updater: GeneIdentifierUpdater,
        disease: Disease,
    ):
        """Create a diagnosis object for known gene-to-phenotype relationships."""
        genomic_interpretations = self.create_gene_genomic_interpretations(
            omim_disease_phenotype_gene_map, gene_identifier_updater
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
        omim_disease_phenotype_gene_map: pl.DataFrame,
        gene_identifier_updater: GeneIdentifierUpdater,
    ):
        """Create an interpretation object for known gene-to-phenotype relationships."""
        phenopacket_util = PhenopacketUtil(self.phenopacket)
        disease = phenopacket_util.return_phenopacket_disease()
        diagnosis = self.create_gene_diagnosis(
            omim_disease_phenotype_gene_map, gene_identifier_updater, disease
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
        omim_disease_phenotype_gene_map: pl.DataFrame,
        gene_identifier_updater: GeneIdentifierUpdater,
    ):
        """Add interpretations to phenopacket."""
        phenopacket_copy = copy(self.phenopacket)
        interpretations = [
            self.create_gene_interpretation(
                omim_disease_phenotype_gene_map, gene_identifier_updater
            )
        ]
        if interpretations is not None:
            phenopacket_copy.interpretations.extend(interpretations)
            return phenopacket_copy
