import unittest
from pathlib import Path
from unittest.mock import Mock, patch

import polars as pl
from phenopackets import (
    Age,
    Diagnosis,
    Disease,
    File,
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
from pheval.utils.phenopacket_utils import (
    GeneIdentifierUpdater,
    create_gene_identifier_map,
)
from polars.testing import assert_frame_equal

from phenotype2phenopacket.utils.phenopacket_utils import (
    OnsetTerm,
    PhenopacketInterpretationExtender,
    PhenopacketUtil,
    PhenotypeAnnotationToPhenopacketConverter,
    SyntheticPatientGenerator,
    create_phenopacket_file_name_from_disease,
)
from phenotype2phenopacket.utils.utils import load_ontology

disease_df = pl.from_dicts(
    [
        {
            "database_id": "OMIM:612567",
            "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
            "qualifier": None,
            "hpo_id": "HP:0000143",
            "reference": "PMID:19890111",
            "evidence": "PCS",
            "onset": "HP:0003593",
            "frequency": "1/1",
            "sex": "FEMALE",
            "modifier": None,
            "aspect": "P",
            "biocuration": "HPO:probinson[2013-03-12];HPO:probinson[2020-11-01]",
        },
        {
            "database_id": "OMIM:612567",
            "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
            "qualifier": None,
            "hpo_id": "HP:0004387",
            "reference": "PMID:19890111",
            "evidence": "PCS",
            "onset": "HP:0003593",
            "frequency": "2/2",
            "sex": None,
            "modifier": None,
            "aspect": "P",
            "biocuration": "HPO:probinson[2013-03-12];HPO:probinson[2020-11-01]",
        },
        {
            "database_id": "OMIM:612567",
            "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
            "qualifier": None,
            "hpo_id": "HP:0033279",
            "reference": "PMID:19890111",
            "evidence": "PCS",
            "onset": None,
            "frequency": "1/2",
            "sex": None,
            "modifier": None,
            "aspect": "P",
            "biocuration": "HPO:probinson[2020-12-07]",
        },
        {
            "database_id": "OMIM:612567",
            "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
            "qualifier": None,
            "hpo_id": "HP:0033256",
            "reference": "PMID:21519361",
            "evidence": "PCS",
            "onset": None,
            "frequency": "1/1",
            "sex": None,
            "modifier": None,
            "aspect": "P",
            "biocuration": "HPO:probinson[2022-02-27]",
        },
        {
            "database_id": "OMIM:612567",
            "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
            "qualifier": None,
            "hpo_id": "HP:0002837",
            "reference": "PMID:21519361",
            "evidence": "PCS",
            "onset": None,
            "frequency": "1/1",
            "sex": None,
            "modifier": None,
            "aspect": "P",
            "biocuration": "HPO:probinson[2022-02-27]",
        },
        {
            "database_id": "OMIM:612567",
            "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
            "qualifier": None,
            "hpo_id": "HP:0009789",
            "reference": "PMID:19890111",
            "evidence": "PCS",
            "onset": "HP:0003593",
            "frequency": "1/2",
            "sex": None,
            "modifier": None,
            "aspect": "P",
            "biocuration": "HPO:skoehler[2015-04-05];HPO:probinson[2020-11-01]",
        },
        {
            "database_id": "OMIM:612567",
            "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
            "qualifier": None,
            "hpo_id": "HP:0025084",
            "reference": "PMID:19890111",
            "evidence": "PCS",
            "onset": "HP:0003593",
            "frequency": "2/2",
            "sex": None,
            "modifier": None,
            "aspect": "P",
            "biocuration": "HPO:skoehler[2018-10-08];HPO:probinson[2020-11-01]",
        },
        {
            "database_id": "OMIM:612567",
            "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
            "qualifier": None,
            "hpo_id": "HP:0025084",
            "reference": "PMID:21519361",
            "evidence": "PCS",
            "onset": None,
            "frequency": "1/1",
            "sex": None,
            "modifier": None,
            "aspect": "P",
            "biocuration": "HPO:probinson[2022-02-27]",
        },
    ]
)

disease_df_with_frequency = pl.from_dicts(
    [
        {
            "database_id": "OMIM:612567",
            "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
            "qualifier": None,
            "hpo_id": "HP:0000143",
            "reference": "PMID:19890111",
            "evidence": "PCS",
            "onset": "HP:0003593",
            "frequency": "1/1",
            "sex": "FEMALE",
            "modifier": None,
            "aspect": "P",
            "biocuration": "HPO:probinson[2013-03-12];HPO:probinson[2020-11-01]",
        },
        {
            "database_id": "OMIM:612567",
            "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
            "qualifier": None,
            "hpo_id": "HP:0004387",
            "reference": "PMID:19890111",
            "evidence": "PCS",
            "onset": "HP:0003593",
            "frequency": "2/2",
            "sex": None,
            "modifier": None,
            "aspect": "P",
            "biocuration": "HPO:probinson[2013-03-12];HPO:probinson[2020-11-01]",
        },
        {
            "database_id": "OMIM:612567",
            "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
            "qualifier": None,
            "hpo_id": "HP:0033279",
            "reference": "PMID:19890111",
            "evidence": "PCS",
            "onset": None,
            "frequency": "1/2",
            "sex": None,
            "modifier": None,
            "aspect": "P",
            "biocuration": "HPO:probinson[2020-12-07]",
        },
    ]
)

gene_to_phenotypes = pl.from_dicts(
    [
        {"disease_id": "OMIM:612567", "ncbi_gene_id": "NCBIGene:3588", "gene_symbol": "IL10RB"},
        {
            "disease_id": "OMIM:612567",
            "ncbi_gene_id": "NCBIGene:1",
            "gene_symbol": "ADA",
        },
    ]
)
phenotypic_features_with_excluded = [
    PhenotypicFeature(type=OntologyClass(id="HP:0000256", label="Macrocephaly")),
    PhenotypicFeature(type=OntologyClass(id="HP:0002059", label="Cerebral atrophy")),
    PhenotypicFeature(type=OntologyClass(id="HP:0100309", label="Subdural hemorrhage")),
    PhenotypicFeature(type=OntologyClass(id="HP:0003150", label="Glutaric aciduria")),
    PhenotypicFeature(type=OntologyClass(id="HP:0001332", label="Dystonia")),
    PhenotypicFeature(
        type=OntologyClass(id="HP:0008494", label="Inferior lens subluxation"), excluded=True
    ),
]

phenopacket_files = [
    File(
        uri="test/path/to/test_1.vcf",
        file_attributes={"fileFormat": "VCF", "genomeAssembly": "GRCh37"},
    ),
    File(
        uri="test_1.ped",
        file_attributes={"fileFormat": "PED", "genomeAssembly": "GRCh37"},
    ),
]
phenopacket_metadata = MetaData(
    created_by="pheval-converter",
    resources=[
        Resource(
            id="hp",
            name="human phenotype ontology",
            url="http://purl.obolibrary.org/obo/hp.owl",
            version="hp/releases/2019-11-08",
            namespace_prefix="HP",
            iri_prefix="http://purl.obolibrary.org/obo/HP_",
        )
    ],
    phenopacket_schema_version="2.0",
)
phenopacket = Phenopacket(
    id="test-subject",
    subject=Individual(id="test-subject-1", sex=1),
    phenotypic_features=phenotypic_features_with_excluded,
    diseases=[
        Disease(
            term=OntologyClass(
                id="OMIM:612567",
                label="Inflammatory bowel disease 25, early onset, autosomal recessive",
            )
        )
    ],
    files=phenopacket_files,
    meta_data=phenopacket_metadata,
)


class TestCreatePhenopacketFileNameFromDisease(unittest.TestCase):
    def test_create_phenopacket_file_name_from_disease(self):
        self.assertEqual(
            create_phenopacket_file_name_from_disease(
                "Developmental and epileptic encephalopathy 96"
            ),
            Path("Developmental_and_epileptic_encephalopathy_96.json"),
        )
        self.assertEqual(
            create_phenopacket_file_name_from_disease(
                "Inflammatory bowel disease 25, early onset, autosomal recessive"
            ),
            Path("Inflammatory_bowel_disease_25_early_onset_autosomal_recessive.json"),
        )
        self.assertEqual(
            create_phenopacket_file_name_from_disease("Williams-Beuren syndrome"),
            Path("Williams_Beuren_syndrome.json"),
        )


class TestSyntheticPatientGenerator(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.synthetic_patient_generator = SyntheticPatientGenerator(
            disease_df=disease_df,
            ontology=load_ontology(),
        )

    def tearDown(self) -> None:
        self.synthetic_patient_generator.filtered_df = []

    def test_get_number_of_terms(self):
        self.assertEqual(type(self.synthetic_patient_generator.get_number_of_terms()), int)

    def test_shuffle_dataframe(self):
        self.assertEqual(
            len(self.synthetic_patient_generator.shuffle_dataframe(disease_df)), len(disease_df)
        )

    def test_add_frequency(self):
        assert_frame_equal(
            self.synthetic_patient_generator.add_frequency().select(
                pl.col("frequency").is_null().any()
            ),
            pl.DataFrame([{"frequency": False}]),
        )

    def test_get_onset_range(self):
        self.assertEqual(
            self.synthetic_patient_generator.get_onset_range(),
            OnsetTerm(lower_age=0, upper_age=1.0),
        )

    def test_check_hpo_frequency_passed(self):
        mock_rand = Mock()
        mock_rand.uniform.return_value = 50.0
        with patch.object(self.synthetic_patient_generator, "secret_rand", mock_rand):
            self.synthetic_patient_generator.check_hpo_frequency(
                {
                    "database_id": "OMIM:612567",
                    "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                    "qualifier": None,
                    "hpo_id": "HP:0025084",
                    "reference": "PMID:19890111",
                    "evidence": "PCS",
                    "onset": "HP:0003593",
                    "frequency": "HP:0040284",
                    "sex": None,
                    "modifier": None,
                    "aspect": "P",
                    "biocuration": "HPO:skoehler[2018-10-08];HPO:probinson[2020-11-01]",
                }
            )
        self.assertEqual(len(self.synthetic_patient_generator.filtered_df), 1)

    def test_check_hpo_frequency_failed(self):
        mock_rand = Mock()
        mock_rand.uniform.return_value = 50.0
        with patch.object(self.synthetic_patient_generator, "secret_rand", mock_rand):
            self.synthetic_patient_generator.check_hpo_frequency(
                {
                    "database_id": "OMIM:612567",
                    "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                    "qualifier": None,
                    "hpo_id": "HP:0025084",
                    "reference": "PMID:19890111",
                    "evidence": "PCS",
                    "onset": "HP:0003593",
                    "frequency": "HP:0040285",
                    "sex": None,
                    "modifier": None,
                    "aspect": "P",
                    "biocuration": "HPO:skoehler[2018-10-08];HPO:probinson[2020-11-01]",
                }
            )
        self.assertEqual(len(self.synthetic_patient_generator.filtered_df), 0)

    def test_check_frequency_threshold_percentage_passed(self):
        self.synthetic_patient_generator.check_frequency_threshold(
            73,
            {
                "database_id": "OMIM:612567",
                "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                "qualifier": None,
                "hpo_id": "HP:0025084",
                "reference": "PMID:19890111",
                "evidence": "PCS",
                "onset": "HP:0003593",
                "frequency": "73%",
                "sex": None,
                "modifier": None,
                "aspect": "P",
                "biocuration": "HPO:skoehler[2018-10-08];HPO:probinson[2020-11-01]",
            },
            17,
        )
        self.assertEqual(len(self.synthetic_patient_generator.filtered_df), 1)

    def test_check_frequency_threshold_percentage_failed(self):
        self.synthetic_patient_generator.check_frequency_threshold(
            17,
            {
                "database_id": "OMIM:612567",
                "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                "qualifier": None,
                "hpo_id": "HP:0025084",
                "reference": "PMID:19890111",
                "evidence": "PCS",
                "onset": "HP:0003593",
                "frequency": "17%",
                "sex": None,
                "modifier": None,
                "aspect": "P",
                "biocuration": "HPO:skoehler[2018-10-08];HPO:probinson[2020-11-01]",
            },
            73,
        )
        self.assertEqual(len(self.synthetic_patient_generator.filtered_df), 0)

    def test_check_frequency_threshold_fraction_passed(self):
        self.synthetic_patient_generator.check_frequency_threshold(
            0.834,
            {
                "database_id": "OMIM:612567",
                "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                "qualifier": None,
                "hpo_id": "HP:0025084",
                "reference": "PMID:19890111",
                "evidence": "PCS",
                "onset": "HP:0003593",
                "frequency": "17%",
                "sex": None,
                "modifier": None,
                "aspect": "P",
                "biocuration": "HPO:skoehler[2018-10-08];HPO:probinson[2020-11-01]",
            },
            0.2345,
        )
        self.assertEqual(len(self.synthetic_patient_generator.filtered_df), 1)

    def test_check_frequency_threshold_fraction_failed(self):
        self.synthetic_patient_generator.check_frequency_threshold(
            0.2345,
            {
                "database_id": "OMIM:612567",
                "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                "qualifier": None,
                "hpo_id": "HP:0025084",
                "reference": "PMID:19890111",
                "evidence": "PCS",
                "onset": "HP:0003593",
                "frequency": "17%",
                "sex": None,
                "modifier": None,
                "aspect": "P",
                "biocuration": "HPO:skoehler[2018-10-08];HPO:probinson[2020-11-01]",
            },
            0.9462,
        )
        self.assertEqual(len(self.synthetic_patient_generator.filtered_df), 0)

    def test_check_percentage_frequency_passed(self):
        mock_rand = Mock()
        mock_rand.uniform.return_value = 50.0
        with patch.object(self.synthetic_patient_generator, "secret_rand", mock_rand):
            self.synthetic_patient_generator.check_percentage_frequency(
                {
                    "database_id": "OMIM:612567",
                    "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                    "qualifier": None,
                    "hpo_id": "HP:0025084",
                    "reference": "PMID:19890111",
                    "evidence": "PCS",
                    "onset": "HP:0003593",
                    "frequency": "90%",
                    "sex": None,
                    "modifier": None,
                    "aspect": "P",
                    "biocuration": "HPO:skoehler[2018-10-08];HPO:probinson[2020-11-01]",
                }
            )
        self.assertEqual(len(self.synthetic_patient_generator.filtered_df), 1)

    def test_check_percentage_frequency_failed(self):
        mock_rand = Mock()
        mock_rand.uniform.return_value = 50.0
        with patch.object(self.synthetic_patient_generator, "secret_rand", mock_rand):
            self.synthetic_patient_generator.check_percentage_frequency(
                {
                    "database_id": "OMIM:612567",
                    "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                    "qualifier": None,
                    "hpo_id": "HP:0025084",
                    "reference": "PMID:19890111",
                    "evidence": "PCS",
                    "onset": "HP:0003593",
                    "frequency": "10%",
                    "sex": None,
                    "modifier": None,
                    "aspect": "P",
                    "biocuration": "HPO:skoehler[2018-10-08];HPO:probinson[2020-11-01]",
                }
            )
        self.assertEqual(len(self.synthetic_patient_generator.filtered_df), 0)

    def test_check_fraction_frequency_passed(self):
        mock_rand = Mock()
        mock_rand.uniform.return_value = 0.21345
        with patch.object(self.synthetic_patient_generator, "secret_rand", mock_rand):
            self.synthetic_patient_generator.check_fraction_frequency(
                {
                    "database_id": "OMIM:612567",
                    "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                    "qualifier": None,
                    "hpo_id": "HP:0025084",
                    "reference": "PMID:19890111",
                    "evidence": "PCS",
                    "onset": "HP:0003593",
                    "frequency": "1/3",
                    "sex": None,
                    "modifier": None,
                    "aspect": "P",
                    "biocuration": "HPO:skoehler[2018-10-08];HPO:probinson[2020-11-01]",
                }
            )
        self.assertEqual(len(self.synthetic_patient_generator.filtered_df), 1)

    def test_check_fraction_frequency_failed(self):
        mock_rand = Mock()
        mock_rand.uniform.return_value = 0.93452
        with patch.object(self.synthetic_patient_generator, "secret_rand", mock_rand):
            self.synthetic_patient_generator.check_fraction_frequency(
                {
                    "database_id": "OMIM:612567",
                    "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                    "qualifier": None,
                    "hpo_id": "HP:0025084",
                    "reference": "PMID:19890111",
                    "evidence": "PCS",
                    "onset": "HP:0003593",
                    "frequency": "1/3",
                    "sex": None,
                    "modifier": None,
                    "aspect": "P",
                    "biocuration": "HPO:skoehler[2018-10-08];HPO:probinson[2020-11-01]",
                }
            )
        self.assertEqual(len(self.synthetic_patient_generator.filtered_df), 0)

    def test_check_float_frequency_passed(self):
        mock_rand = Mock()
        mock_rand.uniform.return_value = 0.93452
        with patch.object(self.synthetic_patient_generator, "secret_rand", mock_rand):
            self.synthetic_patient_generator.check_float_frequency(
                {
                    "database_id": "OMIM:612567",
                    "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                    "qualifier": None,
                    "hpo_id": "HP:0025084",
                    "reference": "PMID:19890111",
                    "evidence": "PCS",
                    "onset": "HP:0003593",
                    "frequency": 0.98234,
                    "sex": None,
                    "modifier": None,
                    "aspect": "P",
                    "biocuration": "HPO:skoehler[2018-10-08];HPO:probinson[2020-11-01]",
                }
            )
        self.assertEqual(len(self.synthetic_patient_generator.filtered_df), 1)

    def test_check_float_frequency_failed(self):
        mock_rand = Mock()
        mock_rand.uniform.return_value = 0.93452
        with patch.object(self.synthetic_patient_generator, "secret_rand", mock_rand):
            self.synthetic_patient_generator.check_float_frequency(
                {
                    "database_id": "OMIM:612567",
                    "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                    "qualifier": None,
                    "hpo_id": "HP:0025084",
                    "reference": "PMID:19890111",
                    "evidence": "PCS",
                    "onset": "HP:0003593",
                    "frequency": 0.2345,
                    "sex": None,
                    "modifier": None,
                    "aspect": "P",
                    "biocuration": "HPO:skoehler[2018-10-08];HPO:probinson[2020-11-01]",
                }
            )
        self.assertEqual(len(self.synthetic_patient_generator.filtered_df), 0)

    def test_check_frequency_hpo(self):
        mock_rand = Mock()
        mock_rand.uniform.return_value = 50.0
        with patch.object(self.synthetic_patient_generator, "secret_rand", mock_rand):
            self.synthetic_patient_generator.check_frequency(
                {
                    "database_id": "OMIM:612567",
                    "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                    "qualifier": None,
                    "hpo_id": "HP:0025084",
                    "reference": "PMID:19890111",
                    "evidence": "PCS",
                    "onset": "HP:0003593",
                    "frequency": "HP:0040284",
                    "sex": None,
                    "modifier": None,
                    "aspect": "P",
                    "biocuration": "HPO:skoehler[2018-10-08];HPO:probinson[2020-11-01]",
                }
            )
        self.assertEqual(len(self.synthetic_patient_generator.filtered_df), 1)

    def test_check_frequency_percentage(self):
        mock_rand = Mock()
        mock_rand.uniform.return_value = 50.0
        with patch.object(self.synthetic_patient_generator, "secret_rand", mock_rand):
            self.synthetic_patient_generator.check_frequency(
                {
                    "database_id": "OMIM:612567",
                    "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                    "qualifier": None,
                    "hpo_id": "HP:0025084",
                    "reference": "PMID:19890111",
                    "evidence": "PCS",
                    "onset": "HP:0003593",
                    "frequency": "90%",
                    "sex": None,
                    "modifier": None,
                    "aspect": "P",
                    "biocuration": "HPO:skoehler[2018-10-08];HPO:probinson[2020-11-01]",
                }
            )
        self.assertEqual(len(self.synthetic_patient_generator.filtered_df), 1)

    def test_check_frequency_fraction(self):
        mock_rand = Mock()
        mock_rand.uniform.return_value = 0.21345
        with patch.object(self.synthetic_patient_generator, "secret_rand", mock_rand):
            self.synthetic_patient_generator.check_fraction_frequency(
                {
                    "database_id": "OMIM:612567",
                    "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                    "qualifier": None,
                    "hpo_id": "HP:0025084",
                    "reference": "PMID:19890111",
                    "evidence": "PCS",
                    "onset": "HP:0003593",
                    "frequency": "1/3",
                    "sex": None,
                    "modifier": None,
                    "aspect": "P",
                    "biocuration": "HPO:skoehler[2018-10-08];HPO:probinson[2020-11-01]",
                }
            )
        self.assertEqual(len(self.synthetic_patient_generator.filtered_df), 1)

    def test_check_frequency_float(self):
        mock_rand = Mock()
        mock_rand.uniform.return_value = 0.93452
        with patch.object(self.synthetic_patient_generator, "secret_rand", mock_rand):
            self.synthetic_patient_generator.check_float_frequency(
                {
                    "database_id": "OMIM:612567",
                    "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                    "qualifier": None,
                    "hpo_id": "HP:0025084",
                    "reference": "PMID:19890111",
                    "evidence": "PCS",
                    "onset": "HP:0003593",
                    "frequency": 0.98234,
                    "sex": None,
                    "modifier": None,
                    "aspect": "P",
                    "biocuration": "HPO:skoehler[2018-10-08];HPO:probinson[2020-11-01]",
                }
            )
        self.assertEqual(len(self.synthetic_patient_generator.filtered_df), 1)

    def test_filter_phenotype_entries(self):
        self.assertEqual(
            len(
                self.synthetic_patient_generator.filter_phenotype_entries(
                    disease_df_with_frequency, 1
                )
            ),
            1,
        )

    def test_get_children_of_term_no_child(self):
        self.assertEqual(
            self.synthetic_patient_generator.get_children_of_term(
                {
                    "database_id": "OMIM:612567",
                    "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                    "qualifier": None,
                    "hpo_id": "HP:0020084",
                    "reference": "PMID:21519361",
                    "evidence": "PCS",
                    "onset": None,
                    "frequency": "1/1",
                    "sex": None,
                    "modifier": None,
                    "aspect": "P",
                    "biocuration": "HPO:probinson[2022-02-27]",
                },
                5,
            ),
            {
                "database_id": "OMIM:612567",
                "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                "qualifier": None,
                "hpo_id": "HP:0020084",
                "reference": "PMID:21519361",
                "evidence": "PCS",
                "onset": None,
                "frequency": "1/1",
                "sex": None,
                "modifier": None,
                "aspect": "P",
                "biocuration": "HPO:probinson[2022-02-27]",
            },
        )

    def test_get_children_of_term(self):
        def mock_choice(lst):
            return lst[0]

        with patch.object(self.synthetic_patient_generator.secret_rand, "choice", mock_choice):
            self.assertEqual(
                self.synthetic_patient_generator.get_children_of_term(
                    {
                        "database_id": "OMIM:612567",
                        "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                        "qualifier": None,
                        "hpo_id": "HP:0011123",
                        "reference": "PMID:21519361",
                        "evidence": "PCS",
                        "onset": None,
                        "frequency": "1/1",
                        "sex": None,
                        "modifier": None,
                        "aspect": "P",
                        "biocuration": "HPO:probinson[2022-02-27]",
                    },
                    2,
                ),
                {
                    "database_id": "OMIM:612567",
                    "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                    "qualifier": None,
                    "hpo_id": "HP:0001047",
                    "reference": "PMID:21519361",
                    "evidence": "PCS",
                    "onset": None,
                    "frequency": "1/1",
                    "sex": None,
                    "modifier": None,
                    "aspect": "P",
                    "biocuration": "HPO:probinson[2022-02-27]",
                },
            )

    def test_get_parents_of_terms(self):
        self.assertEqual(
            self.synthetic_patient_generator.get_parents_of_terms(
                {
                    "database_id": "OMIM:612567",
                    "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                    "qualifier": None,
                    "hpo_id": "HP:0410218",
                    "reference": "PMID:21519361",
                    "evidence": "PCS",
                    "onset": None,
                    "frequency": "1/1",
                    "sex": None,
                    "modifier": None,
                    "aspect": "P",
                    "biocuration": "HPO:probinson[2022-02-27]",
                },
                4,
            ),
            {
                "database_id": "OMIM:612567",
                "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                "qualifier": None,
                "hpo_id": "HP:0011821",
                "reference": "PMID:21519361",
                "evidence": "PCS",
                "onset": None,
                "frequency": "1/1",
                "sex": None,
                "modifier": None,
                "aspect": "P",
                "biocuration": "HPO:probinson[2022-02-27]",
            },
        )

    def test_remove_terms_to_be_randomised(self):
        self.assertEqual(
            self.synthetic_patient_generator.remove_terms_to_be_randomised(
                pl.from_dicts(
                    [
                        {
                            "database_id": "OMIM:612567",
                            "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                            "qualifier": None,
                            "hpo_id": "HP:0000143",
                            "reference": "PMID:19890111",
                            "evidence": "PCS",
                            "onset": "HP:0003593",
                            "frequency": "1/1",
                            "sex": "FEMALE",
                            "modifier": None,
                            "aspect": "P",
                            "biocuration": "HPO:probinson[2013-03-12];HPO:probinson[2020-11-01]",
                        },
                        {
                            "database_id": "OMIM:612567",
                            "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                            "qualifier": None,
                            "hpo_id": "HP:0004387",
                            "reference": "PMID:19890111",
                            "evidence": "PCS",
                            "onset": "HP:0003593",
                            "frequency": "2/2",
                            "sex": None,
                            "modifier": None,
                            "aspect": "P",
                            "biocuration": "HPO:probinson[2013-03-12];HPO:probinson[2020-11-01]",
                        },
                    ]
                ),
                pl.from_dicts(
                    [
                        {
                            "database_id": "OMIM:612567",
                            "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                            "qualifier": None,
                            "hpo_id": "HP:0000143",
                            "reference": "PMID:19890111",
                            "evidence": "PCS",
                            "onset": "HP:0003593",
                            "frequency": "1/1",
                            "sex": "FEMALE",
                            "modifier": None,
                            "aspect": "P",
                            "biocuration": "HPO:probinson[2013-03-12];HPO:probinson[2020-11-01]",
                        }
                    ]
                ),
            ).to_dicts(),
            [
                {
                    "database_id": "OMIM:612567",
                    "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                    "qualifier": None,
                    "hpo_id": "HP:0004387",
                    "reference": "PMID:19890111",
                    "evidence": "PCS",
                    "onset": "HP:0003593",
                    "frequency": "2/2",
                    "sex": None,
                    "modifier": None,
                    "aspect": "P",
                    "biocuration": "HPO:probinson[2013-03-12];HPO:probinson[2020-11-01]",
                }
            ],
        )

    def test_alter_term_specificity_less_specific(self):
        mock_return_value = 0.4
        mock_steps_value = 1
        with (
            patch.object(
                self.synthetic_patient_generator,
                "return_less_or_more_specific",
                return_value=mock_return_value,
            ),
            patch.object(
                self.synthetic_patient_generator,
                "get_number_of_steps_for_randomisation",
                side_effect=[mock_steps_value, mock_steps_value],
            ),
        ):
            altered_phenotype = self.synthetic_patient_generator.alter_term_specificity(
                [],
                {
                    "database_id": "OMIM:612567",
                    "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                    "qualifier": None,
                    "hpo_id": "HP:0004387",
                    "reference": "PMID:19890111",
                    "evidence": "PCS",
                    "onset": "HP:0003593",
                    "frequency": "2/2",
                    "sex": None,
                    "modifier": None,
                    "aspect": "P",
                    "biocuration": "HPO:probinson[2013-03-12];HPO:probinson[2020-11-01]",
                },
            )
        self.assertEqual(
            altered_phenotype,
            [
                {
                    "database_id": "OMIM:612567",
                    "disease_name": "Inflammatory bowel disease 25, "
                    "early onset, autosomal recessive",
                    "qualifier": None,
                    "hpo_id": "HP:0002583",
                    "reference": "PMID:19890111",
                    "evidence": "PCS",
                    "onset": "HP:0003593",
                    "frequency": "2/2",
                    "sex": None,
                    "modifier": None,
                    "aspect": "P",
                    "biocuration": "HPO:probinson[2013-03-12];HPO:probinson[2020-11-01]",
                }
            ],
        )

    def test_alter_term_specificity_more_specific(self):
        mock_return_value = 0.8
        mock_steps_value = 1
        with (
            patch.object(
                self.synthetic_patient_generator,
                "return_less_or_more_specific",
                return_value=mock_return_value,
            ),
            patch.object(
                self.synthetic_patient_generator,
                "get_number_of_steps_for_randomisation",
                side_effect=[mock_steps_value, mock_steps_value],
            ),
        ):
            altered_phenotype = self.synthetic_patient_generator.alter_term_specificity(
                [
                    {
                        "database_id": "OMIM:612567",
                        "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                        "qualifier": None,
                        "hpo_id": "HP:0004387",
                        "reference": "PMID:19890111",
                        "evidence": "PCS",
                        "onset": "HP:0003593",
                        "frequency": "2/2",
                        "sex": None,
                        "modifier": None,
                        "aspect": "P",
                        "biocuration": "HPO:probinson[2013-03-12];HPO:probinson[2020-11-01]",
                    }
                ],
                {
                    "database_id": "OMIM:612567",
                    "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                    "qualifier": None,
                    "hpo_id": "HP:0011862",
                    "reference": "PMID:19890111",
                    "evidence": "PCS",
                    "onset": "HP:0010669",
                    "frequency": "2/2",
                    "sex": None,
                    "modifier": None,
                    "aspect": "P",
                    "biocuration": "HPO:probinson[2013-03-12];HPO:probinson[2020-11-01]",
                },
            )
        self.assertEqual(
            altered_phenotype,
            [
                {
                    "database_id": "OMIM:612567",
                    "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
                    "qualifier": None,
                    "hpo_id": "HP:0004387",
                    "reference": "PMID:19890111",
                    "evidence": "PCS",
                    "onset": "HP:0003593",
                    "frequency": "2/2",
                    "sex": None,
                    "modifier": None,
                    "aspect": "P",
                    "biocuration": "HPO:probinson[2013-03-12];HPO:probinson[2020-11-01]",
                },
                {
                    "database_id": "OMIM:612567",
                    "disease_name": "Inflammatory bowel disease 25, "
                    "early onset, autosomal recessive",
                    "qualifier": None,
                    "hpo_id": "HP:0003784",
                    "reference": "PMID:19890111",
                    "evidence": "PCS",
                    "onset": "HP:0010669",
                    "frequency": "2/2",
                    "sex": None,
                    "modifier": None,
                    "aspect": "P",
                    "biocuration": "HPO:probinson[2013-03-12];HPO:probinson[2020-11-01]",
                },
            ],
        )


class TestPhenotypeAnnotationToPhenopacketConverter(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.converter = PhenotypeAnnotationToPhenopacketConverter(load_ontology())
        cls.phenotype_entry = {
            "database_id": "OMIM:612567",
            "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
            "qualifier": None,
            "hpo_id": "HP:0004387",
            "reference": "PMID:19890111",
            "evidence": "PCS",
            "onset": "HP:0003593",
            "frequency": "2/2",
            "sex": None,
            "modifier": "HP:0012828",
            "aspect": "P",
            "biocuration": "HPO:probinson[2013-03-12];HPO:probinson[2020-11-01]",
        }
        cls.basic_phenotype_entry = {
            "database_id": "OMIM:612567",
            "disease_name": "Inflammatory bowel disease 25, early onset, autosomal recessive",
            "qualifier": None,
            "hpo_id": "HP:0004387",
            "reference": "PMID:19890111",
            "evidence": "PCS",
            "onset": None,
            "frequency": None,
            "sex": None,
            "modifier": None,
            "aspect": "P",
            "biocuration": "HPO:probinson[2013-03-12];HPO:probinson[2020-11-01]",
        }

    def test_create_individual_no_onset(self):
        self.assertEqual(self.converter.create_individual(), Individual(id="patient1"))

    def test_create_individual_with_onset(self):
        mock_rand = Mock()
        mock_rand.randint.return_value = 65
        with patch.object(self.converter, "secrets_random_num", mock_rand):
            self.assertEqual(
                self.converter.create_individual(onset_range=OnsetTerm(upper_age=80, lower_age=40)),
                Individual(
                    id="patient1",
                    time_at_last_encounter=TimeElement(age=Age(iso8601duration="P65Y")),
                ),
            )

    def test_create_onset(self):
        self.assertEqual(
            self.converter.create_onset(self.phenotype_entry),
            TimeElement(ontology_class=OntologyClass(id="HP:0003593", label="Infantile onset")),
        )

    def test_create_onset_none(self):
        self.assertEqual(self.converter.create_onset(self.basic_phenotype_entry), None)

    def test_create_modifier(self):
        self.assertEqual(
            self.converter.create_modifier(self.phenotype_entry),
            [OntologyClass(id="HP:0012828", label="Severe")],
        )

    def test_create_modifier_none(self):
        self.assertEqual(self.converter.create_modifier(self.basic_phenotype_entry), None)

    def test_create_phenotypic_feature(self):
        self.assertEqual(
            self.converter.create_phenotypic_feature(self.phenotype_entry),
            PhenotypicFeature(
                type=OntologyClass(id="HP:0004387", label="Enterocolitis"),
                onset=TimeElement(
                    ontology_class=OntologyClass(id="HP:0003593", label="Infantile onset")
                ),
                modifiers=[OntologyClass(id="HP:0012828", label="Severe")],
                excluded=False,
            ),
        )

    def test_create_phenotypic_feature_basic(self):
        self.assertEqual(
            self.converter.create_phenotypic_feature(self.basic_phenotype_entry),
            PhenotypicFeature(
                type=OntologyClass(id="HP:0004387", label="Enterocolitis"),
                onset=None,
                modifiers=None,
                excluded=False,
            ),
        )

    def test_create_phenotypic_features(self):
        self.assertEqual(
            self.converter.create_phenotypic_features(disease_df_with_frequency),
            [
                PhenotypicFeature(
                    type=OntologyClass(id="HP:0000143", label="Rectovaginal fistula"),
                    onset=TimeElement(
                        ontology_class=OntologyClass(id="HP:0003593", label="Infantile onset")
                    ),
                    modifiers=None,
                    excluded=False,
                ),
                PhenotypicFeature(
                    type=OntologyClass(id="HP:0004387", label="Enterocolitis"),
                    onset=TimeElement(
                        ontology_class=OntologyClass(id="HP:0003593", label="Infantile onset")
                    ),
                    modifiers=None,
                    excluded=False,
                ),
                PhenotypicFeature(
                    type=OntologyClass(id="HP:0033279", label="Enterocutaneous fistula"),
                    onset=None,
                    modifiers=None,
                    excluded=False,
                ),
            ],
        )

    def test_create_disease(self):
        self.assertEqual(
            self.converter.create_disease(self.phenotype_entry),
            Disease(
                term=OntologyClass(
                    id="OMIM:612567",
                    label="Inflammatory bowel disease 25, early onset, autosomal recessive",
                )
            ),
        )

    def test_create_omim_resource(self):
        self.assertEqual(
            self.converter.create_omim_resource(),
            Resource(
                id="omim",
                name="Online Mendelian Inheritance in Man",
                url="https://www.omim.org",
                version="hp/releases/2023-04-18",
                namespace_prefix="OMIM",
                iri_prefix="https://omim.org/entry/",
            ),
        )

    def test_create_human_phenotype_ontology_resource(self):
        self.assertEqual(
            self.converter.create_human_phenotype_ontology_resource("2023-04-05"),
            Resource(
                id="hp",
                name="human phenotype ontology",
                url="http://purl.obolibrary.org/obo/hp.owl",
                version="hp/releases/2023-04-05",
                namespace_prefix="HP",
                iri_prefix="http://purl.obolibrary.org/obo/HP_",
            ),
        )


class TestPhenopacketUtil(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.phenopacket = PhenopacketUtil(phenopacket)

    def test_return_phenopacket_disease(self):
        self.assertEqual(
            self.phenopacket.return_phenopacket_disease(),
            Disease(
                term=OntologyClass(
                    id="OMIM:612567",
                    label="Inflammatory bowel disease 25, early onset, autosomal recessive",
                )
            ),
        )


class TestPhenopacketInterpretationExtender(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.phenopacket = PhenopacketInterpretationExtender(phenopacket)
        cls.phenotype_to_gene_entry = {
            "disease_id": "OMIM:612567",
            "ncbi_gene_id": "NCBIGene:3588",
            "gene_symbol": "IL10RB",
        }
        cls.gene_identifier_updater = GeneIdentifierUpdater(
            gene_identifier="ensembl_id",
            identifier_map=create_gene_identifier_map(),
        )

    def test_create_gene_genomic_interpretation(self):
        self.assertEqual(
            self.phenopacket.create_gene_genomic_interpretation(
                self.phenotype_to_gene_entry, self.gene_identifier_updater
            ),
            GenomicInterpretation(
                subject_or_biosample_id="patient1",
                interpretation_status=4,
                gene=GeneDescriptor(
                    value_id="ENSG00000243646",
                    symbol="IL10RB",
                ),
            ),
        )

    def test_create_gene_genomic_interpretations(self):
        self.assertEqual(
            self.phenopacket.create_gene_genomic_interpretations(
                gene_to_phenotypes, self.gene_identifier_updater
            ),
            [
                GenomicInterpretation(
                    subject_or_biosample_id="patient1",
                    interpretation_status=4,
                    gene=GeneDescriptor(
                        value_id="ENSG00000243646",
                        symbol="IL10RB",
                    ),
                ),
                GenomicInterpretation(
                    subject_or_biosample_id="patient1",
                    interpretation_status=4,
                    gene=GeneDescriptor(value_id="ENSG00000196839", symbol="ADA"),
                ),
            ],
        )

    def test_create_gene_diagnosis(self):
        self.assertEqual(
            self.phenopacket.create_gene_diagnosis(
                gene_to_phenotypes,
                self.gene_identifier_updater,
                Disease(
                    term=OntologyClass(
                        id="OMIM:612567",
                        label="Inflammatory bowel disease 25, early onset, autosomal recessive",
                    )
                ),
            ),
            Diagnosis(
                disease=OntologyClass(
                    id="OMIM:612567",
                    label="Inflammatory bowel disease 25, early onset, autosomal recessive",
                ),
                genomic_interpretations=[
                    GenomicInterpretation(
                        subject_or_biosample_id="patient1",
                        interpretation_status=4,
                        gene=GeneDescriptor(
                            value_id="ENSG00000243646",
                            symbol="IL10RB",
                        ),
                    ),
                    GenomicInterpretation(
                        subject_or_biosample_id="patient1",
                        interpretation_status=4,
                        gene=GeneDescriptor(value_id="ENSG00000196839", symbol="ADA"),
                    ),
                ],
            ),
        )

    def test_create_gene_interpretation(self):
        self.assertEqual(
            self.phenopacket.create_gene_interpretation(
                gene_to_phenotypes, self.gene_identifier_updater
            ),
            Interpretation(
                id="Inflammatory bowel disease 25, early onset, autosomal recessive"
                + "-interpretation",
                progress_status=0,
                diagnosis=Diagnosis(
                    disease=OntologyClass(
                        id="OMIM:612567",
                        label="Inflammatory bowel disease 25, early onset, autosomal recessive",
                    ),
                    genomic_interpretations=[
                        GenomicInterpretation(
                            subject_or_biosample_id="patient1",
                            interpretation_status=4,
                            gene=GeneDescriptor(
                                value_id="ENSG00000243646",
                                symbol="IL10RB",
                            ),
                        ),
                        GenomicInterpretation(
                            subject_or_biosample_id="patient1",
                            interpretation_status=4,
                            gene=GeneDescriptor(value_id="ENSG00000196839", symbol="ADA"),
                        ),
                    ],
                ),
            ),
        )
