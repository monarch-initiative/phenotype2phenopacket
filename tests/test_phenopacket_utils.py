import unittest
from pathlib import Path
from unittest.mock import Mock, patch

import polars as pl
from polars.testing import assert_frame_equal

from phenotype2phenopacket.utils.phenopacket_utils import (
    OnsetTerm,
    SyntheticPatientGenerator,
    create_phenopacket_file_name_from_disease,
)
from phenotype2phenopacket.utils.utils import load_ontology, load_ontology_factory

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
            ontology_factory=load_ontology_factory(),
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
                    "hpo_id": "HP:0009789",
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
        with patch.object(self.synthetic_patient_generator, "return_less_or_more_specific",
                          return_value=mock_return_value), \
                patch.object(self.synthetic_patient_generator, 'get_number_of_steps_for_randomisation',
                             side_effect=[mock_steps_value, mock_steps_value]):
            altered_phenotype = self.synthetic_patient_generator.alter_term_specificity([], {
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
            })
        self.assertEqual(altered_phenotype, [{'database_id': 'OMIM:612567',
                                              'disease_name': 'Inflammatory bowel disease 25, '
                                                              'early onset, autosomal recessive',
                                              'qualifier': None, 'hpo_id': 'HP:0002583', 'reference': 'PMID:19890111',
                                              'evidence': 'PCS', 'onset': 'HP:0003593', 'frequency': '2/2', 'sex': None,
                                              'modifier': None, 'aspect': 'P',
                                              'biocuration': 'HPO:probinson[2013-03-12];HPO:probinson[2020-11-01]'}]
                         )

    def test_alter_term_specificity_more_specific(self):
        mock_return_value = 0.8
        mock_steps_value = 1
        with patch.object(self.synthetic_patient_generator, "return_less_or_more_specific",
                          return_value=mock_return_value), \
                patch.object(self.synthetic_patient_generator, 'get_number_of_steps_for_randomisation',
                             side_effect=[mock_steps_value, mock_steps_value]):
            altered_phenotype = self.synthetic_patient_generator.alter_term_specificity([{
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
            }], {
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
            })
        self.assertEqual(altered_phenotype, [{
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
        }, {'database_id': 'OMIM:612567',
            'disease_name': 'Inflammatory bowel disease 25, '
                            'early onset, autosomal recessive',
            'qualifier': None, 'hpo_id': 'HP:0003784', 'reference': 'PMID:19890111',
            'evidence': 'PCS', 'onset': 'HP:0010669', 'frequency': '2/2', 'sex': None,
            'modifier': None, 'aspect': 'P',
            'biocuration': 'HPO:probinson[2013-03-12];HPO:probinson[2020-11-01]'}])
