# Phenotype2Phenopacket

_Phenotype2Phenopacket_ is a command-line interface (CLI) application for the construction of phenopackets
from a phenotype annotation file.

## Installation

Phenotype2Phenopacket can be installed from PyPi.

```shell
pip install phenotype2phenopacket
```

## Usages

To convert all OMIM diseases in a phenotype annotation file to disease phenopackets, where all phenotypes are retained:

```shell
p2p convert --phenotype-annotation /path/to/phenotype.hpoa --output-dir /path/to/output-dir
```

To create synthetic patient disease phenopackets, where the dataset is more variable and frequencies are taken
into account and constrained noise is applied :

```shell
p2p create --phenotype-annotation /path/to/phenotype.hpoa --output-dir /path/to/output-dir
```

You can also limit the number of disease phenopackets converted/created:

```shell
p2p convert --phenotype-annotation /path/to/phenotype.hpoa --output-dir /path/to/output-dir --num-diseases 100
```

Or limit for a specific OMIM disease:

```shell
p2p create --phenotype-annotation /path/to/phenotype.hpoa --output-dir /path/to/output-dir --omim-id OMIM:619340
```

Or limit for a list of OMIM IDs specified in a text file, with each ID separated by a new line:

```shell
p2p create --phenotype-annotation /path/to/phenotype.hpoa --output-dir /path/to/output-dir --omim-id-list /path/to/list.txt
```

To add known gene-to-phenotype relationships to phenopackets:

```shell
p2p add-genes --phenopacket-dir /path/to/synthetic-phenopackets --genes-to-disease /path/to/genes_to_disease.txt --hgnc-data /path/to/hgnc_complete_set.txt --output-dir /path/to/output-dir
```

> **_NOTE:_** To add known gene-to-phenotype the genes_to_disease.txt is expected. It can be downloaded [here](https://hpo.jax.org/data/annotations).