from collections import defaultdict
from pathlib import Path

from phenotype2phenopacket.utils.utils import read_hgnc_data


def create_hgnc_dict(hgnc_data_file: Path) -> defaultdict:
    """Creates reference for updating gene symbols and identifiers."""
    hgnc_df = read_hgnc_data(hgnc_data_file)
    hgnc_data = defaultdict(dict)
    for _index, row in hgnc_df.iterrows():
        previous_names = []
        hgnc_data[row["symbol"]]["ensembl_id"] = row["ensembl_gene_id"]
        hgnc_data[row["symbol"]]["hgnc_id"] = row["hgnc_id"]
        hgnc_data[row["symbol"]]["entrez_id"] = row["entrez_id"]
        hgnc_data[row["symbol"]]["refseq_accession"] = row["refseq_accession"]
        previous = str(row["prev_symbol"]).split("|")
        for p in previous:
            previous_names.append(p.strip('"'))
        hgnc_data[row["symbol"]]["previous_symbol"] = previous_names

    return hgnc_data


def create_gene_identifier_map(hgnc_data_file: Path) -> dict:
    hgnc_df = read_hgnc_data(hgnc_data_file)
    identifier_map = {}
    for _index, row in hgnc_df.iterrows():
        identifier_map[row["ensembl_gene_id"]] = row["symbol"]
        identifier_map[row["hgnc_id"]] = row["symbol"]
        identifier_map[row["entrez_id"]] = row["symbol"]
        identifier_map[row["refseq_accession"]] = row["symbol"]
    return identifier_map


class GeneIdentifierUpdater:
    def __init__(self, gene_identifier: str, hgnc_data: dict = None, identifier_map: dict = None):
        self.hgnc_data = hgnc_data
        self.gene_identifier = gene_identifier
        self.identifier_map = identifier_map

    def find_identifier(self, gene_symbol: str) -> str:
        """Finds the specified gene identifier for a gene symbol."""
        if gene_symbol in self.hgnc_data.keys():
            return self.hgnc_data[gene_symbol][self.gene_identifier]
        else:
            for _symbol, data in self.hgnc_data.items():
                for prev_symbol in data["previous_symbol"]:
                    if prev_symbol == gene_symbol:
                        return data[self.gene_identifier]

    def obtain_gene_symbol_from_identifier(self, query_gene_identifier: str) -> str:
        """
        Obtain gene symbol from a gene identifier. (e.g.)
        "
        obtain_gene_symbol_from_identifier(query_gene_identifier="HGNC:5")
        "
        """
        return self.identifier_map[query_gene_identifier]
