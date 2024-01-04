from collections import defaultdict
from pathlib import Path

from phenotype2phenopacket.utils.utils import read_hgnc_data


def create_hgnc_dict(hgnc_data_file: Path) -> defaultdict:
    """
    Create a dictionary as a reference for updating gene symbols and identifiers based on HGNC data.

    Args:
        hgnc_data_file (Path): Path to hgnc complete set file.

    Returns:
        defaultdict: A dictionary containing gene symbols as keys and their associated gene information.

    Notes:
        The dictionary structure:
        {
            'gene_symbol': {
                'ensembl_id': str,
                'hgnc_id': str,
                'entrez_id': str,
                'refseq_accession': str,
                'previous_symbol': [str, ...]
            },
            ...
        }
    """
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
    """
    Create a mapping of gene identifiers to gene symbols using HGNC data.

    Args:
        hgnc_data_file (Path): Path to the HGNC complete set file.

    Returns:
        dict: A mapping of gene identifiers to gene symbols.

    Notes:
        The dictionary structure:
        {
            'identifier': 'gene_symbol',
            ...
        }
    """
    hgnc_df = read_hgnc_data(hgnc_data_file)
    identifier_map = {}
    for _index, row in hgnc_df.iterrows():
        identifier_map[row["ensembl_gene_id"]] = row["symbol"]
        identifier_map[row["hgnc_id"]] = row["symbol"]
        identifier_map[row["entrez_id"]] = row["symbol"]
        identifier_map[row["refseq_accession"]] = row["symbol"]
    return identifier_map


class GeneIdentifierUpdater:
    """Class for updating gene identifiers within genomic interpretations."""

    def __init__(self, gene_identifier: str, hgnc_data: dict = None, identifier_map: dict = None):
        """
        Initialise the GeneIdentifierUpdater.

        Args:
            gene_identifier (str): The gene identifier to update to.
            hgnc_data (dict): A dictionary containing HGNC data (default: None).
            identifier_map (dict): A dictionary mapping gene identifiers (default: None).
        """
        self.hgnc_data = hgnc_data
        self.gene_identifier = gene_identifier
        self.identifier_map = identifier_map

    def find_identifier(self, gene_symbol: str) -> str:
        """
        Find the specified gene identifier for a gene symbol.

        Args:
            gene_symbol (str): The gene symbol to find the identifier for.

        Returns:
            str: The identified gene identifier.
        """
        if gene_symbol in self.hgnc_data.keys():
            return self.hgnc_data[gene_symbol][self.gene_identifier]
        else:
            for _symbol, data in self.hgnc_data.items():
                for prev_symbol in data["previous_symbol"]:
                    if prev_symbol == gene_symbol:
                        return data[self.gene_identifier]

    def obtain_gene_symbol_from_identifier(self, query_gene_identifier: str) -> str:
        """
        Obtain gene symbol from a gene identifier.

        Args:
            query_gene_identifier (str): The gene identifier.

        Returns:
            str: The gene symbol corresponding to the identifier.
        """
        return self.identifier_map[query_gene_identifier]
