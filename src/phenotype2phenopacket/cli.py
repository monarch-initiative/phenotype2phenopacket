import click

from phenotype2phenopacket.cli_add import add_genes_command, add_variants_command
from phenotype2phenopacket.cli_convert import convert_to_phenopackets_command


@click.group()
def main():
    """phenotype2phenopacket converter"""


main.add_command(convert_to_phenopackets_command)
main.add_command(add_variants_command)
main.add_command(add_genes_command)
if __name__ == "__main__":
    main()
