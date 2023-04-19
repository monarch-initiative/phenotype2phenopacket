import click

from phenotype2phenopacket.convert import convert_to_phenopackets_command


@click.group()
def main():
    """phenotype2phenopacket converter"""


main.add_command(convert_to_phenopackets_command)

if __name__ == "__main__":
    main()
