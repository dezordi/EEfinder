"""Console entrypoint: the top-level ``eefinder`` command group."""

import click
from eefinder.scripts.get_databases import get_databases
from eefinder.scripts.screening import screening
from eefinder import __version__


@click.group()
@click.version_option(__version__)
def cli():
    """EEfinder: find Endogenous Elements in eukaryote genomes.

    Use ``screening`` to run the identify endogenous elements and ``get-databases``
    to download the RefSeq protein databases it needs.
    """


cli.add_command(get_databases)
cli.add_command(screening)
