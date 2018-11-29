import os

import click

from taxies import gen


@click.group()
def cli():
    """Taxies - rapid and accurate taxonomy predictions for microbial genomics
    """
    pass

@cli.command()
@click.argument('genome', type=click.Path(exists=True))
@click.argument('db', type=click.Path(exists=True, file_okay=False))
@click.argument('mark', type=click.Path(exists=False, writable=True))
@click.option('-e', '--max-evalue', type=click.FLOAT, default=0.001,
              show_default=True, help="E-value threshold")
@click.option('-t', '--threads', type=click.INT, default=1,
              show_default=True, help="number of threads")
def genmark(genome, db, mark, max_evalue, threads):
    """Extract markers from assembled genomes"""
    gen.mark(genome, db, mark, max_evalue=max_evalue, threads=threads)


@cli.command()
@click.argument('mark', type=click.Path(exists=True))
@click.argument('db', type=click.Path(exists=True, file_okay=False))
@click.argument('toph', type=click.Path(exists=False, writable=True))
@click.option('-p', '--pnames', type=click.Path(exists=True),
    help="custom protein names file")
@click.option('-d', '--min-ident', type=click.FLOAT, default=0.60,
              show_default=True, help="minimum top hit identity")
@click.option('-t', '--threads', type=click.INT, default=1,
              show_default=True, help="number of threads")

def gentoph(mark, db, toph, pnames, min_ident, threads):
    """Identify the top-hit allele for each marker"""

    gen.toph(mark, db, toph, custom_pnames_fn=pnames, min_ident=min_ident,
        threads=threads)    


@cli.command()
@click.argument('toph', type=click.Path(exists=True))
@click.argument('tax', type=click.Path(exists=False, writable=True))
@click.option('-c', '--min-conf', type=click.FLOAT, default=0.80,
              show_default=True, help="minimum confidence level")
@click.option('-n', '--name', type=click.STRING, help="custom name")
def gentax(toph, tax, min_conf, name):
    """Taxonomic classification from the allele top-hits"""

    gen.tax(toph, tax, min_conf=min_conf, name=name)
    