import os

import click

from taxies import asm


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
def asmmark(genome, db, mark, max_evalue, threads):
    """Extract markers from assembled genomes"""
    asm.mark(genome, db, mark, max_evalue=max_evalue, threads=threads)


@cli.command()
@click.argument('mark', type=click.Path(exists=True))
@click.argument('db', type=click.Path(exists=True, file_okay=False))
@click.argument('tax', type=click.Path(exists=False, writable=True))
@click.option('-P', '--pnames', type=click.Path(exists=True),
    help="custom protein names file")
@click.option('-p', '--toph', type=click.Path(exists=False, writable=True),
    help="write top hits in tab-delimited format")
@click.option('-c', '--min-conf', type=click.FLOAT, default=0.80,
              show_default=True, help="minimum confidence level")
@click.option('-d', '--min-ident', type=click.FLOAT, default=0.60,
              show_default=True, help="minimum top hit identity")
@click.option('-t', '--threads', type=click.INT, default=1,
              show_default=True, help="number of threads")

def asmtax(mark, db, tax, pnames, toph, min_conf, min_ident, threads):
    """Taxonomic classification from markers"""

    tax_prefix = os.path.splitext(tax)[0]
    toph_tmp = tax_prefix + ".tophits.txt"

    asm.toph(mark, db, toph_tmp, custom_pnames_fn=pnames, 
        min_ident=min_ident, threads=threads)
    asm.lca(toph_tmp, tax, min_conf=min_conf)
    
    if toph is not None:
        os.rename(toph_tmp, toph)
    else:
        os.remove(toph_tmp)
    
    