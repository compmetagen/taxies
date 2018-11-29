import sys
import os
import csv
from io import StringIO
import subprocess

import click

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# Suppress experimental warnings in SearchIO
import warnings
from Bio import BiopythonExperimentalWarning
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio import SearchIO

from . import db
from . import utils

  

def mark(genome_fn, db_dir, mark_fn, max_evalue=0.001, threads=1):

    def write(nucleotide_sequences_fn, hmm_target_fn, mark_fn):
        nucleotide_sequences = SeqIO.to_dict( \
            SeqIO.parse(nucleotide_sequences_fn, "fasta"))
        
        n_mark = 0
        mark_handle = open(mark_fn, 'w')
        for qres in SearchIO.parse(hmm_target_fn, 'hmmer3-tab'):

            hit = qres.hits[0]
            target_best, evalue_best = hit.id, hit.evalue
            for hit in qres.hits[1:]:
                if hit.evalue < evalue_best:
                    target_best = hit.id
                    evalue_best = hit.evalue

            rec = nucleotide_sequences[target_best]
            rec.id = qres.id
            rec.description = ""
            SeqIO.write(rec, mark_handle, "fasta")
            n_mark += 1

        mark_handle.close()

        return n_mark


    mark_name = os.path.basename(os.path.splitext(mark_fn)[0])
    mark_dir = os.path.dirname(mark_fn)
    mark_prefix = os.path.join(mark_dir, mark_name)
    profiles_fn = os.path.join(db_dir, db.PROFILES_HMM_FN)
    pnames_fn = os.path.join(db_dir, db.PROTEIN_NAMES_FN)
    pnames = utils.read_protein_names(pnames_fn)

    prodigal_cmd = ['prodigal',
                    '-i', genome_fn,
                    '-o', mark_prefix + '.gene.coords.gbk',
                    '-a', mark_prefix + '.protein.translations.faa',
                    '-d', mark_prefix + '.nucleotide.sequences.fna']
    prodigal_proc = subprocess.run(prodigal_cmd, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)

    hmmsearch_cmd = ['hmmsearch',
                     '--tblout', mark_prefix + '.hmm.target.txt',
                     '-E', str(max_evalue),                    
                     '--cpu', str(threads),
                     profiles_fn,
                     mark_prefix + '.protein.translations.faa']
    hmmsearch_proc = subprocess.run(hmmsearch_cmd, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)

    n_mark = write(mark_prefix + '.nucleotide.sequences.fna', 
                   mark_prefix + '.hmm.target.txt',
                   mark_fn)

    os.remove(mark_prefix + '.gene.coords.gbk')
    os.remove(mark_prefix + '.protein.translations.faa')
    os.remove(mark_prefix + '.nucleotide.sequences.fna')
    os.remove(mark_prefix + '.hmm.target.txt')

    click.echo("{:d}/{:d} markers found".format(n_mark, len(pnames)))


def classify(mark_fn, db_dir, cls_fn, custom_pnames_fn=None, min_supp=0.80,
    threads=1):

    if custom_pnames_fn is None:
        pnames_fn = os.path.join(db_dir, db.PROTEIN_NAMES_FN)
        pnames = utils.read_protein_names(pnames_fn)
    else:
        pnames = utils.read_protein_names(custom_pnames_fn)

    mark_rec_dict = SeqIO.to_dict(SeqIO.parse(mark_fn, "fasta"))

    cls_handle = open(cls_fn, 'w')
    cls_writer = csv.writer(cls_handle, delimiter='\t', lineterminator='\n')

    with click.progressbar(pnames, label="Classification") as bar:
        for pname in bar:
            
            alleles_fn = os.path.join(db_dir, db.DATA_DIR, "{}{}".\
                format(pname, db.ALLELES_SUFFIX))

            tax, status = "", ""

            try:
                marker_str = mark_rec_dict[pname].format("fasta")
            except KeyError:
                pass
            else:
                status = "*"

                sintax_cmd = ["vsearch",
                              "--sintax", "-",
                              "--db", alleles_fn,
                              "--tabbedout", "-",
                              "--sintax_cutoff", str(min_supp),
                              "--strand", "both",
                              "--threads", str(threads)]

                sintax_proc = subprocess.run(sintax_cmd, input=marker_str,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                    universal_newlines=True)

                stdout_handle = StringIO(sintax_proc.stdout)
                stdout_reader = csv.reader(stdout_handle, delimiter='\t')

                try:
                    tax = stdout_reader.__next__()[-1]
                except StopIteration:
                    pass

            cls_writer.writerow([pname, tax, status])

        cls_handle.close()


def toph(mark_fn, db_dir, toph_fn, custom_pnames_fn=None, min_ident=0.60,
    threads=1):

    if custom_pnames_fn is None:
        pnames_fn = os.path.join(db_dir, db.PROTEIN_NAMES_FN)
        pnames = utils.read_protein_names(pnames_fn)
    else:
        pnames = utils.read_protein_names(custom_pnames_fn)
    
    mark_rec_dict = SeqIO.to_dict(SeqIO.parse(mark_fn, "fasta"))
    click.echo("{:d} markers will be used". \
        format(len(set(pnames) & set(mark_rec_dict.keys()))))

    toph_handle = open(toph_fn, 'w')
    toph_writer = csv.writer(toph_handle, delimiter='\t',
        lineterminator='\n')

    with click.progressbar(pnames, label="") as bar:
        for pname in bar:
            
            alleles_fn = os.path.join(db_dir, db.DATA_DIR, "{}{}".\
                format(pname, db.ALLELES_SUFFIX))

            tax, status, allele_id, ident, qcov = "", "", "", "", ""
            
            try:
                marker_str = mark_rec_dict[pname].format("fasta")
            except KeyError:
                pass
            else:
                status = "*"

                vsearch_cmd = ["vsearch",
                               "--usearch_global", "-",
                               "--db", alleles_fn,
                               "--userout", "-",
                               "--id", str(min_ident),
                               "--userfields", "query+target+id+qcov",
                               "--strand", "both",
                               "--threads", str(threads)]

                vsearch_proc = subprocess.run(vsearch_cmd, input=marker_str,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                    universal_newlines=True)

                stdout_handle = StringIO(vsearch_proc.stdout)
                stdout_reader = csv.reader(stdout_handle, delimiter='\t')

                try:
                    Q, T_tax, ident, qcov = stdout_reader.__next__()
                    T, tax = T_tax.split(";tax=")
                    allele_id = T.rsplit("_", 1)[1]
                except StopIteration:
                    pass

            toph_writer.writerow([pname, tax, status, allele_id,
                ident, qcov])

        toph_handle.close()


def tax(input_fn, tax_fn, min_conf=0.8, name=None):
    taxtable_df = utils.read_taxtable(input_fn)
    cladogram = utils.taxtable2cladogram(taxtable_df)
    tot_count = sum([n.count for n in cladogram.seed_node.child_node_iter()])
    
    nodes_path = []
    node = cladogram.seed_node
    while node.num_child_nodes() > 0:
        count_max, child_node_max = 0, None
        for child_node in node.child_node_iter():
            if child_node.count > count_max:
                count_max = child_node.count
                child_node_max = child_node
        nodes_path.append(child_node_max)
        node = child_node_max

    tax_conf = ";".join( \
        ["{}({:.2f})".format(n.taxon.label, n.count/tot_count) \
        for n in nodes_path])

    tax = ";".join( \
        ["{}".format(n.taxon.label, n.count/tot_count) \
        for n in nodes_path if (n.count/tot_count) > min_conf])

    perc_hits = 100*tot_count/cladogram.seed_node.count

    if name is None:
        name = os.path.basename(input_fn)

    with open(tax_fn, 'w') as tax_handle:
        tax_handle.write("{}\t{}\t{}\t{:.2f}\n". \
            format(name, tax, tax_conf, perc_hits))
