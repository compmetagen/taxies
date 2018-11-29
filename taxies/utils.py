import os
import re
import csv

import pandas as pd
import dendropy


def read_protein_names(input_fn):
        
    protein_names = set()
    with open(input_fn, 'r') as input_handle:
        for row in csv.reader(input_handle, delimiter='\t'):
            if len(row) > 0:
                protein_names.add(row[0])
                
    return list(protein_names)


class TaxReader:
    def __init__(self, handle):
        self.__handle = handle
        self.__reader = csv.reader(handle, delimiter='\t')

    def __iter__(self):
        return self

    def __parse_tax(self, s):
        """Parse a taxonomy string and returns a list.
        """
        tax = []
        for elem in re.split(';|,', s):
            elem = re.sub(r"^\S__|^D_\d+__|^\S:", "", elem.strip())
            elem = re.sub(r"__SPACE__", " ", elem.strip())
            if elem == "":
                break
            tax.append(elem)
        return tax

    def __next__(self):
        row = self.__reader.__next__()
        return row[0], self.__parse_tax(row[1])


def read_taxtable(input_fn):
    """
    """

    tax_dict = dict()
    with open(input_fn, 'rU') as input_handle:
        taxreader = TaxReader(input_handle)
        for seqid, tax in taxreader:
            tax_dict[seqid] = tax

    return pd.DataFrame.from_dict(tax_dict, orient='index')


def taxtable2cladogram(df):
    def find_child_nodes_with_taxon_label(node, label):
        return [child for child in node.child_node_iter() \
            if child.taxon.label == label]

    tree = dendropy.Tree()   
    tree.is_rooted = True
    taxon_namespace = dendropy.TaxonNamespace()
   
    seed_taxon = dendropy.Taxon("Root")
    taxon_namespace.add_taxon(seed_taxon)
    tree.seed_node.taxon = seed_taxon

    tree.seed_node.count = 0
    tree.seed_node.rank = 0
    tree.max_rank = 0

    for index, row in df.iterrows():
        tree.seed_node.count += 1
        node_curr = tree.seed_node
        
        for i, taxon in enumerate(row):
            if taxon is None:
                break
            else:
                nodes_taxon = \
                    find_child_nodes_with_taxon_label(node_curr, taxon)
                
                if len(nodes_taxon) == 0:
                    taxon_taxon = dendropy.Taxon(taxon)
                    taxon_namespace.add_taxon(taxon_taxon)
                    child = dendropy.Node(taxon=taxon_taxon)
                    child.count = 1
                    child.rank = i+1
                    if child.rank > tree.max_rank:
                        tree.max_rank = child.rank
                    node_curr.add_child(child)
                else:
                    child = nodes_taxon[0]
                    child.count += 1

                node_curr = child

    tree.taxon_namespace = taxon_namespace

    return tree
