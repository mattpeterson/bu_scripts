#!/usr/bin/env python
"""
get_genes_for_terms.py

Given a set of enriched genes, a mapping between genes and terms (pathways, GO
terms, etc.), returns the set of genes matching terms

Usage:
    get_genes_for_terms [gene list] [mapping] [comma-separated term list]

TODO:  This should get merged into a FunctionalMapping class or something, as
part of transnet.
"""
import sys
from collections import defaultdict

def read_gene_set(gene_file):
    """
    Reads in a set of gene names from a file
    """
    genes = [g.rstrip() for g in gene_file]
    return genes

def read_mapping(mapping_file):
    """
    Reads in a mapping from gene to term.  Expected to be flattened - one term
    per line.

    :param mapping_handle: handle to be read from
    :type mapping_handle: file
    """
    mapping = defaultdict(list)

    for line in mapping_file:
        tokens = line.rstrip("\r\n").split()
        mapping[tokens[0]].append(tokens[1])

    return mapping

def get_matching_genes(gene_list, mapping, term):
    """
    Gets the set of genes from the gene set.
    """
    matched = []

    for g in gene_list:
        terms = mapping[g]
        if term in terms:
            matched.append(g)

    return matched

if __name__ == "__main__":
    with open(sys.argv[1]) as geneset_handle:
        genes = read_gene_set(geneset_handle)

    with open(sys.argv[2]) as mapping_handle:
        mapping = read_mapping(mapping_handle)


    terms_to_match = sys.argv[3].split(",")

    for term in terms_to_match:
        matching_genes = get_matching_genes(genes, mapping, term)
        print("%s\t%s" % (term, ",".join(matching_genes)))

