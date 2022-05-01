# vim: fdm=indent
'''
author:     Fabio Zanini
date:       30/04/22
content:    Validate gene names et al.
'''
import pandas as pd

from models import (
    gene_order,
    )
genes_idx = gene_order.index
gene_matrix = (pd.Series(gene_order.index)
                 .str
                 .split('', expand=True)
                 .iloc[:, 1:-1]
                 .fillna('')
                 .values
                 .astype('U1'))
gene_maxlen = gene_matrix.shape[1]


def convert_numbers_in_gene_name(gene):
    '''Convert numbers in gene name into digits'''
    from text_recognition.assets import numbers

    # It's typical for gene names to end with a number (e.g. Car4), catch
    # those cases here
    for ntext, ndigit in numbers[::-1]:
        if gene.endswith(ntext):
            gene = gene[:-len(ntext)]+str(ndigit)
            return gene

    # NOTE: more complex cases, e.g. Col1a1, are more tricky to catch
    # and may require ML or whitelists


def validate_correct_gene(gene):
    '''Validate and correct misspellings for a single gene name'''

    # Murine genes are capitalized
    gene = gene.capitalize()

    if gene in genes_idx:
        return gene

    gene = convert_numbers_in_gene_name(gene)

    # Not found in whitelist, try to correct
    gene_array = list(gene)[:gene_maxlen]
    hamming = (gene_array != gene_matrix[:, :len(gene_array)]).sum(axis=1)

    # If the uncorrected gene is shorter (e.g. Col1) it can be a perfect match
    # for multiple whitelist genes (e.g. Col1a1, Col1a2), then ask for confirmation
    idx = (hamming == 0).nonzero()[0]
    # TODO: build mismatch scoring table based on natural English (e.g. p-t)
    if len(idx) == 0:
        idx = (hamming == 1).nonzero()[0]
    # If there is only one (partial) perfect match, take it. If multiple
    # perfect matches, ask. If no perfect and one imperfect match, take it.
    # If multiple imperfect or no imperfect, ask.
    if len(idx) == 1:
        gene_closest = genes_idx[idx[0]]
        print(gene, gene_closest)
        return gene_closest
    else:
        # At least one gene was not understood, ask for a written confirmation
        return None


def validate_correct_genestr(genestr):
    '''Validate gene names and correct misspellings if possible'''

    # TODO: check punctuation more accurately
    genes = genestr.replace('.', ',').replace(';', ',').split(',')

    # Validate
    genesv = []
    for gene in genes:
        genev = validate_correct_gene(gene)
        if genev is None:
            return None
        genesv.append(genev)

    genestr = ','.join(genesv)
    return genestr
