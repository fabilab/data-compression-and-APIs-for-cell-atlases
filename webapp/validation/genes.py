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


def validate_correct_genestr(genestr):
    '''Validate gene names and correct misspellings if possible'''

    # TODO: check punctuation more accurately
    genes = genestr.replace('.', ',').replace(';', ',').split(',')

    # Murine genes are capitalized
    genes = [x.capitalize() for x in genes]

    # Validate
    genesv = []
    for gene in genes:
        if gene in genes_idx:
            genesv.append(gene)
            continue

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
            genesv.append(gene_closest)
            continue
        else:
            # At least one gene was not understood, ask for a written confirmation
            return None

    genestr = ','.join(genesv)
    return genestr
