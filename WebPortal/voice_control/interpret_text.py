# vim: fdm=indent
'''
author:     Fabio Zanini
date:       23/04/22
content:    Interpret text from Google TTS into a redirect URL.
'''
import numpy as np
import pandas as pd
# FIXME: Absolute import??
from helper import gene_order
genes_idx = gene_order.index
gene_matrix = (pd.Series(gene_order.index)
                 .str
                 .split('', expand=True)
                 .iloc[:, 1:-1]
                 .fillna('')
                 .values
                 .astype('U1'))
gene_maxlen = gene_matrix.shape[1]


phrase_dict = {
    'expression_by_celltype': {
        'prefixes': [
            'what is the expression of',
            'what\'s the expression of',
            'expression of',
            'gene expression of',
            'show expression of',
            'show gene expression of',
            'show the expression of',
            'show the gene expression of',
        ],
        'url_func': lambda sfx: f'/celltype/{sfx}'
    },
    'expression_unified_heatmap': {
        'prefixes': [
            'what is the developmental expression of',
            'what\'s the developmental expression of',
            'progression of',
            'developmental expression of',
            'developmental progression of',
            'expression across development of',
            'show developmental progression of',
            'show developmental expression of',
            'show the developmental expression of',
            'show the developmental gene expression of',
        ],
        'url_func': lambda sfx: f'/heatmap_unified/{sfx}'
    },
    'gene_friends': {
        'prefixes': [
            'what are gene friends of',
            'gene friends of',
            'correlates of',
            'show the gene friends of',
            'show gene friends of',
            'show correlates of',
        ],
        'url_func': 'TODO',  # TODO
    },
}
phrase_dict_inv = {}
for key, val in phrase_dict.items():
    for prefix in val['prefixes']:
        phrase_dict_inv[prefix] = key



def validate_correct_genestr(genestr):
    '''Validate gene names and correct misspellings if possible'''

    # TODO: check punctuation more accurately
    genes = genestr.replace('.', ',').replace(';', ',').split(',')

    # Murine genes are capitalized
    genes = [x.capitalize() for x in genes]

    # Validate and check
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

    genestr = ','.join(genes)
    return genestr


def text_to_url(text):
    '''Interpret natural language into categories

    Of course ML would be useful here, but we just do a few rules for now.
    '''
    print(type(text))
    text = text.lower()

    for prefix, category in phrase_dict_inv.items():
        if text.startswith(prefix):
            genestr = text[len(prefix):].replace(' ', '')

            genestr_corrected = validate_correct_genestr(genestr)

            if genestr_corrected is None:
                return {
                    'outcome': 'question',
                    'question': 'gene_string',
                    'gene_string': genestr,
                    'url_prefix': phrase_dict[category]['url_func'](''),
                }

            url = phrase_dict[category]['url_func'](genestr_corrected)
            return {
                'outcome': 'success',
                'url': url,
            }

    # Default answer raises an alert in the frontend
    return {'outcome': 'fail'}
