# vim: fdm=indent
'''
author:     Fabio Zanini
date:       23/04/22
content:    Interpret text from Google TTS into a redirect URL.
'''
import numpy as np
import pandas as pd

from models import (
    gene_order,
    celltypes as celltypes_all,
    get_marker_genes,
    )
from validation import (
    validate_correct_genestr,
    validate_correct_celltypestr,
    )


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
        'suffix_type': 'genestring',
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
        'suffix_type': 'genestring',
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
        'suffix_type': 'genestring',
        'url_func': 'TODO',  # TODO
    },
    'marker_genes': {
        'prefixes': [
            'what are the markers of',
            'what are markers of',
            'what are marker genes of',
            'markers for',
            'markers of',
            'marker genes of',
            'marker genes for',
            'show markers of',
            'show markers for',
            'show marker genes of',
            'show marker genes for',
            'show the markers of',
            'show the markers for',
            'show the marker genes of',
            'show the marker genes for',
            'upregulated in',
            'upregulated genes in',
            'genes upregulated in',
        ],
        'suffix_type': 'celltypestring',
        'url_func': lambda sfx: '/celltype/'+get_marker_genes(sfx)
    }
}
phrase_dict_inv = {}
for key, val in phrase_dict.items():
    for prefix in val['prefixes']:
        phrase_dict_inv[prefix] = key


def text_to_url(text):
    '''Interpret natural language into categories

    Of course ML would be useful here, but we just do a few rules for now.
    '''
    print(type(text))
    text = text.lower()

    for prefix, category in phrase_dict_inv.items():
        if text.startswith(prefix):
            suffix = text[len(prefix):].replace(' ', '')

            if phrase_dict[category]['suffix_type'] == 'genestring':
                suffix_corrected = validate_correct_genestr(suffix)
                question = 'gene_string'
            elif phrase_dict[category]['suffix_type'] == 'celltypestring':
                suffix_corrected = validate_correct_celltypestr(suffix)
                question = 'celltype_string'
            else:
                raise ValueError('Category not implemented')

            if suffix_corrected is None:
                return {
                    'outcome': 'question',
                    'question': question,
                    question: suffix,
                    'url_prefix': phrase_dict[category]['url_func'](''),
                }

            url = phrase_dict[category]['url_func'](suffix_corrected)
            return {
                'outcome': 'success',
                'url': url,
            }

    # Default answer raises an alert in the frontend
    return {'outcome': 'fail'}
