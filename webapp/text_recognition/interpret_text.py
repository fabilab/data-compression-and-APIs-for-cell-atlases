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
    )
from validation import (
    validate_correct_genestr,
    validate_correct_celltypestr,
    validate_correct_celltypedatasettimepoint,
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
            'exp of',
            'e of',
            '!e',
        ],
        'suffix_type': 'genestring',
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
            'prog of',
            'p of',
            '!p',
        ],
        'suffix_type': 'genestring',
    },
    'gene_friends': {
        'prefixes': [
            'what are gene friends of',
            'gene friends of',
            'correlates of',
            'show the gene friends of',
            'show gene friends of',
            'show correlates of',
            'friends of',
            'f of',
            '!f',
        ],
        'suffix_type': 'genestring',
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
            'show the markers for','hyperoxia'
            'show the marker genes of',
            'show the marker genes for',
            '!m',
        ],
        'suffix_type': 'celltypestring',
    },
    'differentially_expressed_genes': {
        'prefixes': [
            'differentially expressed genes in',
            'degs in',
            '!d',
        ],
        'suffix_type': 'celltype_dataset_timepoint_string',
    },
    'upregulated_genes': {
        'prefixes': [
            'genes upregulated in',
            'upregulated genes in',
            'upregulated in',
            'up in',
            '!u',
        ],
        'suffix_type': 'celltype_dataset_timepoint_string',
    },
    'downregulated_genes': {
        'prefixes': [
            'genes downregulated in',
            'downregulated genes in',
            'downregulated in',
            'down in',
            '!d',
        ],
        'suffix_type': 'celltype_dataset_timepoint_string',
    }
}
phrase_dict_inv = {}
for key, val in phrase_dict.items():
    for prefix in val['prefixes']:
        phrase_dict_inv[prefix] = key


def infer_command_from_text(text_raw):
    '''Figure out category of command'''
    text = text_raw.lower()

    for prefix, category in phrase_dict_inv.items():
        if not text.startswith(prefix):
            continue

        # Figure out type of suffix
        suffix_type = phrase_dict[category]['suffix_type']

        # Cut prefix
        if suffix_type == 'celltype_dataset_timepoint_string':
            suffix = text_raw[len(prefix):]
        else:
            suffix = text[len(prefix):]
            # Remove whitespace from suffix
            suffix = suffix.replace(' ', '')

        # Cut punctuation at the end of the command
        suffix = suffix.rstrip('?!.')

        return {
            'prefix': prefix,
            'suffix': suffix,
            'category': category,
            }
    return None


def interpret_text(text):
    inferred_dict = infer_command_from_text(text)
    if inferred_dict is None:
        return None

    prefix = inferred_dict['prefix']
    suffix = inferred_dict['suffix']
    category = inferred_dict['category']

    suffix_type = phrase_dict[category]['suffix_type']
    if suffix_type == 'genestring':
        suffix_corrected = validate_correct_genestr(suffix)
        question = 'gene_string'
    elif suffix_type == 'celltypestring':
        suffix_corrected = validate_correct_celltypestr(suffix)
        question = 'celltype_string'
    elif suffix_type == 'celltype_dataset_timepoint_string':
        suffix_corrected = validate_correct_celltypedatasettimepoint(suffix)
        question = 'celltype_dataset_timepoint_string'
    else:
        raise ValueError('Category not implemented')

    return {
        'prefix': prefix,
        'suffix': suffix,
        'suffix_corrected': suffix_corrected,
        'category': category,
        'question': question,
    }
