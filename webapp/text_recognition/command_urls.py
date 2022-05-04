# vim: fdm=indent
'''
author:     Fabio Zanini
date:       23/04/22
content:    Interpret text from Google TTS into a redirect URL.
'''
import numpy as np
import pandas as pd
from flask import url_for

from models import get_marker_genes, get_degs



command_dict = {
    'expression_by_celltype': {
        'url_func': lambda sfx: url_for('heatmap_by_celltype_genes', genestring=sfx),
    },
    'expression_unified_heatmap': {
        'url_func': lambda sfx: f'/heatmap_unified/{sfx}',
    },
    'gene_friends': {
        'url_func': 'TODO',  # TODO
    },
    'marker_genes': {
        'url_func': lambda sfx: url_for('heatmap_by_celltype_genes',
            genestring=get_marker_genes(sfx)),
    },
    'differentially_expressed_genes': {
        'url_func': lambda sfx: '/heatmap_differential/'+get_degs(sfx, kind='both'),
    },
    'upregulated_genes': {
        'url_func': lambda sfx: '/heatmap_differential/'+get_degs(sfx, kind='up'),
    },
    'downregulated_genes': {
        'url_func': lambda sfx: '/heatmap_differential/'+get_degs(sfx, kind='down'),
    },
    'list_cell_types': {
        'url_func': lambda sfx: url_for('list_celltypes_timepoint', timepoint=sfx),
    },
}


def get_command_response(text_dict):
    '''Format response to each given command'''
    if text_dict is None:
        # Default answer raises an alert in the frontend
        return {'outcome': 'fail'}

    suffix_corrected = text_dict['suffix_corrected']
    question = text_dict['question']
    category = text_dict['category']

    if suffix_corrected is None:
        return {
            'outcome': 'question',
            'question': question,
            question: text_dict['suffix'],
            'url_prefix': command_dict[category]['url_func'](''),
        }

    url = command_dict[category]['url_func'](suffix_corrected)
    return {
        'outcome': 'success',
        'url': url,
    }

