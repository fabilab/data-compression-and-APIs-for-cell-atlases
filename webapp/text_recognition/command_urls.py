# vim: fdm=indent
'''
author:     Fabio Zanini
date:       23/04/22
content:    Interpret text from Google TTS into a redirect URL.
'''
import numpy as np
import pandas as pd

from models import get_marker_genes



command_dict = {
    'expression_by_celltype': {
        'url_func': lambda sfx: f'/celltype/{sfx}',
    },
    'expression_unified_heatmap': {
        'url_func': lambda sfx: f'/heatmap_unified/{sfx}',
    },
    'gene_friends': {
        'url_func': 'TODO',  # TODO
    },
    'marker_genes': {
        'url_func': lambda sfx: '/celltype/'+get_marker_genes(sfx),
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

