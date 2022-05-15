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
from validation.genes import validate_correct_genestr
from validation.celltypes import (
    validate_correct_celltypestr,
    validate_correct_celltypedatasettimepoint,
    )
from validation.timepoints import validate_correct_timepoint


phrase_dict = {
    'expression_by_celltype': {
        'prefixes': [
            'what is the expression of',
            'what\'s the expression of',
            'what is the gene expression of',
            'what\'s the gene expression of',
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
            'what is the developmental progression of',
            'what\'s the developmental progression of',
            'what is the gene progression of',
            'what\'s the gene progression of',
            'progression of',
            'developmental expression of',
            'developmental progression of',
            'expression across development of',
            'show developmental progression of',
            'show developmental expression of',
            'show the developmental expression of',
            'show the developmental gene expression of',
            'show the progression of',
            'show progression of',
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
            'show the markers for',
            'show the marker genes of',
            'show the marker genes for',
            '!m',
        ],
        'suffix_type': 'celltypestring',
    },
    'differentially_expressed_genes': {
        'prefixes': [
            'differentially expressed genes',
            'degs',
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
    },
    'list_cell_types': {
        'prefixes': [
            'list cell types',
            'list the cell types',
            'list cell type ',
            'list the cell type '
            'show cell types',
            'what are the cell types',
            'what cell types are there',
            'what cell types are present',
            'cell types',
            '!ct',
        ],
        'suffix_type': 'timepoint',
    },
    'celltype_abundance': {
        'prefixes': [
            'show cell type abundance at',
            'cell type abundance at',
            'abundance of cell types at',
        ],
        'suffix_type': 'timepoint',
    },
}
phrase_dict_inv = {}
for key, val in phrase_dict.items():
    for prefix in val['prefixes']:
        phrase_dict_inv[prefix] = key


def infer_command_from_text(text_raw):
    '''Figure out category of command'''
    # Cut punctuation at the end of the command
    text = text_raw.rstrip('?!.')

    # Prefixes are checked case-insensitive
    text_lower = text.strip().lower()

    for prefix, category in phrase_dict_inv.items():
        if not text_lower.startswith(prefix):
            continue

        # Figure out type of suffix
        suffix_type = phrase_dict[category]['suffix_type']

        # Cut prefix... there are two common cases for how to treat the suffix
        cats_keep_whitespace = (
            'celltype_dataset_timepoint_string',
            'timepoint',
            'celltypestring',
            )
        if suffix_type in cats_keep_whitespace:
            suffix = text[len(prefix):]
        else:
            suffix = text_lower[len(prefix):]

        # Extract species from suffix
        suffix, species = excise_species_from_suffix(suffix)

        if suffix_type not in cats_keep_whitespace:
            # Remove whitespace from suffix
            suffix = suffix.replace(' ', '')

        return {
            'prefix': prefix,
            'suffix': suffix,
            'category': category,
            'species': species,
            }
    return None


def excise_phrase(text, phrase):
    '''Excise a phrase from a text'''
    if phrase not in text:
        return text
    if text.endswith(phrase):
        # Assume space before
        return text[:-len(phrase)-1]
    # Assume space after
    idx = text.find(phrase+' ')
    return text[:idx]+text[idx+len(phrase)+1:]


def excise_species_from_suffix(suffix):
    '''Excise species from suffix if found'''
    phrases = {
        'human': ['in humans', 'in human'],
        'mouse': ['in mouse', 'in mice'],
    }
    for species, phrases_species in phrases.items():
        for phrase in phrases_species:
            if phrase in suffix:
                suffix = excise_phrase(suffix, phrase)
                return suffix, species
    return suffix, 'mouse'


def interpret_text(text):
    '''Interpret natural language text as command'''
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
        suffix_corrected = validate_correct_celltypedatasettimepoint(
                suffix)
        question = 'celltype_dataset_timepoint_string'
    elif suffix_type == 'timepoint':
        suffix_corrected = validate_correct_timepoint(suffix)
        question = 'timepoint'
    else:
        raise ValueError('Category not implemented')

    return {
        'prefix': prefix,
        'suffix': suffix,
        'suffix_corrected': suffix_corrected,
        'category': category,
        'question': question,
        'species': inferred_dict['species'],
    }
