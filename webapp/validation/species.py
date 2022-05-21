# vim: fdm=indent
'''
author:     Fabio Zanini
date:       22/05/22
content:    Validate and correct questions related to comparisons between
            species.
'''
from .genes import validate_correct_genestr


def validate_correct_species(text):
    '''Validate species names'''
    text = text.lower().rstrip('s')
    if text in ('human', 'mouse lemur', 'lemur'):
        return text

    if text == 'mouse':
        return text

    if text == 'monkey':
        return 'mouse lemur'

    return None


def validate_correct_species_genestr(species, text):
    '''Validate comparison of expression between species'''
    # The suffix starts with a genestr, then comes something like 'vs <species>'
    # the target species has already been excised
    sep = ' vs '
    if sep not in text:
        return None

    idx = text.index(sep)
    genestr = text[:idx]
    genestr = validate_correct_genestr(genestr)

    # TODO: expand to more complex syntax, e.g. specific time points
    species_baseline = text[idx+len(sep):]
    species_baseline = validate_correct_species(species_baseline)

    return {
        'species': species,
        'species_baseline': species_baseline,
        'genes': genestr,
    }
