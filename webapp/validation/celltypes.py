# vim: fdm=indent
'''
author:     Fabio Zanini
date:       27/04/22
content:    Some quirks of the celltypes in the data.
'''
import numpy as np


celltype_tuples = [
    ('Adventitial fibroblast', 'Adventitial FB'),
    ('Early adventitial fibroblast', 'Early adventitial FB'),
    ('Fibroblast precursor', 'FB precursor'),
    ('Early alveolar fibroblast', 'Early alveolar FB'),
    ('Alveolar fibroblast', 'Alveolar FB'),
    ('Proliferating fibroblast', 'Proliferating FB'),
    ('Proliferating myofibroblast', 'Proliferating MyoF'),
    ('Myofibroblast', 'MyoF'),
    ('Myofibroblast and smooth muscle precursor', 'MyoF/ASM precursor'),
    ('Early airway smooth muscle', 'Early ASM'),
    ('Airway smooth muscle', 'ASM'),
    ('Vascular smooth muscle', 'VSM'),
    'Pericyte',
    'Proliferating pericyte',
    'Striated muscle',
    ('Lymphatic EC', 'Lymphatic'),
    ('Arterial EC II', 'Arterial I'),
    ('Arterial EC I', 'Arterial II'),
    ('Venous EC', 'Venous'),
    ('Nonproliferative embryonic EC', 'Embryonic cap'),
    ('Proliferative EC', 'Proliferative cap'),
    ('Early Car4- capillaries', 'Early gCap'),
    ('Late Car4- capillaries', 'Late gCap'),
    ('Car4+ capillaries', 'Aerocyte'),
    'B cell',
    'NK cell',
    'T cell',
    'IL cell',
    'DC I',
    'DC II',
    'DC III',
    'Mac I',
    'Mac II',
    'Mac III',
    'Mac IV',
    ('Mac V', 'Monocyte'),
    ('basophil', 'Basophil'),
    ('mast cell', 'Mast cell'),
    ('neutrophil', 'Neutrophil'),
    'Alveolar type I',
    'Alveolar type II',
]
celltype_dict = {}
for ct in celltype_tuples:
    if isinstance(ct, str):
        celltype_dict[ct] = ct
    else:
        celltype_dict[ct[0]] = ct[1]


def adjust_celltypes(celltypes_raw):
    # TODO: reorder even when it's stratified by dataset/timepoint
    for ct in celltypes_raw:
        if '_ACZ' in ct:
            return celltypes_raw, np.arange(len(celltypes_raw))

    celltypes_raw = list(celltypes_raw)
    ct_adj = []
    idx = []
    for ct in celltype_tuples:
        if isinstance(ct, str):
            ct1, ct2 = ct, ct
        else:
            ct1, ct2 = ct
        if ct1 in celltypes_raw:
            idx.append(celltypes_raw.index(ct1))
            ct_adj.append(ct2)
    return np.array(ct_adj), np.array(idx)


def rename_celltypes(celltypes_raw):
    '''Rename celltypes according to the standard table above, no order change'''
    return [celltype_dict[x] for x in celltypes_raw]


def validate_correct_celltypestr(celltypestr):
    '''Validate cell type names and correct misspellings if possible'''
    # TODO: check punctuation more accurately
    celltypes = celltypestr.replace('.', ',').replace(';', ',').split(',')

    # Capitalization is not great, lowercase check
    celltypes = [x.lower() for x in celltypes]
    celltypesd = {x.lower(): x for x in celltypes_all}

    # Validate
    celltypesv = []
    for celltype in celltypes:
        if celltype in celltypesd:
            celltypesv.append(celltypesd[celltype])
            continue

        # Cut plural if it is found
        if celltype.endswith('s'):
            celltype = celltype[:-1]
            if celltype in celltypesd:
                celltypesv.append(celltypesd[celltype])
                continue

        # TODO: implement more spelling correction
        return None

    celltypestr = ','.join(celltypesv)
    return celltypestr


def validate_correct_celltypedatasettimepoint(suffix):
    '''Validate celltype + dataset + timepoint'''
    #TODO
    return suffix

