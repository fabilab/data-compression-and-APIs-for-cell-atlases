# vim: fdm=indent
'''
author:     Fabio Zanini
date:       16/05/22
content:    Compress mouse lemur atlas.
'''
import os
import sys
import numpy as np
import pandas as pd
import h5py
import anndata



if __name__ == '__main__':

    # Load data
    fn_atlas = '../data/tabula_sapiens/TS_Lung.h5ad'
    adata = anndata.read_h5ad(fn_atlas)

    # Find cell type column
    column = 'cell_ontology_class'
    # NOTE: there is a free annotation column that appears to be exactly the
    # same information yet again, and very much not "free"... ask Bob?

    # Set cell type order
    #celltypes = adata.obs[column].unique()
    celltypes = [
        'Adventitial FB',
        'Alveolar FB',
        'Mesothelial',
        'MyoF',
        'ASM',
        'VSM',
        'Pericyte',
        'Venous',
        'Cap',
        'Aerocyte',
        'Art',
        'Lymphatic',
        'Bronchial EC',
        'Plasmablast',
        'B cell',
        'NK cell',
        'pDC',
        'cDC',
        'Alv/interst mac',
        'Monocyte',
        'Neutrophil',
        'Basophil',
        'Alveolar type I',
        'Alveolar type II',
        'Basal',
        'Ciliated',
        'Club',
        'Goblet',
        'Ionocyte',
        'Serous',
        'Mucous',
    ]

    # Merge some annotations together
    annotation_merges = {
        'Adventitial FB': ['adventitial cell'],
        'Alveolar FB': ['alveolar fibroblast'],
        'Mesothelial': ['mesothelial cell'],
        'MyoF': ['myofibroblast cell'],
        'ASM': ['bronchial smooth muscle cell'],
        'VSM': ['vascular associated smooth muscle cell'],
        'Pericyte': ['pericyte cell'],
        'Venous': ['vein endothelial cell'],
        'Cap': [
            'capillary endothelial cell',
            'lung microvascular endothelial cell',
        ],
        'Aerocyte': ['capillary aerocyte'],
        'Art': ['endothelial cell of artery'],
        'Lymphatic': ['endothelial cell of lymphatic vessel'],
        'Bronchial EC': ['bronchial vessel endothelial cell'],
        'Ionocyte': ['pulmonary ionocyte'],
        'Plasmablast': ['plasma cell'],
        'NK cell': ['nk cell'],
        'IL cell': ['innate lymphoid cell'],
        'B cell': ['b cell'],
        'T cell': [
            'T cell',
            'cd8-positive, alpha-beta t cell',
            'cd4-positive, alpha-beta t cell',
            'cd8-positive alpha-beta t cell',
            'cd4-positive alpha-beta t cell',
        ],
        'Alveolar type I': ['type i pneumocyte'],
        'Alveolar type II': ['type ii pneumocyte'],
        'Goblet': ['respiratory goblet cell'],
        'Ciliated': ['lung ciliated cell'],
        'Club': ['club cell'],
        'Basal': ['basal cell'],
        'Serous': ['serous cell of epithelium of bronchus'],
        'Mucous': ['respiratory mucous cell'],
        'cDC': ['dendritic cell'],
        'pDC': ['plasmacytoid dendritic cell'],
        'Alv/interst mac': ['macrophage'],
        'Monocyte': [
            'classical monocyte',
            'non-classical monocyte',
            'intermediate monocyte',
            ],
        'Neutrophil': ['neutrophil'],
        'Basophil': ['basophil'],
    }
    adata.obs[column+'_compressed'] = ''
    for ct, csts in annotation_merges.items():
        adata.obs.loc[adata.obs[column].isin(csts), column+'_compressed'] = ct

    # Average, proportion expressing, number of cells
    genes = adata.var_names
    avg_exp = pd.DataFrame(
            np.zeros((len(genes), len(celltypes)), np.float32),
            index=genes,
            columns=celltypes)
    frac_exp = avg_exp.copy()
    ncells = pd.Series(np.zeros(len(celltypes), np.int64), index=celltypes)
    for ct in celltypes:
        # Avg gene expression (compute cpm after averaging)
        avg_exp[ct] = np.asarray(adata[adata.obs[column+'_compressed'] == ct].X.mean(axis=0))[0]
        avg_exp[ct] *= 1e6 / avg_exp[ct].sum()
        # Proportion expressing
        frac_exp[ct] = np.asarray((adata[adata.obs[column+'_compressed'] == ct].X > 0).mean(axis=0))[0]
        # Number of cells
        ncells[ct] = (adata.obs[column+'_compressed'] == ct).sum()

    # Store to file
    fn_out = '../webapp/static/scData/human_condensed_lung_atlas_in_cpm.h5'
    with pd.HDFStore(fn_out) as h5_data:
        h5_data.put('celltype/gene_expression_average', avg_exp.T)
        h5_data.put('celltype/gene_proportion_expression', frac_exp.T)
        h5_data.put('celltype/cell_count', ncells)

