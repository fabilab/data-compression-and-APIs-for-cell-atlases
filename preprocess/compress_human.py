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


data_fdn = '../webapp/static/scData/'


if __name__ == '__main__':

    # Load data
    print('Load single cell data')
    fn_atlas = '../data/tabula_sapiens/TS_Lung.h5ad'
    adata = anndata.read_h5ad(fn_atlas)

    # NOTE: the human data is in some weird normalization between 0 and 10,
    # use adata.raw.X for computations to avoid log trasformations and whatnot

    print('Compress')
    # Find cell type column
    column = 'cell_ontology_class'
    # NOTE: there is a free annotation column that appears to be exactly the
    # same information yet again, and very much not "free"... ask Bob?

    # Set cell type order
    #celltypes = adata.obs[column].unique()
    celltypes = [
        'Mesothelial',
        'Adventitial FB',
        'Alveolar FB',
        'MyoF',
        'ASM',
        'VSM',
        'Pericyte',
        'Venous',
        'Cap',
        'Aerocyte',
        'Art II',
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
        'Art II': ['endothelial cell of artery'],
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
        # Avg gene expression (use adata.raw.X, compute cpm after averaging)
        avg_exp[ct] = np.asarray(adata[adata.obs[column+'_compressed'] == ct].raw.X.mean(axis=0))[0]
        avg_exp[ct] *= 1e6 / avg_exp[ct].sum()
        # Proportion expressing (use adata.raw.X)
        frac_exp[ct] = np.asarray((adata[adata.obs[column+'_compressed'] == ct].raw.X > 0).mean(axis=0))[0]
        # Number of cells
        ncells[ct] = (adata.obs[column+'_compressed'] == ct).sum()

    # Store to file
    print('Store compressed atlas')
    fn_out = data_fdn + 'human_condensed_lung_atlas_in_cpm.h5'
    with pd.HDFStore(fn_out) as h5_data:
        h5_data.put('celltype/gene_expression_average', avg_exp.T)
        h5_data.put('celltype/gene_proportion_expression', frac_exp.T)
        h5_data.put('celltype/cell_count', ncells)

    print('Compute gene friends')
    # Compute top correlated genes
    counts = avg_exp.values.T
    genes = avg_exp.index

    # focus on highly expressed genes
    idx = counts.max(axis=0) > 30
    counts = counts[:, idx]
    genes = genes[idx]

    ngenes = len(genes)
    mu = counts.mean(axis=0)
    countsc = counts - mu
    sigma = counts.std(axis=0)

    friends = {}
    for i in range(ngenes):
        if ((i+1) % 100) == 0:
            print(f'Gene {i+1} out of {ngenes}', end='\r')
        corr_i = (countsc[:, i] * counts.T).mean(axis=1) / (sigma[i] * sigma.T)
        corr_i[np.isnan(corr_i)] = 0
        # Self is always the highest
        itop = np.argsort(corr_i)[::-1][1:6]

        # Cutoff at 20%
        tmp = corr_i[itop]
        tmp = tmp[tmp >= 0.2]
        itop = itop[:len(tmp)]

        friends_i = genes[itop]
        friends[genes[i]] = np.array(friends_i)
    print(f'Gene {i+1} out of {ngenes}')

    print('Store gene friends')
    output_fn = data_fdn + 'mouselemur_gene_friends.h5'
    with h5py.File(output_fn, 'w') as f:
        for gene, friends_i in friends.items():
            lmax = max(len(x) for x in friends_i)
            friends_i = friends_i.astype('S'+str(lmax))
            dset = f.create_dataset(gene, data=friends_i)


