# vim: fdm=indent
'''
author:     Fabio Zanini
date:       31/05/22
content:    Compress Want et al 2020 (snRNA-seq).
'''
import os
import sys
import numpy as np
import pandas as pd
import h5py
import anndata


data_fdn = '../webapp/static/scData/'


def correct_annotations(adata, column):
    '''Overwrite annotations for macrophages and DCs'''
    annotations = pd.read_csv(
        '../data/tabula_microcebus/reannotations.tsv',
        sep='\t',
        index_col=0,
    ).squeeze(axis=1)
    annotations = pd.Series(annotations, dtype='category')

    cats_old = adata.obs[column].cat.categories
    cats_new = list(set(annotations.cat.categories) - set(cats_old))
    adata.obs[column] = adata.obs[column].cat.add_categories(cats_new)
    annotations = annotations.cat.set_categories(adata.obs[column].cat.categories)

    adata.obs.loc[annotations.index, column] = annotations
    adata.obs[column] = adata.obs[column].cat.remove_unused_categories()


if __name__ == '__main__':

    # Load data
    print('Load single cell data')
    fn_atlas = '../data/Wang_et_al_2020/GSE161382_atlas.h5ad'
    adata = anndata.read_h5ad(fn_atlas)

    # Find cell type column
    column = 'celltype'

    # Exclude unassigned/doublets
    unwanted_types = [
        'unclassified',
        'AT1/AT2-like',
        'AT2/Club-like',
        'PNECs',
    ]
    adata = adata[~adata.obs[column].isin(unwanted_types)]

    # Set cell type order
    celltypes = [
        'Adventitial FB',
        'Alveolar FB',
        'MyoF',
        'ASM',
        'VSM',
        'Pericyte',
        'Lymphatic',
        'Arterial II',
        'Venous',
        'gCap',
        'Aerocyte',
        'B cell',
        'NK cell',
        'T cell',
        'DC',
        'Alveolar mac',
        'Interstitial mac',
        'Monocyte',
        'Mast',
        'Erythrocytes',
        'Alveolar type I',
        'Alveolar type II',
        'Ciliated',
        'Club',
        'Chondrocyte',
        'Goblet',
        'Basal',
    ]
    annotation_merges = {
        'Adventitial FB': ['matrix fibroblast 2'],
        'Alveolar FB': ['matrix fibroblast 1'],
        'MyoF': ['myofibroblasts'],
        'ASM': ['airway smooth muscle'],
        'VSM': ['vascular smooth muscle'],
        'Pericyte': ['pericytes'],
        'Lymphatic': ['lymphatics'],
        'Arterial II': ['arteries'],
        'Venous': ['veins'],
        'gCap': ['Cap1'],
        'Aerocyte': ['Cap2'],
        'B cell': ['B cells'],
        'NK cell': ['NK cells'],
        'T cell': ['T cells'],
        'DC': ['dendritic cells'],
        'Alveolar mac': ['Alveolar macrophages'],
        'Interstitial mac': ['Interstitial macrophages'],
        'Monocyte': ['monocytes'],
        'Mast': ['mast cells'],
        'Erythrocytes': ['enucleated erythrocytes'],
        'Alveolar type I': ['alveolar type 1 cells'],
        'Alveolar type II': ['alveolar type 2 cells'],
        'Ciliated': ['ciliated cells'],
        'Club': ['club cells'],
        'Chondrocyte': ['chondrocytes'],
        'Goblet': ['goblet cells'],
        'Basal': ['basal cells'],
    }
    adata.obs[column+'_compressed'] = ''
    for ct, csts in annotation_merges.items():
        adata.obs.loc[adata.obs[column].isin(csts), column+'_compressed'] = ct

    # Average, proportion expressing, number of cells
    genes = adata.var_names
    avgs = []
    ages = ['31wk', '3yr', '31yr']
    for age in ages:
        adata_age = adata[adata.obs['age'] == age]
        columns = [ct+'_Wang et al 2020_'+age for ct in celltypes]
        avg_exp = pd.DataFrame(
                np.zeros((len(genes), len(celltypes)), np.float32),
                index=genes,
                columns=columns)
        frac_exp = avg_exp.copy()
        ncells = pd.Series(np.zeros(len(celltypes), np.int64), index=columns)
        for ct in celltypes:
            # Avg gene expression (compute cpm after averaging)
            label = ct+'_Wang et al 2020_'+age
            idx = adata_age.obs[column+'_compressed'] == ct
            avg_exp[label] = np.asarray(adata_age[idx].X.mean(axis=0))[0]
            avg_exp[label] *= 1e6 / avg_exp[label].sum()
            # Proportion expressing
            frac_exp[label] = np.asarray((adata_age[idx].X > 0).mean(axis=0))[0]
            # Number of cells
            ncells[label] = (adata_age.obs[column+'_compressed'] == ct).sum()
        avgs.append({
            'age': age,
            'gene_expression_average': avg_exp,
            'gene_proportion_expression': frac_exp,
            'cell_count': ncells,
        })

    # Restructure as a single dataframe for each category
    # celltype_dataset_timepoint
    avg_exp = pd.concat([x['gene_expression_average'] for x in avgs]).fillna(0).astype(np.float32)
    frac_exp = pd.concat([x['gene_proportion_expression'] for x in avgs]).fillna(0).astype(np.float32)
    ncells = pd.concat([x['cell_count'] for x in avgs]).fillna(0)

    # Store to file (pandas can read/write from a previously touched file)
    print('Store compressed atlas')
    fn_out = data_fdn + 'human_condensed_lung_atlas_in_cpm.h5'
    with pd.HDFStore(fn_out) as h5_data:
        h5_data.put('celltype_dataset_timepoint/gene_expression_average', avg_exp)
        h5_data.put('celltype_dataset_timepoint/gene_proportion_expression', frac_exp)
        h5_data.put('celltype_dataset_timepoint/cell_count', ncells)

