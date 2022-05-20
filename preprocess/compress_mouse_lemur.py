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
    fn_atlas = '../data/tabula_microcebus/Lung_FIRM_hvg.h5ad'
    adata = anndata.read_h5ad(fn_atlas)

    # Find cell type column
    column = 'cell_ontology_class_v1'
    # NOTE: there is a free annotation column that is informative, basically
    # a higher resolution clustering but with some questionable things like
    # CD4+ CD8+ T cells. Let's go trad for now

    # Exclude unassigned/doublets
    unwanted_types = [
        'lymphocyte',
        'unassigned',
        'epithelial cell of uterus',
        'mature NK T cell',
        'endothelial cell',
        'granulocyte monocyte progenitor cell',
    ]
    adata = adata[~adata.obs[column].isin(unwanted_types)]

    # Set cell type order
    #celltypes = adata.obs[column].unique()
    celltypes = [
        'Fibroblast I',
        'Fibroblast II',
        'Mesothelial',
        'MyoF',
        'VSM',
        'Pericyte',
        'Venous',
        'Cap',
        'Art',
        'Lymphatic',
        'Schwann cell',
        'Plasmablast',
        'B cell',
        'NK cell',
        'IL cell',
        'pDC',
        'cDC',
        'DC X',
        'Alveolar mac',
        'Interstitial mac',
        'Monocyte',
        'Neutrophil',
        'Basophil',
        'Eosinophil',
        'Platelet',
        'Alveolar type I',
        'Alveolar type II',
        'Ciliated',
        'Club',
        'Brush',
        'Basal',
    ]

    # Merge some annotations together
    annotation_merges = {
        'Fibroblast I': ['fibroblast'],
        'Fibroblast II': ['fibroblast of lung'],
        'Mesothelial': ['mesothelial cell'],
        'MyoF': ['myofibroblast cell'],
        'VSM': ['vascular associated smooth muscle cell'],
        'Pericyte': ['pericyte cell'],
        'Schwann cell': [
            'myelinating Schwann cell',
            'non-myelinating Schwann cell',
        ],
        'Venous': ['vein endothelial cell'],
        'Cap': ['capillary endothelial cell'],
        'Art': ['endothelial cell of artery'],
        'Lymphatic': ['endothelial cell of lymphatic vessel'],
        'Plasmablast': ['plasma cell'],
        'NK cell': ['natural killer cell'],
        'IL cell': ['innate lymphoid cell'],
        'B cell': ['B cell'],
        'T cell': [
            'T cell',
            'CD4-positive, alpha-beta T cell',
            'CD8-positive, alpha-beta T cell',
        ],
        'red blood cell lineage': [
            'erythroid lineage cell',
            'hematopoietic precursor cell',
            'erythroid progenitor cell',
            'megakaryocyte progenitor cell',
        ],
        'Alveolar type I': ['type I pneumocyte'],
        'Alveolar type II': ['type II pneumocyte'],
        'Ciliated': ['lung ciliated cell'],
        'Club': ['club cell'],
        'Brush': ['brush cell of bronchus'],
        'Basal': ['basal cell of epithelium of bronchus'],
        'cDC': ['conventional dendritic cell'],
        'DC X': ['dendritic cell'],
        'pDC': ['plasmacytoid dendritic cell'],
        'Alveolar mac': ['alveolar macrophage'],
        'Interstitial mac': ['macrophage'],
        'Monocyte': ['monocyte'],
        'Neutrophil': ['neutrophil'],
        'Basophil': ['basophil'],
        'Eosinophil': ['eosinophil'],
        'Platelet': ['platelet'],
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
    fn_out = '../webapp/static/scData/mouselemur_condensed_lung_atlas_in_cpm.h5'
    with pd.HDFStore(fn_out) as h5_data:
        h5_data.put('celltype/gene_expression_average', avg_exp.T)
        h5_data.put('celltype/gene_proportion_expression', frac_exp.T)
        h5_data.put('celltype/cell_count', ncells)

