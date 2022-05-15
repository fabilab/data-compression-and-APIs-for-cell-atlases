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

    # Average
    column = 'XXX'
    celltypes = adata.obs[column].unique()
    genes = adata.var_names
    avg_exp = pd.DataFrame(
            np.zeros((len(genes), len(celltypes)), np.float32),
            index=genes,
            columns=celltypes)
    frac_exp = avg_exp.copy()
    for ct in celltypes:
        # Avg gene expression (compute cpm after averaging)
        avg_exp[ct] = adata[adata.obs[column] == ct].X.mean(axis=0)
        avg_exp[ct] *= 1e6 / avg_exp[ct].sum()
        # Proportion expressing
        frac_exp[ct] = (adata[adata.obs[column] == ct].X > 0).mean(axis=0)

    # Store to file
    fn_out = '../webapp/static/scData/mouselemur_condensed_lung_atlas_in_cpm.h5'
    with pd.HDFStore(fn_out) as h5_data:
        h5_data.put('celltype/gene_expression_average', avg_exp)
        h5_data.put('celltype/gene_proportion_expression', frac_exp)

