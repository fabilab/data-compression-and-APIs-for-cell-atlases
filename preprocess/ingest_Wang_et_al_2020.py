# vim: fdm=indent
'''
author:     Fabio Zanini
date:       16/05/22
content:    Take a look at Wang et al. 2020, snRNA-Seq on 3 ages: 30wGA, 3y, 30y,
            with three lungs per age group.
'''
import os
import sys
import numpy as np
import pandas as pd
import h5py
import anndata

import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append('/home/fabio/university/postdoc/singlet/build/lib')
import singlet
import scanpy as sc


data_fdn = '../webapp/static/scData/'


if __name__ == '__main__':

    # Load data
    print('Load single cell data')
    fn_atlas_pfx = '../data/Wang_et_al_2020/GSE161382_'
    adata = anndata.read_mtx(f'{fn_atlas_pfx}counts.mtx.gz').T
    obs = pd.read_csv(f'{fn_atlas_pfx}barcodes.tsv.gz', sep='\t', compression='gzip')
    adata.obs_names = obs['x'].values
    obs = pd.read_csv(f'{fn_atlas_pfx}metadata.txt.gz', sep=' ', compression='gzip')
    adata.obs = obs
    var = pd.read_csv(f'{fn_atlas_pfx}features.tsv.gz', sep='\t', compression='gzip')
    adata.var_names = var['x'].values

    umap = pd.read_csv(f'{fn_atlas_pfx}UMAP_coord.tsv.gz', sep='\t', compression='gzip')
    adata.obsm['X_umap'] = umap.values

    adata.write(f'{fn_atlas_pfx}atlas.h5ad')
