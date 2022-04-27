# vim: fdm=indent
'''
author:     Fabio Zanini
date:       16/04/22
content:    Precompute list of top correlated genes for "friends" button.
'''
import os
import sys
import h5py
import numpy as np
import pandas as pd
from scipy import stats


data_fdn = '../WebPortal/static/scData/'


def plot_fcs(fcmax, fcmean, ct, cands):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.scatter(fcmax, fcmean)
    ax.set_xlabel('fold change from max')
    ax.set_xlabel('fold change from mean')
    ax.set_title(ct)
    ax.grid(True)

    for gene, i in cands.items():
        ax.text(fcmax[i], fcmean[i], gene, ha='center', va='center')

    fig.tight_layout()

    return fig, ax


if __name__ == '__main__':

    fn_atlas = data_fdn + 'condensed_lung_atlas_in_cpm.h5'
    with h5py.File(fn_atlas) as f:
        dic = f['celltype']['gene_expression_average']
        genes = np.array(dic['axis0'].asstr())
        celltypes = np.array(dic['axis1'].asstr())
        counts = np.array(dic['block0_values']).astype(np.float32)

    # Fix this gene that has "inf" counts: CT010467.1
    counts[:, 5370] = 0
    counts = 1e6 * (counts.T / counts.sum(axis=1)).T

    # FIXME: focus on highly expressed genes
    idx = counts.max(axis=0) > 30
    counts = counts[:, idx]
    genes = genes[idx]

    ncts = len(celltypes)
    markerd = {}
    for i, ct in enumerate(celltypes):
        exp_ct = np.log10(0.5 + counts[i])
        exp_other = np.log10(0.5 + counts[[j for j in range(ncts) if j != i]])
        fcmax = exp_ct - exp_other.max(axis=0)
        fcmean = exp_ct - exp_other.mean(axis=0)

        idx_cand = list(set(np.argsort(fcmax)[-10:]) | set(np.argsort(fcmean)[-10:]))
        cands = pd.Series(idx_cand, index=genes[idx_cand])
        markerd[ct] = list(cands.index)

        #fig, ax = plot_fcs(fcmax, fcmean, ct, cands)

    output_fn = data_fdn + 'marker_genes.h5'
    with h5py.File(output_fn, 'w') as f:
        group = f.create_group('celltype')
        for i, ct in enumerate(celltypes):
            genes_i = np.asarray(markerd[ct])
            lmax = max(len(x) for x in genes_i)
            genes_i = genes_i.astype('S'+str(lmax))
            group.create_dataset(ct, data=genes_i)
