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
from scipy import stats


data_fdn = '../WebPortal/static/scData/'


if __name__ == '__main__':

    fn_atlas = data_fdn + 'condensed_lung_atlas_in_cpm.h5'
    with h5py.File(fn_atlas) as f:
        dic = f['celltype']['gene_expression_average']
        genes = dic['axis0']
        genes = np.array([x.decode() for x in genes])
        counts = np.array(dic['block0_values']).astype(np.float32)

    # Fix this gene that has "inf" counts: CT010467.1
    counts[:, 5370] = 0
    counts = 1e6 * (counts.T / counts.sum(axis=1)).T

    # FIXME: focus on highly expressed genes
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
    print()

    output_fn = data_fdn + 'gene_friends.h5'
    with h5py.File(output_fn, 'w') as f:
        for gene, friends_i in friends.items():
            lmax = max(len(x) for x in friends_i)
            friends_i = friends_i.astype('S'+str(lmax))
            dset = f.create_dataset(gene, data=friends_i)
