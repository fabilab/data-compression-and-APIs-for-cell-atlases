# vim: fdm=indent
'''
author:     Fabio Zanini
date:       19/04/22
content:    Test plotly dotplots in python, exporting the JSON to the frontend.
'''
import os
import sys
import h5py
import numpy as np
import pandas as pd
from plotly import graph_objects as go
from plotly import colors as pcolors


data_fdn = '../WebPortal/static/scData/'


if __name__ == '__main__':


    # Get the data
    fn_atlas = data_fdn + 'condensed_lung_atlas_in_cpm.h5'
    with h5py.File(fn_atlas) as f:
        dic = f['celltype_dataset_timepoint']

        labels = np.array(dic['cell_count']['axis0'].asstr())
        values = np.array(dic['cell_count']['block0_values'])[0]
        ncells = pd.Series(data=values, index=labels)

        genes = np.array(dic['gene_expression_average']['axis0'].asstr())
        labels = np.array(dic['gene_expression_average']['axis1'].asstr())
        counts = np.array(dic['gene_expression_average']['block0_values']).astype(np.float32)

    # Fix this gene that has "inf" counts: CT010467.1
    counts[:, 5370] = 0
    counts = 1e6 * (counts.T / counts.sum(axis=1)).T

    counts = pd.DataFrame(counts.T, index=genes, columns=labels)

    # NOTE: Now both are in pandas dataframes

    # Choose a gene and visualization arguments
    gene = 'Car4'
    use_log = False
    use_hierarchical = False

    countg = counts.loc[gene]
    if use_log:
        countg = np.log10(0.1 + countg)

    # Sort the rows
    timepoint_order = ['E18.5', 'P1', 'P3', 'P7', 'P14', 'P21', '3m', '18m', '24m']
    timepoint_orderd = {x: i for i, x in enumerate(timepoint_order)}
    dataset_order = ['ACZ', 'Hurskainen2021', 'TMS']
    dataset_orderd = {x: i for i, x in enumerate(dataset_order)}
    ncells.index = ncells.index.str.split('_', 1, expand=True)
    ncells = ncells.unstack(0, fill_value=0)
    ncells[['dataset', 'timepoint']] = ncells.index.str.split('_', expand=True).tolist()
    ncells['timepoint_order'] = ncells['timepoint'].map(timepoint_orderd)
    ncells['dataset_order'] = ncells['dataset'].map(dataset_orderd)
    ncells.sort_values(['timepoint_order', 'dataset_order'], inplace=True)
    countg.index = countg.index.str.split('_', 1, expand=True)
    countg = countg.unstack(0, fill_value=-1).loc[ncells.index]
    ncells = ncells.loc[:, countg.columns]

    # Sort the columns
    if use_hierarchical:
        # TODO
        pass

    ## Normalize ncells by row (decide what to do, but roughly ok)
    #ncells = (ncells.astype(np.float32).T / ncells.sum(axis=1)).T

    # Construct data structures (lists) for plotly
    plot_data = {'x': [], 'y': [], 'marker': {'color': [], 'size': [], 'opacity': []}}
    colormax = countg.values.max()
    color_steps = pcolors.colorbrewer.Reds
    ncolors = len(color_steps)
    for i, label in enumerate(ncells.index):
        for j, celltype in enumerate(countg.columns):
            plot_data['x'].append(j)
            plot_data['y'].append(i)
            nc = ncells.at[label, celltype]
            if nc < 5:
                ms = 1
                opacity = 0.1
            elif nc < 20:
                ms = 10
                opacity = 0.5
            else:
                ms = 30
                opacity = 0.9
            plot_data['marker']['size'].append(ms)
            plot_data['marker']['opacity'].append(opacity)
            ge = min(ncolors - 1, int(ncolors * countg.at[label, celltype] / colormax))
            plot_data['marker']['color'].append(ge)

    plot_data['mode'] = 'markers'
    plot_data['marker_symbol'] = 'square'

    layout = {
        'title': {
            'text': f'{gene} expression over time',
            'x': 0.5,
            'xanchor': 'center',
        },
        'xaxis': {
            'title': 'Cell subtype',
            #'type': 'category',
            #'nticks': countg.shape[1],
            #'tick0': 0,
            #'dtick': 1,
            'ticktext': list(countg.columns),
            'tickvals': list(range(countg.shape[1])),
            'tickangle': 60,
        },
        'yaxis': {
            #'nticks': countg.shape[0],
            #'tick0': 0,
            #'dtick': 1,
            'ticktext': list(countg.index),
            'tickvals': list(range(countg.shape[0])),
            'autorange': 'reversed',
        },
    }

    fig = go.Figure(
        data=[go.Scatter(**plot_data)],
        layout=layout,
    )
    plot_json = fig.to_plotly_json()
