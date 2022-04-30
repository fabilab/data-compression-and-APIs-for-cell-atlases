from xml.dom import NotFoundErr
import h5py
import numpy as np
import pandas as pd
import json
import re
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

from validation.celltypes import adjust_celltypes


fdn_data = "./static/scData/"
# FIXME
fn_atlas = fdn_data + "condensed_lung_atlas_in_cpm.h5"
fn_atlas = fdn_data + "condensed_lung_atlas_ordered.h5"


def read_gene_order(data_type='celltype'):
    with h5py.File(fn_atlas, "r") as h5_data:
        genes = np.array(h5_data[data_type]["gene_expression_average"]["axis0"].asstr())
    gene_order = pd.Series(data=np.arange(len(genes)), index=genes)
    return gene_order


def read_cell_types():
    with h5py.File(fn_atlas, "r") as h5_data:
        celltypes = np.array(h5_data['celltype']["gene_expression_average"]["axis1"].asstr())
    celltypes, _ = adjust_celltypes(celltypes)
    return celltypes


gene_order = read_gene_order()
celltypes = read_cell_types()


def read_counts_from_file(df_type, genes=None):
    '''Read the h5 file with the compressed atlas

    gives the name of dataset we want as an input
    celltype / celltype_dataset / celltype_dataset_timepoint
    '''
    gene_order = read_gene_order(data_type=df_type)

    with h5py.File(fn_atlas, "r") as h5_data:
        columns = np.array(h5_data[df_type]["gene_expression_average"]["axis1"].asstr())
        columns, idx_cols = adjust_celltypes(columns)
        counts = h5_data[df_type]["gene_expression_average"]["block0_values"]
        if genes is not None:
            index = []
            for gene in genes:
                candidates = gene_order.index[gene_order.index.str.match(gene)]
                if len(candidates) == 0:
                    raise KeyError('Gene match not found: '+gene)
                for cand in candidates:
                    if cand not in index:
                        index.append(cand)
                
            idx = gene_order.loc[index].values
            # Slices of h5 must be ordered...
            tmp = pd.Series(idx, index=np.arange(len(idx))).to_frame(name='idx')
            tmp = tmp.sort_values('idx')
            tmp['new_order'] = np.arange(len(tmp))
            counts = np.array(counts[:, tmp['idx'].values])
            # ... now reorder them according to the user's needs
            counts = counts[:, tmp.sort_index()['new_order'].values]
        else:
            index = gene_order.index
        counts = np.array(counts).astype(np.float32)

    # FIXME: fix this in the h5 file...
    #Fix this gene that has "inf" counts: CT010467.1
    if genes is None:
        counts[:, 5370] = 0
        counts = 1e6 * (counts.T / counts.sum(axis=1)).T

    df = pd.DataFrame(
            data=counts[idx_cols].T,
            index=index,
            columns=columns,
            )
    return df


def read_number_cells_from_file(df_type):
    '''Read the number of cells per time point per dataset from file'''
    with h5py.File(fn_atlas) as f:
        dic = f[df_type]

        labels = np.array(dic['cell_count']['axis0'].asstr())
        values = np.array(dic['cell_count']['block0_values'])[0]
        ncells = pd.Series(data=values, index=labels)

    return ncells


def dataset_by_timepoint(genename, df_type, datatype, plottype):
    '''get the cell type dataset timepoint as a dictionary

    user input a gene name.
    for each unique dataset of this gene, we plot a heatmap of all timepoint vs celltypes
    select and pre-preprocessing data
    '''
    df_filtered = read_counts_from_file(df_type, genes=[genename]).iloc[0]

    # split into separate unstacked dataframes using a multi-index
    df_filtered.index = df_filtered.index.str.split('_', expand=True).swaplevel(0, 1)

    # for each dataset name, we construct a new dataframe for it (celltypes x timepoints)
    datasets = set(df_filtered.index.get_level_values(0))
    dic_per_dataset = {}
    timepoint_order = {
        'ACZ': ['E18.5', 'P1', 'P7', 'P21'],
        'TMS': ['3m', '18m', '24m'],
        'Hurskainen2021': ['P3', 'P7', 'P14'],
    }
    for dsname in datasets:
        gene_exp_df = df_filtered.loc[dsname].unstack(1, fill_value=-1)
        gene_exp_df = gene_exp_df.loc[:, timepoint_order[dsname][::-1]]

        if datatype == "log10":
            gene_exp_df = np.log10(0.1 + gene_exp_df)
        if plottype == "hierachical":
            distance = pdist(gene_exp_df.values)
            # print(distance)
            Z = linkage(distance, optimal_ordering=True)
            new_order = leaves_list(Z)
            gene_exp_df = gene_exp_df.iloc[new_order]

        # convert dataframe into json format
        dic_per_dataset[dsname] = json.loads(gene_exp_df.to_json())
    return dic_per_dataset


def get_big_heatmap(gene, use_log, use_hierarchical):
    '''Get JS plotly code for big heatmap

    NOTE: this technically breaks the API concept. Let's keep it in mind
    and see what we can do. The positive is that a bunch of computations
    happen on the server... wait, is that good?
    '''
    from plotly import colors as pcolors

    ncells = read_number_cells_from_file('celltype_dataset_timepoint')
    countg = read_counts_from_file('celltype_dataset_timepoint', genes=[gene]).iloc[0]

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

    # Sort the columns
    if use_hierarchical:
        distance = pdist(countg.T.values)
        # print(distance)
        Z = linkage(distance, optimal_ordering=True)
        new_order = leaves_list(Z)
        countg = countg.iloc[:, new_order]

    ncells = ncells.loc[:, countg.columns]

    # NOTE: this is the part that we could outsource to the frontend
    # Construct data structures (lists) for plotly
    plot_data = {'x': [], 'y': [], 'text': [], 'marker': {'color': [], 'size': [], 'opacity': []}}
    colormax = countg.values.max()
    color_steps = pcolors.colorbrewer.Reds
    ncolors = len(color_steps)
    for i, label in enumerate(ncells.index):
        for j, celltype in enumerate(countg.columns):
            nc = ncells.at[label, celltype]
            ge = countg.at[label, celltype]
            plot_data['x'].append(celltype)
            plot_data['y'].append(label)
            plot_data['text'].append(str(ge))
            if nc < 5:
                ms = 8
                opacity = 0.5
            elif nc < 40:
                ms = 13
                opacity = 0.7
            else:
                ms = 20
                opacity = 0.9
            plot_data['marker']['size'].append(ms)
            plot_data['marker']['opacity'].append(opacity)
            color = min(ncolors - 1, int(ncolors * ge / colormax))
            plot_data['marker']['color'].append(color)

    plot_data['mode'] = 'markers'
    plot_data['marker_symbol'] = 'square'

    xticks = list(countg.columns)
    yticks = list(ncells.index)
    yticktext = []
    for i, label in enumerate(ncells.index):
        tp = label.split('_')[1]
        if (i == 0) or (tp != yticktext[-1]):
            yticktext.append(tp)
        else:
            yticktext.append('')

    return {
        'data': plot_data,
        'xticks': xticks,
        'yticks': yticks,
        'yticktext': yticktext,
        }


def get_friends(genes):
    '''Get genes that are "friends" (correlated) with a list of genes'''
    friends = []
    with h5py.File(fdn_data + "gene_friends.h5", "r") as h5_data:
        for gene in genes:
            friends.append(gene)
            # We only took decently expressed genes, but this might appear
            # as a bug if we are not careful...
            if gene not in h5_data.keys():
                continue
            friends_i = [x.decode() for x in h5_data[gene]]
            for friend in friends_i:
                # If a friend is already in the list of genes, leave it only
                # at the user-specified location
                if (friend not in friends) and (friend not in genes):
                    friends.append(friend)
    return ",".join(friends)


def get_marker_genes(celltypes):
    '''Get markers for cell types'''
    # Sometimes we get a string
    if isinstance(celltypes, str):
        celltypes = celltypes.split(',')

    markers = []
    with h5py.File(fdn_data + "marker_genes.h5", "r") as h5_data:
        group = h5_data['celltype']
        for celltype in celltypes:
            if celltype not in group:
                return ''
            genes_i = list(np.array(group[celltype].asstr()))
            for gene in genes_i:
                if gene not in markers:
                    markers.append(gene)

    return ",".join(markers)
    

def get_data_hyperoxia(data_type, genes=None):
    '''Get heatmap data for hyperoxia'''
    df_ho = read_counts_from_file(
        'celltype_dataset_timepoint_hyperoxia',
        genes=genes,
        )

    if data_type == "log2FC":
        df_normal = read_counts_from_file(
                'celltype_dataset_timepoint',
                genes=genes,
            )
        # Restrict to hyperoxia celltypes, datasets, and timepoints
        # NOTE: no dataset took hyperoxia from a timepoint that has no normal.
        # However, come cell types are hyperoxia specific
        for key in df_ho:
            if key not in df_normal.columns:
                # Default to zero expression
                df_normal[key] = 0
        df_normal = df_normal.loc[df_ho.index, df_ho.columns]
        df_ho = np.log2(df_ho + 0.5) - np.log2(df_normal + 0.5)

    elif data_type == "log10":
        df_ho = np.log10(df_ho + 0.5)

    # Split by dataset and timepoint, and unstack
    df_ho = df_ho.T
    datasets = ['ACZ', 'Hurskainen2021']
    timepointd = {'ACZ': ['P7'], 'Hurskainen2021': ['P3', 'P7', 'P14']}
    result = []
    for ds in datasets:
        for tp in timepointd[ds]:
            item = {
                'dataset': ds,
                'timepoint': tp,
                'data_scale': data_type,
            }
            dfi = df_ho.loc[df_ho.index.str.endswith(f'{ds}_{tp}')]
            dfi.index = dfi.index.str.split('_', expand=True).get_level_values(0)
            item['data'] = dfi
            result.append(item)
    return result

