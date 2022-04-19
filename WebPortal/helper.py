from xml.dom import NotFoundErr
import h5py
import numpy as np
import pandas as pd
import json
import re
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist


fdn_data = "./static/scData/"
fn_atlas = fdn_data + "condensed_lung_atlas_in_cpm.h5"


def read_gene_order():
    with h5py.File(fn_atlas, "r") as h5_data:
        genes = np.array(h5_data['celltype']["gene_expression_average"]["axis0"].asstr())
    gene_order = pd.Series(data=np.arange(len(genes)), index=genes)
    return gene_order


gene_order = read_gene_order()


def read_counts_from_file(df_type, genes=None):
    '''Read the h5 file with the compressed atlas

    gives the name of dataset we want as an input
    celltype / celltype_dataset / celltype_dataset_timepoint
    '''
    with h5py.File(fn_atlas, "r") as h5_data:
        columns = np.array(h5_data[df_type]["gene_expression_average"]["axis1"].asstr())
        counts = h5_data[df_type]["gene_expression_average"]["block0_values"]
        if genes is not None:
            index = genes
            idx = gene_order.loc[genes].values
            counts = counts[:, idx]
        else:
            index = gene_order.index
        counts = np.array(counts).astype(np.float32)
        # NOTE: this is the same as gene order
        #index = np.array(h5_data[df_type]["gene_expression_average"]["axis0"].asstr())

    # FIXME: fix this in the h5 file...
    #Fix this gene that has "inf" counts: CT010467.1
    if genes is None:
        counts[:, 5370] = 0
        counts = 1e6 * (counts.T / counts.sum(axis=1)).T

    df = pd.DataFrame(
            data=counts.T,
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


def data_preprocessing(input_gene_names, df_type):
    '''Preprocess compressed atlas data

    this function is used when more multiple gene names are searched by the user
    always run the read_counts_from_file() to get dataframe in the right format before running this
    '''
    df = read_counts_from_file(df_type)
    df_genes = df.index  # all the genes that available in the current dataframe
    for search_gene in input_gene_names.split(","):
        if search_gene not in df_genes:
            return None

    # if user does not search for any specific genes, then plot everything
    if input_gene_names is None:
        plot_data = df.T
    else:
        a_gene_names = input_gene_names.split(",")
        plot_df = df.filter(items=a_gene_names, axis=0)
        plot_data = plot_df.T
    return plot_data


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
    for i in datasets:
        gene_exp_df = df_filtered.loc[i].unstack(1, fill_value=-1)
        gene_exp_df = gene_exp_df.loc[:, timepoint_order[i][::-1]]

        if datatype == "log10":
            gene_exp_df = np.log10(0.1 + gene_exp_df)
        if plottype == "hierachical":
            distance = pdist(gene_exp_df.values)
            # print(distance)
            Z = linkage(distance, optimal_ordering=True)
            new_order = leaves_list(Z)
            gene_exp_df = gene_exp_df.iloc[new_order]

        dic_per_dataset[i] = json.loads(gene_exp_df.to_json())
    return dic_per_dataset
    # for each of the dataframe in dic_per_dataset, we convert it into json format


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
                ms = 1
                opacity = 0.1
            elif nc < 20:
                ms = 8
                opacity = 0.5
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
            if gene not in h5_data.keys():
                continue
            friends_i = [x.decode() for x in h5_data[gene]]
            for friend in friends_i:
                # If a friend is already in the list of genes, leave it only
                # at the user-specified location
                if (friend not in friends) and (friend not in genes):
                    friends.append(friend)
    return ",".join(friends)
