from xml.dom import NotFoundErr
import h5py
import numpy as np
import pandas as pd
import json
import pickle
import re
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

from validation.celltypes import (
        adjust_celltypes,
        rename_celltypes,
        validate_correct_celltypestr,
        )
from validation.differential_expression import get_deg_conditions


fdn_data = "./static/scData/"
fn_atlasd = {
    'mouse': fdn_data + "condensed_lung_atlas_ordered.h5",
    'lemur': fdn_data + "mouselemur_condensed_lung_atlas_in_cpm.h5",
    'human': fdn_data + "human_condensed_lung_atlas_in_cpm.h5",
}
fn_GO = fdn_data + 'mouse_GO_tables.pkl'


def read_gene_order(data_type='celltype', species='mouse'):
    fn_atlas = fn_atlasd[species]
    with h5py.File(fn_atlas, "r") as h5_data:
        genes = np.array(h5_data[data_type]["gene_expression_average"]["axis0"].asstr())
    gene_order = pd.Series(data=np.arange(len(genes)), index=genes)
    return gene_order


def read_cell_types(species='mouse'):
    fn_atlas = fn_atlasd[species]
    with h5py.File(fn_atlas, "r") as h5_data:
        celltypes = np.array(h5_data['celltype']["gene_expression_average"]["axis1"].asstr())
    celltypes, _ = adjust_celltypes(celltypes)
    return celltypes


gene_orderd = {key: read_gene_order(species=key) for key in fn_atlasd}
celltypes = read_cell_types()


def read_counts_from_file(df_type, genes=None, species='mouse', missing='throw'):
    '''Read the h5 file with the compressed atlas

    gives the name of dataset we want as an input
    celltype / celltype_dataset / celltype_dataset_timepoint
    '''
    gene_order = gene_orderd[species]

    fn_atlas = fn_atlasd[species]
    with h5py.File(fn_atlas, "r") as h5_data:
        columns = np.array(h5_data[df_type]["gene_expression_average"]["axis1"].asstr())
        # For mouse, the cell types used for compression are not the final
        # ones we want to show
        if species == "mouse":
            columns, idx_cols = adjust_celltypes(columns)
        else:
            idx_cols = np.arange(len(columns))

        counts = h5_data[df_type]["gene_expression_average"]["block0_values"]
        if genes is not None:
            index = []
            for gene in genes:
                # Try perfect match
                if gene in gene_order.index:
                    if gene not in index:
                        index.append(gene)
                    continue

                # Try imperfect match
                candidates = gene_order.index[gene_order.index.str.match(gene)]
                if len(candidates) == 0:
                    if missing == 'throw':
                        raise KeyError('Gene match not found: '+gene)
                    else:
                        continue
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

    df = pd.DataFrame(
            data=counts[idx_cols].T,
            index=index,
            columns=columns,
            )

    return df


def read_number_cells_from_file(df_type, species='mouse'):
    '''Read the number of cells per time point per dataset from file'''
    fn_atlas = fn_atlasd[species]
    with h5py.File(fn_atlas) as f:
        dic = f[df_type]

        labels = np.array(dic['cell_count']['axis0'].asstr())
        values = np.array(dic['cell_count']['block0_values'])[0]
        ncells = pd.Series(data=values, index=labels)

    return ncells


def dataset_by_timepoint(genename, df_type, datatype, plottype, species='mouse'):
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
            Z = linkage(distance, optimal_ordering=True)
            new_order = leaves_list(Z)
            gene_exp_df = gene_exp_df.iloc[new_order]

        # convert dataframe into json format
        dic_per_dataset[dsname] = json.loads(gene_exp_df.to_json())
    return dic_per_dataset


def dataset_unified(gene, species='mouse'):
    '''Get JS plotly code for big heatmap

    NOTE: this technically breaks the API concept. Let's keep it in mind
    and see what we can do. The positive is that a bunch of computations
    happen on the server... wait, is that good?
    '''
    from plotly import colors as pcolors

    ncells = read_number_cells_from_file('celltype_dataset_timepoint')
    countg = read_counts_from_file('celltype_dataset_timepoint', genes=[gene]).iloc[0]

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

    # Set canonical celltype order
    celltypes_adj, idx = adjust_celltypes(ncells.columns)
    ncells = ncells.iloc[:, idx]
    countg = countg.iloc[:, idx]
    ncells.columns = celltypes_adj
    countg.columns = celltypes_adj

    # Sort the columns
    distance = pdist(countg.T.values)
    Z = linkage(distance, optimal_ordering=True)
    new_order = leaves_list(Z)

    celltypes = countg.columns.tolist()
    celltypes_hierarchical = countg.columns[new_order].tolist()
    row_labels = ncells.index.tolist()

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
        'gene': gene,
        'gene_expression': countg.T.to_dict(),
        'ncells': ncells.T.to_dict(),
        'celltypes': celltypes,
        'celltypes_hierarchical': celltypes_hierarchical,
        'row_labels': row_labels,
        'xticks': xticks,
        'yticks': yticks,
        'yticktext': yticktext,
        }


def get_friends(genes, species="mouse"):
    '''Get genes that are "friends" (correlated) with a list of genes'''
    friends = []
    if len(genes) == 0:
        return ''

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


def get_marker_genes(celltypes, species='mouse'):
    '''Get markers for cell types'''
    from validation.celltypes import celltype_dict
    
    celltype_dict_inv = {val: key for key, val in celltype_dict.items()}
    celltype_dict_inv.update({ct: ct for ct in celltype_dict_inv.values()})

    # Sometimes we get a string
    if isinstance(celltypes, str):
        celltypes = celltypes.split(',')

    markers = []
    with h5py.File(fdn_data + "marker_genes.h5", "r") as h5_data:
        group = h5_data['celltype']
        for celltype in celltypes:
            if celltype not in celltype_dict_inv:
                return ''
            raw_ct = celltype_dict_inv[celltype]
            genes_i = list(np.array(group[raw_ct].asstr()))
            for gene in genes_i:
                if gene not in markers:
                    markers.append(gene)

    return ",".join(markers)


def get_de_url(suffix, kind='both'):
    '''Get differential expression url for up-, downregulated genes, or both'''
    from urllib.parse import urlencode

    deg_dict = get_deg_conditions(suffix)
    deg_request = {
        'comparison': deg_dict['kind'],
        'ct1': deg_dict['conditions'][0]['celltype'],
        'ct2': deg_dict['conditions'][1]['celltype'],
        'ds1': deg_dict['conditions'][0]['dataset'],
        'ds2': deg_dict['conditions'][1]['dataset'],
        'tp1': deg_dict['conditions'][0]['timepoint'],
        'tp2': deg_dict['conditions'][1]['timepoint'],
        'dis1': deg_dict['conditions'][0].get('disease', 'normal'),
        'dis2': deg_dict['conditions'][1].get('disease', 'normal'),
        'kind': kind,
    }
    deg_request['n_genes'] = deg_request.get('n_genes', 20)
    url = urlencode(deg_request)
    return url


def get_data_differential(conditions, kind='both', n_genes=20, genes=None, species='mouse'):
    '''Get differentially expressed, up- or downregulated genes

    the special phease " in " separated the split (e.g. hyperoxia) from
    the celltype/timepoint/dataset specification (e.g. basophil ACZ P7)
    '''
    # If requested, compute DEGs on the fly
    if genes is None:
        dfs = []
    # In any case, collect data to show in the heatmap
    dfs_alltypes = []
    for condition in conditions:
        if condition.get('disease', 'normal') == "hyperoxia":
            dfi = read_counts_from_file(
             'celltype_dataset_timepoint_hyperoxia',
             genes=genes,
             ) 
        else:
            dfi = read_counts_from_file(
                'celltype_dataset_timepoint',
                genes=genes,
                )
        
        # Get data for only this condition if computation of DEGs is requested
        if genes is None:
            idx_col = '_'.join([
                condition['celltype'],
                condition['dataset'],
                condition['timepoint'],
                ])
            if idx_col not in dfi.columns:
                raise ValueError(f'Condition not found: {idx_col}')
            dfi_cond = dfi.loc[:, idx_col]
            dfs.append(dfi_cond)

        # In any case, get data for this condition, but all cell types
        idx_col = dfi.columns.str.contains(
                '_'+condition['dataset']+'_'+condition['timepoint'])
        dfi_alltypes = dfi.loc[:, idx_col]
        # Rename the columns as the cell types
        dfi_alltypes.columns = [x.split('_')[0] for x in dfi_alltypes.columns]
        dfs_alltypes.append(dfi_alltypes)

    if genes is None:
        # Get log2 fold change
        log2fc = np.log2(dfs[0] + 0.5) - np.log2(dfs[1] + 0.5)

        # Get DEGs
        deg_up = log2fc.nlargest(n_genes).index.tolist()
        deg_down = log2fc.nsmallest(n_genes).index.tolist()
        genes = deg_up + deg_down[::-1]

    # Use or recycle data to make heatmap
    # FIXME: make an expected total ordering including disease
    celltypes = dfs_alltypes[0].columns.tolist()
    for ct in dfs_alltypes[1].columns:
        if ct not in celltypes:
            celltypes.append(ct)
    dfs_alltypes = [dfi.loc[genes] for dfi in dfs_alltypes]
    for i, dfi in enumerate(dfs_alltypes):
        for ct in celltypes:
            if ct not in dfi.columns:
                dfi[ct] = 0
        dfs_alltypes[i] = dfi

    return dfs_alltypes


def get_data_hyperoxia(genes=None):
    '''Get heatmap data for hyperoxia'''
    df_ho = read_counts_from_file(
        'celltype_dataset_timepoint_hyperoxia',
        genes=genes,
        )

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

    # Split by dataset and timepoint, and unstack
    df_dict = {'data': df_ho.T, 'data_baseline': df_normal.T}
    result = []

    datasets = ['ACZ', 'Hurskainen2021']
    timepointd = {'ACZ': ['P7'], 'Hurskainen2021': ['P3', 'P7', 'P14']}
    for ds in datasets:
        for tp in timepointd[ds]:
            item = {
                'dataset': ds,
                'timepoint': tp,
            }
            for key, df in df_dict.items():
                dfi = df.loc[df.index.str.endswith(f'{ds}_{tp}')]
                dfi.index = dfi.index.str.split('_', expand=True).get_level_values(0)
                # Adjust cell type names and order
                celltypes_adj, idx = adjust_celltypes(dfi.index)
                dfi = dfi.iloc[idx]
                item[key] = dfi
            result.append(item)

    return result


def get_celltype_abundances(timepoint, dataset='ACZ', kind='qualitative', species='mouse'):
    '''Get cell type abundances at a certain time point'''
    ncells = read_number_cells_from_file('celltype_dataset_timepoint')

    idx = ncells.index.str.contains('_'+dataset+'_'+timepoint)
    ncells_tp = ncells[idx]

    ncells_tp.index = [x.split('_')[0] for x in ncells_tp.index]

    ntot = ncells_tp.sum()
    fracs = 1.0 * ncells_tp / ntot

    if kind == 'quantitative':
        return fracs

    # Sort them naturally
    fracs.sort_values(ascending=False, inplace=True)

    # Rename to pretty names
    fracs.index = rename_celltypes(fracs.index)

    result = {
        'major': [],
        'minor': [],
        'rare': [],
        'missing': [],
    }
    for ct, fr in fracs.items():
        if fr > 0.1:
            result['major'].append(ct)
        elif fr > 0.03:
            result['minor'].append(ct)
        elif fr > 0:
            result['rare'].append(ct)
        else:
            result['missing'].append(ct)

    return result


def get_gene_ids(genes, species='mouse'):
    '''Get the ids (MGI/GeneCards) of a list of genes'''
    if species == 'mouse':
        fn = fdn_data+'mouse_gene_names.tsv'
        df = pd.read_csv(fn, sep='\t', index_col=0)
        mgi_dict = df['MGI_id'].to_dict()
        human_dict = df['HumanGeneName'].fillna('').to_dict()
        id_dict = {}
        for gene in genes:
            new_name = human_dict.get(gene, '')
            if new_name == '':
                new_name = mgi_dict.get(gene, '')
            id_dict[gene] = new_name
    elif species == 'human':
        # GeneCards uses gene names as URL ids
        id_dict = {g: g for g in genes}
    elif species == 'lemur':
        # Seems like lemur and human share gene names
        return get_gene_ids(genes, species='human')
    else:
        raise NotImplementedError

    return id_dict


def get_gene_ontology_terms(genes, species='mouse'):
    '''Get the GO terms for these genes'''
    genesd = get_orthologs(genes, species, 'mouse')

    print(genesd)

    with open(fn_GO, 'rb') as fin:
        table = pickle.load(fin)['genes_to_GO']

    result = {}
    for gene_mouse, gene_orig in zip(genesd['mouse'], genesd[species]):
        if gene_mouse in table:
            result[gene_orig] = table[gene_mouse]

    return result


def get_genes_in_GO_term(go_term, species='mouse'):
    '''Get the genes in a GO term'''
    with open(fn_GO, 'rb') as fin:
        table = pickle.load(fin)['GO_to_genes']

    if go_term not in table:
        raise KeyError("GO term not found")

    genes = table[go_term]
    if species != 'mouse':
        genes = get_orthologs(genes, 'mouse', species)[species]

    return genes


def get_orthologs(genes, species, new_species):
    '''Connect orthologs from a species to another'''
    # Mouse lemur seems to use same gene names as human... is it a primate thing?
    if species == new_species:
        return {
            species: list(genes),
            new_species: list(genes),
            }
    elif species == 'lemur':
        dic = get_orthologs(genes, 'human', new_species)
        dic[species] = dic.pop('human')
        return dic
    elif new_species == 'lemur':
        dic = get_orthologs(genes, species, 'human')
        dic[new_species] = dic.pop('human')
        return dic
    elif (species == 'mouse') and (new_species == 'human'):
        fn = fdn_data+'mouse_gene_names.tsv'
        conv = pd.read_csv(fn, sep='\t', index_col=0, usecols=[0, 3])
        conv = conv.squeeze("columns")
        genes = [g for g in genes if g in conv.index]
        dic = conv.loc[genes].fillna('')
    elif (species == 'human') and (new_species == 'mouse'):
        fn = fdn_data+'human_mouse_gene_orthologs.tsv'
        conv = pd.read_csv(fn, sep='\t', index_col=0)
        conv = conv.squeeze("columns")
        genes = [g for g in genes if g in conv.index]
        dic = conv.loc[genes].fillna('')
    else:
        raise NotImplementedError

    dic = dic[dic != '']
    return {
        species: list(dic.index),
        new_species: list(dic.values),
    }


def guess_genes_species(genes):
    '''Guess the species from the genes'''
    species = ('mouse',)
    for gene in genes:
        if len(gene) > 1:
            break
    else:
        return species

    if gene == gene.upper():
        species = ('human', 'lemur')
    return species


def get_data_species_comparison(species, species_baseline, genes):
    '''Get dataframe of average expression in two species.

    Genes with no ortholog will be ignored.
    '''
    gene_species = guess_genes_species(genes)

    # Several species could have the same
    if (species in gene_species) and (species_baseline in gene_species):
        genes_baseline = genes
    elif species in gene_species:
        dic = get_orthologs(genes, species, species_baseline)
        genes = dic[species]
        genes_baseline = dic[species_baseline]
    else:
        dic = get_orthologs(genes, species_baseline, species)
        genes = dic[species]
        genes_baseline = dic[species_baseline]

    # NOTE: Now genes are synched via orthology

    # Get both dataframes
    dfs = [
        read_counts_from_file('celltype', genes, species),
        read_counts_from_file('celltype', genes_baseline, species_baseline),
        ]

    # Restrict to common cell types
    celltypes_shared = set(dfs[0].columns.values) & set(dfs[1].columns.values)
    celltypes0 = [ct for ct in dfs[0].columns if ct in celltypes_shared]
    celltypes1 = [ct for ct in dfs[1].columns if ct in celltypes_shared]
    dfs = [dfs[0].loc[:, celltypes0], dfs[1].loc[:, celltypes1]]

    return dfs
