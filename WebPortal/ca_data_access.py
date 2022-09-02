import h5py
import numpy as np
import pandas as pd
import json
import re
import functools
from scipy.cluster.hierarchy import linkage,leaves_list
from scipy.spatial.distance import pdist
import time

with h5py.File('./static/scData/condensed_lung_atlas_in_cpm.h5',"r") as h5_data:
    L =list(np.array(h5_data['celltype']['gene_expression_average']['axis0'].asstr()))


def gene_validation(all_genes,genes):
    """for each gene in the given list, check if it's valid

    :param all_genes: all genes in our dataset
    :type all_genes: list
    :param genes: input genes
    :type genes: list
    :return: invalid genes
    :rtype: list
    """
    invalid_genes = []
    for gene in genes:
        if gene not in all_genes:
            invalid_genes.append(gene)
    return invalid_genes

def read_file_average_exp(df_type,genes=None):
    '''
    read in a .h5 file from the source directory,
    return 
        1. if invalid input, return a list of invalid genes
        2. if all inputs are valid, return a dataframe based on the user's specified dataset type and genes of interest
        parameters:
            dataset type(string):  celltype/celltype_dataset/celltype_dataset_timepoint
            genes name (string)
        Return:
            success (boolean)
            a dataframe (df) / invalid genes (list)
    '''
    genes = genes.split(",")
    invalid_genes = gene_validation(L,genes)
    if len(invalid_genes) > 0:
        return False, invalid_genes
    with h5py.File('./static/scData/condensed_lung_atlas_in_cpm.h5',"r") as h5_data:
        # List of genes
        indexs = []
        if isinstance(genes,list):
            for gene in genes:
                indexs.append(L.index(gene))
        else:
            indexs.append(L.index(genes))
        # sort the indexs
        indexs.sort()
        new_genes = []
        for index in indexs:
            new_gene = L[index]
            new_genes.append(new_gene)

        data = np.array(h5_data[df_type]['gene_expression_average']['block0_values'][:, indexs]).astype(np.float32)
        columns=np.array(h5_data[df_type]['gene_expression_average']['axis1'].asstr())
        df = pd.DataFrame(
            data = data.T,
            index = new_genes,
            columns=columns,
        )

    return True, df

def read_file_proportion_exp(df_type,genes=None):
    '''
    read in a .h5 file from the source directory,
    return a dataframe based on the user's specified dataset type and genes of interest
        parameters:
            dataset type(string):  celltype/celltype_dataset/celltype_dataset_timepoint
            genes name (string)
        Return:
            a dataframe (df)
    '''
    genes = genes.split(",")
    invalid_genes = gene_validation(L,genes)
    if len(invalid_genes) > 0:
        return False, invalid_genes
    
    with h5py.File('./static/scData/condensed_lung_atlas_in_cpm.h5',"r") as h5_data:
        # List of genes
        indexs = []
        if isinstance(genes,list):
            for gene in genes:
                indexs.append(L.index(gene))
        else:
            indexs.append(L.index(genes))
        # sort the indexs
        indexs.sort()
        new_genes = []
        for index in indexs:
            new_gene = L[index]
            new_genes.append(new_gene)
        # expression table (only numbers,no index nor columns)
        data = np.array(h5_data[df_type]['gene_proportion_expression']['block0_values'][:, indexs]).astype(np.float32)
        columns=np.array(h5_data[df_type]['gene_proportion_expression']['axis1'].asstr())
        df = pd.DataFrame(
            data = data.T,
            index = new_genes,
            columns=columns,
        )
    return True, df

def get_marker_genes_list(celltype=None):
    '''
    get a list of marker genes of the given celltype
        Parameters:
            celltype (string)
        Returns:
            a list of gene names (list)
    '''
    with open('./static/scData/celltypeMarkerGeneList_descending.txt') as f:
        data = json.load(f)

    return data[celltype]

def marker_genes_expression(celltype):
    '''
    return a dataframe containing expression level of a list of genes in all celltype
    N rows(marker genes of the specified celltype) x M columns (all celltypes)
        parameter:
            celltype selected by the user (string)
        return:
            json object with genes as keys, value: list of celltype and gene expression values
    '''
    
    original_gene_list = get_marker_genes_list(celltype)
    gene_list = ','.join(original_gene_list)
    res, df_original = read_file_average_exp("celltype",gene_list)
    df_scaled = df_original.copy()
    df_scaled['current'] = df_scaled[celltype]
    for column in df_scaled.columns:
        df_scaled[column] = (df_scaled[column] / df_scaled['current']).round(3) * 100
    df_scaled.drop(['current'], axis=1,inplace=True)
    # also get the proportion data
    _, df_proportion = read_file_proportion_exp("celltype",gene_list)

    return {'data_original': json.loads(df_original.to_json()),'data_scaled':json.loads(df_scaled.to_json()), 'exp_proportion':json.loads(df_proportion.to_json()),'order': original_gene_list}

def timepoint_reorder(tp1, tp2):
    '''
    sort the order of two given timepoints, such that timepoints are in order:
    --> E18.5, P1 P3 P7 P14 P21, 3m 18m 24m -->
        parameters:
            2 timepoints (string)
        return:
            0 if tp_1 == tp_2
            -1 if tp_1 comes before tp_2
             1 if tp_2 comes before tp_1
    '''
    if '_' in tp1:
        tp1 = tp1.split('_')[1]
    
    if '_' in tp2:
        tp2 = tp2.split('_')[1]

    if tp1 == tp2:
        return 0
    
    if 'P' in tp1 and 'P' in tp2:
        value_1 = int(tp1.split('P')[1])
        value_2 = int(tp2.split('P')[1])
        if value_1 < value_2:
            return 1
        else:
            return -1
    
    if 'm' in tp1 and 'm' in tp2:
        value_1 = int(tp1.split('m')[0])
        value_2 = int(tp2.split('m')[0])
        if value_1 < value_2:
            return 1
        else:
            return -1
    
    # E always before others
    if 'E' in tp1 and 'E' not in tp2:
        return 1
    
    if 'E' in tp2 and 'E' not in tp1:
        return -1
    
    # m always after others
    if 'm' in tp1 and 'm' not in tp2:
        return -1

    if 'm' in tp2 and 'm' not in tp1:
        return 1 


#  function that get the cell type dataset timepoint as a dictionary
# user input a gene name.
#  for each unique dataset of this gene, we plot a heatmap of all timepoint vs celltypes 
def dataset_pre_process(df,gene):
# select and pre-preprocessing data

    # select data for a given gene name
    df_filtered=df.loc[[gene]]
    # for each dataset name, we construct a new dataframe for it (celltypes x timepoints)
    datasets = set([name.split("_")[1] for name in df_filtered.columns])
    dic_per_dataset = {}
    new_cell_type_order = {}
    # ACZ,Hur.2021,TMS
    for i in datasets:
        columns_with_this_dataset = []
        for col_name in df_filtered.columns:
             if i in col_name:
                 columns_with_this_dataset.append(col_name)
        dataset_name = df_filtered[columns_with_this_dataset] 
        timepoints = list(set([name.split("_")[2] for name in dataset_name.columns]))
        celltypes = set([name.split("_")[0] for name in dataset_name.columns])
        
        # timepoint sorted function depends on key(0,1,-1), output of the sorting function
        timepoints = sorted(timepoints, key=functools.cmp_to_key(timepoint_reorder))

        # create a new empty dataframe for each dataset   
        gene_exp_df = pd.DataFrame(np.eye(len(celltypes),len(timepoints)), columns=timepoints, index=celltypes)
        for tp in timepoints:
            for ct in celltypes:
                col_name = '_'.join([ct,i,tp])
                # celltype Car4+ cappilaries can not be recognized, as it contains a + sign. which means one or more 4 in regex
                # replace regex '+' with string '\+'
                col_name = re.sub('\+', '\+', col_name)
                try:
                    # try to see if the col_name is in df
                    gene_exp_df[tp][ct] = dataset_name[col_name].values[0]
                except:
                    # col_name not in df
                    pass

        # celltype order for hierarchical clustered set
        distance = pdist(gene_exp_df.values)
        Z = linkage(distance,optimal_ordering=True)
        new_order = leaves_list(Z)
        clustered_gene_exp_df = gene_exp_df.iloc[new_order]
        new_cell_type_order[i] = [ct for ct in clustered_gene_exp_df.index]
        dic_per_dataset[i] = json.loads(gene_exp_df.to_json())
    return {"result":dic_per_dataset, "hierarchicalCelltypeOrder":new_cell_type_order}
    # for each of the dataframe in dic_per_dataset, we convert it into json format

def result_datasets(df_type,gene):
    res, df_avg = read_file_average_exp(df_type,gene)
    res, df_pro = read_file_proportion_exp(df_type,gene)

    exp_avg = dataset_pre_process(df_avg,gene)
    exp_pro = dataset_pre_process(df_pro,gene)

    return {"exp_avg":exp_avg['result'], "exp_pro":exp_pro['result'], "hierarchicalCelltypeOrder":exp_avg['hierarchicalCelltypeOrder']}


    
# giving a dataframe, extract the information needed for generate the corresponded plot
def unified_process(df,genename):

    filtered_df = df.filter(items=[genename],axis=0)
    all_celltypes = []
    dt_combinations = []  # store all the existing dataset and timepoint combinations

    for column_name in filtered_df.columns:
        celltype = column_name.split("_")[0]
        dataset_timepoint = column_name.split(celltype+"_")[1]
        if celltype not in all_celltypes:
            all_celltypes.append(celltype)
        if dataset_timepoint not in dt_combinations:
            dt_combinations.append(dataset_timepoint) 

    gene_exp_df = pd.DataFrame(np.eye(len(all_celltypes),len(dt_combinations)), columns=dt_combinations, index=all_celltypes)
    for dt in dt_combinations:
        for ct in all_celltypes:
            name = "_".join([ct,dt])
            if name not in filtered_df.columns:
                exp_value = -1
            else:
                exp_value = float(filtered_df[name].values[0])

            gene_exp_df[dt][ct] = exp_value
    
    cell_type_label = all_celltypes
    distance = pdist(gene_exp_df.values)
    Z = linkage(distance,optimal_ordering=True)
    new_order = leaves_list(Z)
    clustered_gene_exp_df = gene_exp_df.iloc[new_order]
    new_cell_type_order = [ct for ct in clustered_gene_exp_df.index]
    dt_sorted = sorted([key for key in gene_exp_df.columns], key=functools.cmp_to_key(timepoint_reorder))

    return {'exp': json.loads(gene_exp_df.to_json()), 'dataset_timepoint': dt_sorted, 'cell_type': cell_type_label, 'gene': genename, 'hierarchicalCelltypeOrder':new_cell_type_order}

def dataset_unified(genename):
    
    res, df_avg = read_file_average_exp('celltype_dataset_timepoint',genename)
    res, df_pro = read_file_proportion_exp('celltype_dataset_timepoint',genename)

    exp_avg = unified_process(df_avg,genename)
    exp_pro = unified_process(df_pro,genename)

    return {'exp_avg': exp_avg['exp'] , 'exp_pro': exp_pro['exp'],'dataset_timepoint': exp_avg['dataset_timepoint'], 'cell_type': exp_avg['cell_type'], 'gene': exp_avg['gene'], 'hierarchicalCelltypeOrder':exp_avg['hierarchicalCelltypeOrder']}

def unified_by_cell_process(df,celltype,genes):

    df_filtered = df.filter(regex=celltype)
    dataset_timepoint = []
    for column in df_filtered.columns:
        dt = column.split('_',1)[1]
        dataset_timepoint.append(dt)
    df_filtered.columns = dataset_timepoint
    dt_sorted = sorted([key for key in df_filtered.columns], key=functools.cmp_to_key(timepoint_reorder))

    distance = pdist(df_filtered.values)
    Z = linkage(distance,optimal_ordering=True)
    new_order = leaves_list(Z)
    clustered_gene_exp_df = df_filtered.iloc[new_order]
    new_gene_order = [ct for ct in clustered_gene_exp_df.index]

    return {'exp': json.loads(df_filtered.to_json()), 'dataset_timepoint': dt_sorted, 'cell_type': celltype, 'genes': genes.split(','), 'clustered_gene_order': new_gene_order}

def result_unified_by_cell(celltype,genes):
    
    res, df_avg = read_file_average_exp('celltype_dataset_timepoint',genes)
    res, df_pro = read_file_proportion_exp('celltype_dataset_timepoint',genes)

    exp_avg = unified_by_cell_process(df_avg,celltype,genes)
    exp_pro = unified_by_cell_process(df_pro,celltype,genes)

    return {'exp_avg': exp_avg['exp'] , 'exp_pro': exp_pro['exp'],'dataset_timepoint': exp_avg['dataset_timepoint'], 'cell_type': exp_avg['cell_type'], 'genes': exp_avg['genes'], 'clustered_gene_order': exp_avg['clustered_gene_order']}


