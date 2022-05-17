from xml.dom import NotFoundErr
import h5py
import numpy as np
import pandas as pd
import json
import re
import functools
from scipy.cluster.hierarchy import linkage,leaves_list
from scipy.spatial.distance import pdist

# gives the name of dataset we want as an input
# celltype / celltype_dataset / celltype_dataset_timepoint
def read_file(df_type):
    with h5py.File('./static/scData/condensed_lung_atlas_in_cpm.h5',"r") as h5_data:
        
        df = pd.DataFrame(
            data=np.array(h5_data[df_type]['gene_expression_average']['block0_values']).astype(np.float32),
            index=np.array(h5_data[df_type]['gene_expression_average']['axis1'].asstr()),
            columns=np.array(h5_data[df_type]['gene_expression_average']['axis0'].asstr()),
            ).T
    
    return df

# this function is used when more multiple gene names are searched by the user
# always run the read_file() to get dataframe in the right format before running this 
def data_preprocessing(input_gene_names,df_type):
    df = read_file(df_type)
    df_genes = df.index # all the genes that available in the current dataframe
    for search_gene in input_gene_names.split(","):
        if search_gene not in df_genes:
            return None

    # if user does not search for any specific genes, then plot everything
    if input_gene_names is None:
        plot_data = df.T
    else:
        a_gene_names = input_gene_names.split(",")
        plot_df = df.filter(items = a_gene_names,axis=0)
        plot_data = plot_df.T
    return plot_data

# this is a requirement if you re using this function as a sorting criteria
# return -1 if name_1 comes before name_2
# return 1 if name_2 comes before name_1
# return 0 if name_1 == name_2

# output should be in this order: E18.5, P1 P3 P7 P14 P21, 3m 18m 24m
def timepoint_reorder(tp1, tp2):
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
def dataset_by_timepoint(genename,df_type,datatype,plottype):
    # select and pre-preprocessing data
    df_tp = read_file(df_type)

    # select data for a given gene name
    df_filtered=df_tp.filter(items = [genename], axis=0)
    # for each dataset name, we construct a new dataframe for it (celltypes x timepoints)
    datasets = set([name.split("_")[1] for name in df_filtered.columns])
    dic_per_dataset = {}

    for i in datasets:
        dataset_name = df_filtered.filter(regex=i, axis=1)
        timepoints = list(set([name.split("_")[2] for name in dataset_name.columns]))
        celltypes = set([name.split("_")[0] for name in dataset_name.columns])
        
        # re-arrange the dataframe based on timepoint order:
        # e.g: P21, P3
        # sorted function depends on key(0,1,-1), output of the sorting function
        timepoints = sorted(timepoints, key=functools.cmp_to_key(timepoint_reorder))

        # create a new empty dataframe for each dataset   
        gene_exp_df = pd.DataFrame(np.eye(len(celltypes),len(timepoints)), columns=timepoints, index=celltypes)
        for tp in timepoints:
            for ct in celltypes:
                regex = '_'.join([ct,i,tp])
                # celltype Car4+ cappilaries can not be recognized, as it contains a + sign. which means one or more 4 in regex
                # replace regex '+' with string '\+'
                regex = re.sub('\+', '\+', regex)
                gene_expression = dataset_name.filter(regex=regex,axis=1)
                if gene_expression.empty:
                    gene_exp_df[tp][ct] = 0
                else:
                    gene_exp_df[tp][ct] = gene_expression.iloc[0,0] # same as getting the gene_expression value
        
        if datatype == "log10":
            gene_exp_df = np.log10(0.1+gene_exp_df)  
        if plottype == 'hieracical':
            distance = pdist(gene_exp_df.values)
            # print(distance)
            Z = linkage(distance,optimal_ordering=True)
            new_order = leaves_list(Z)
            gene_exp_df = gene_exp_df.iloc[new_order]
        dic_per_dataset[i] = json.loads(gene_exp_df.to_json())
    return dic_per_dataset
    # for each of the dataframe in dic_per_dataset, we convert it into json format
    

# def dataset_by_timepoint_dataset(genename,datatype,plottype):
''' generate a dictionary for the unified heatmap data'''
def dataset_unified(genename,datatype,plottype):
    # ,datatype,plottype
    df = read_file('celltype_dataset_timepoint')
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

    # create a dictionary of dictionary that contains dataset+timepoint as first key, then celltypes as the second keys
    '''
    {
    "ACZ_E18.5": {
        "Adventitial FB": -1.0,
        "Early adventitial FB": 0.020416259765625,
    }
    "TMS_P21:" {
        "Adventitial FB": 0.23456,
        "Early adventitial FB": -1,
    }
        '''
    expression = {}
    for dt in dt_combinations:
        expression[dt] = {}
        for ct in all_celltypes:
            name = "_".join([ct,dt])
            if name not in filtered_df.columns:
                exp_value = -1
            else:
                exp_value = float(filtered_df[name].values[0])
            expression[dt][ct] = exp_value
        
    if datatype == "log10":
        expression = np.log10(0.1+expression)
    
    if plottype == 'hieracical':
        # print(df.values.shape)  #(41x5)
        distance = pdist(expression.values)
        # print(distance)
        Z = linkage(distance,optimal_ordering=True)
        new_order = leaves_list(Z)
        expression = expression.iloc[new_order]
    
    dt_sorted = sorted([key for key in expression], key=functools.cmp_to_key(timepoint_reorder))

    
    return {'expression': expression, 'dataset_timepoint': dt_sorted, 'gene': genename}
    
