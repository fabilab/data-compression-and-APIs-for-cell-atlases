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
    

def dataset_by_timepoint_dataset(genename,datatype,plottype):
    df = read_file('celltype_dataset_timepoint')
    filtered_df = df.filter(items=[genename],axis=0)

    # get a list of all the celltypes
    all_celltypes = []
    all_timepoints = []
    all_datasets = []
    for column_name in filtered_df.columns:
        # print(column_name)
        all_celltypes.append(column_name.split("_")[0])
        all_datasets.append(column_name.split("_")[1])
        all_timepoints.append(column_name.split("_")[2])

    # unique
    all_celltypes = list(set(all_celltypes))
    all_datasets = list(set(all_datasets))
    all_timepoints = list(set(all_timepoints))
    all_timepoints = sorted(all_timepoints, key=functools.cmp_to_key(timepoint_reorder),reverse=True)
    # to construct this for each cell, e.g
    # T cell = [
    #     TMS->[value1,value2 .....value(length of all timepoints)],
    #     ACZ->[value1,value2 .....value(length of all timepoints)],
    #     Hurskainen2021->[value1,value2 .....value(length of all timepoints)] 
    # ]
    bigData = []
    for celltype in all_celltypes:
        list_all_dataset = []
        for dataset in all_datasets:
            list_all_times = []
            for timepoint in all_timepoints:
                names = (celltype,dataset,timepoint)
                # check if the combination of celltypes_dataset_timepoint exists
                name = "_".join(names)
                if name not in filtered_df.columns:
                    list_all_times.append(np.nan)
                    continue
                else:
                # check if the value is None
                    has_value = filtered_df[name].values[0]
                    if has_value is None:
                        list_all_times.append(0)
                    else:
                        list_all_times.append(has_value)
            list_all_dataset.append(list_all_times)
        bigData += list_all_dataset
    
    all_celltypes_repeated_3 = list(np.repeat(all_celltypes, 3)) # need to repeate for all 3 dataset
    big_df = pd.DataFrame(data=bigData,columns=all_timepoints)
    big_df['celltype'] = all_celltypes_repeated_3
    big_df['dataset'] = all_datasets * len(all_celltypes)

    big_df.index = [x + '_' + y for x, y in zip(big_df['celltype'], big_df['dataset'])]
    df = big_df[all_timepoints].copy().T

    if datatype == "log10":
        df = np.log10(0.1+df)
    
    if plottype == 'hieracical':
        # print(df.values.shape)  #(41x5)
        distance = pdist(df.values)
        # print(distance)
        Z = linkage(distance,optimal_ordering=True)
        new_order = leaves_list(Z)
        df = df.iloc[new_order]
    
    return json.loads(df.to_json())
