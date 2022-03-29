from xml.dom import NotFoundErr
import h5py
import numpy as np
import pandas as pd
import json
import re
from scipy.cluster.hierarchy import linkage,leaves_list
from scipy.spatial.distance import pdist

# gives the name of dataset we want as an input
# celltype / celltype_dataset / celltype_dataset_timepoint
def read_file(df_type):
    h5_data = h5py.File('./static/scData/condensed_lung_atlas_in_cpm.h5',"r")
    
    df = pd.DataFrame(data=np.array(h5_data[df_type]\
    ['gene_expression_average']['block0_values']),\
    index=np.array(h5_data[df_type]['gene_expression_average']['axis1'])\
    ,columns=np.array(h5_data[df_type]['gene_expression_average']['axis0'])).T

    new_index = []
    for i in df.index:
        new_index.append(i.decode('utf-8'))
    # convert column name from binary to string
    new_column_name = []
    for i in df.columns:
        new_column_name.append(i.decode('utf-8'))
    df.index=new_index
    df.columns=new_column_name
    df = df.astype(np.float32)
    
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
        timepoints = set([name.split("_")[2] for name in dataset_name.columns])
        celltypes = set([name.split("_")[0] for name in dataset_name.columns])

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
    
