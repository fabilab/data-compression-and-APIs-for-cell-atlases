from xml.dom import NotFoundErr
import h5py
import numpy as np
import pandas as pd
import json
from scipy.cluster.hierarchy import linkage,leaves_list
from scipy.spatial.distance import pdist

def read_file():
    h5_data = h5py.File('./static/scData/condensed_lung_atlas.h5',"r")

    # Dataset 1:
    df = pd.DataFrame(data=np.array(h5_data['cell_type']\
    ['gene_expression_average']['block0_values']),\
    index=np.array(h5_data['cell_type']['gene_expression_average']['axis1'])\
    ,columns=np.array(h5_data['cell_type']['gene_expression_average']['axis0'])).T

    # Dataset 2:
    df_2 = pd.DataFrame(data=np.array(h5_data['cell_type_dataset']\
    ['gene_expression_average']['block0_values']),\
    index=np.array(h5_data['cell_type_dataset']['gene_expression_average']['axis1'])\
    ,columns=np.array(h5_data['cell_type_dataset']['gene_expression_average']['axis0'])).T

    # Dataset 3:
    df_3 = pd.DataFrame(data=np.array(h5_data['cell_type_dataset_timepoint']\
    ['gene_expression_average']['block0_values']),\
    index=np.array(h5_data['cell_type_dataset_timepoint']['gene_expression_average']['axis1'])\
    ,columns=np.array(h5_data['cell_type_dataset_timepoint']['gene_expression_average']['axis0'])).T

    return [df,df_2,df_3]

def data_preprocessing(gene_names=None):
    data = read_file()
    plot_data_sets = []
    for df in data:
    # current index in the dataframe is writtern as binary string.
    # We need to convert it into normal string
        new_index=[]
        for i in df.index:
            new_index.append(i.decode('utf-8'))
        # Similarly for columns name
        new_columns=[]
        for i in df.columns:
            new_columns.append(i.decode('utf-8'))

        df.index = new_index
        df.columns = new_columns

        all_genes = df.index
        for user_gene in gene_names.split(","):
            if user_gene not in all_genes:
                return None

        if gene_names is None:
            plot_data = df.T
        else:
            a_gene_names = gene_names.split(",")
            plot_df = df.filter(items = a_gene_names,axis=0)
            plot_data = plot_df.T
        plot_data_sets.append(plot_data)
    return plot_data_sets

#  function that get the cell type dataset timepoint as a dictionary
# user input a gene name.
#  for each unique dataset of this gene, we plot a heatmap of all timepoint vs celltypes 
def dataset_by_timepoint(genename,datatype,plottype):
    # select and pre-preprocessing data
    df_tp = read_file()[2]
    new_index = []
    for i in df_tp.index:
        new_index.append(i.decode('utf-8'))
    new_column_name = []
    for i in df_tp.columns:
        new_column_name.append(i.decode('utf-8'))
    df_tp.index=new_index
    df_tp.columns=new_column_name

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
                gene_expression = dataset_name.filter(regex=regex,axis=1)
                if gene_expression.empty:
                    gene_expression = 0
                else:
                    gene_exp_df[tp][ct] = gene_expression = gene_expression.iloc[0,0] # same as getting the gene_expression value
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
    

# def dataset_by_timepoint():
