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
    df = read_file()[2]
    # convert index from binary to string
    new_index = []
    for i in df.index:
        new_index.append(i.decode('utf-8'))
    # convert column name from binary to string
    new_column_name = []
    for i in df.columns:
        new_column_name.append(i.decode('utf-8'))
    df.index=new_index
    df.columns=new_column_name

    df_filtered=df.filter(items = [genename], axis=0)
    # modify data base on the selected datatype and plottype (cpm?log10?hierachical?)
    if datatype == "log10":
        df_filtered = np.log10(0.1+df_filtered)
    
    if plottype == 'hieracical':
        print(df_filtered.values.shape) #df_filtered.values.shape should have (331x3)
        # distance = pdist(df_filtered.values)
        #  print(distance)
        # Z = linkage(distance,optimal_ordering=True)
        # new_order = leaves_list(Z)
        # df_filtered = df_filtered.iloc[new_order]

    datasets = set([name.split("_")[1] for name in df_filtered.columns])
    result = {}
    # for each dataset, we give it a unique heatmap:
    # ACZ, Hurskainen2021, TMS ...
    for dataset in datasets:
        result[dataset] = {}
        current_data = df_filtered.filter(regex=dataset, axis=1)
        celltypes = set([name.split("_")[0] for name in current_data.columns])
        timepoints = set([name.split("_")[2] for name in current_data.columns])
        for timepoint in timepoints:
            timepoint_data = current_data.filter(regex=timepoint, axis=1)
            result[dataset][timepoint] = {}
            for cell in celltypes:
                celltype_data = timepoint_data.filter(regex="^{}_".format(cell), axis=1)
                # if the cell is missing tin thsi time point
                if len(celltype_data.values[0]) == 0:
                    result[dataset][timepoint][cell] = 0
                else:
                    result[dataset][timepoint][cell] = float(celltype_data.values[0][0])

    return result