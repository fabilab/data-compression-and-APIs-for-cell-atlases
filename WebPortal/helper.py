from xml.dom import NotFoundErr
import h5py
import numpy as np
import pandas as pd

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

