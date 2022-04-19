# Web imports
from flask import request
from flask_restful import Resource, Api
import json

# Data import
import h5py
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
import plotly


# Helper functions
from helper import (
        data_preprocessing,
        dataset_by_timepoint,
        get_big_heatmap,
        get_friends,
    )


class geneExpTimeUnified(Resource):
    def get(self):
        gene = request.args.get("gene")
        use_log = request.args.get("use_log") == '1'
        use_hierarchical = request.args.get("use_hierarchical") == '1'
        data = get_big_heatmap(gene, use_log, use_hierarchical)
        return data


class geneExpTime(Resource):
    def get(self):
        genename = request.args.get("gene")
        datatype = request.args.get("datatype")
        plottype = request.args.get("plottype")

        data = dataset_by_timepoint(
            genename,
            "celltype_dataset_timepoint",
            datatype,
            plottype,
        )
        return data


class geneExp(Resource):
    '''API for the compressed atlas data, to be visualized as a heatmap'''
    def get(self):
        gene_names = request.args.get("gene_names")
        plot_type = request.args.get("plot_type")
        data_type = request.args.get("data_type")
        df = data_preprocessing(gene_names, "celltype")
        if df is None:
            return None

        if data_type == "log10":
            df = np.log10(0.1 + df)

        if plot_type == "hierachical":
            print(df.values.shape)  # (41x5)
            distance = pdist(df.values)
            # print(distance)
            Z = linkage(distance, optimal_ordering=True)
            new_order = leaves_list(Z)
            df = df.iloc[new_order]

        return json.loads(df.to_json())


class plotsForSeachGenes(Resource):
    def get(self):
        gene_names = request.args.get("gene_names")
        df = data_preprocessing(gene_names, "celltype")

        if df is None:
            return None

        df = df.T
        a_gene_names = gene_names.split(",")

        if len(a_gene_names) == 2:
            result = {}
            plot_df = df.filter(items=a_gene_names, axis=0)
            gene1 = plot_df.index[0]
            gene2 = plot_df.index[1]

            gene1_expr = list(plot_df.loc[gene1])
            gene2_expr = list(plot_df.loc[gene2])
            # plot_data = plot_df.to_json()
            result["gene1_name"] = gene1
            result["gene2_name"] = gene2
            result["gene1_expr"] = gene1_expr
            result["gene2_expr"] = gene2_expr
            result["cell_types"] = list(plot_df.columns)

        return result


class geneFriends(Resource):
    def get(self):
        genenames = request.args.get("gene_names")
        genenames = genenames.replace(" ", "").split(",")
        data = get_friends(genenames)
        return data
