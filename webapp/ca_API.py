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
from models import (
        read_counts_from_file,
        dataset_by_timepoint,
        get_big_heatmap,
        get_friends,
        get_marker_genes,
    )
from validation import (
        validate_correct_genestr,
        validate_correct_celltypestr,
    )


class geneExp(Resource):
    '''API for the compressed atlas data, to be visualized as a heatmap'''
    def get(self):
        genestring = request.args.get("gene_names")
        plot_type = request.args.get("plot_type")
        data_type = request.args.get("data_type")

        gene_names = genestring.replace(' ', '').split(',')
        try:
            df = read_counts_from_file("celltype", genes=gene_names).T
        except KeyError:
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


class geneExpTimeUnified(Resource):
    def get(self):
        gene = request.args.get("gene")
        use_log = request.args.get("use_log") == '1'
        use_hierarchical = request.args.get("use_hierarchical") == '1'
        data = get_big_heatmap(gene, use_log, use_hierarchical)
        return data


class plotsForSeachGenes(Resource):
    def get(self):
        genestring = request.args.get("gene_names")
        gene_names = genestring.replace(' ', '').split(",")
        try:
            df = read_counts_from_file("celltype", genes=gene_names).T
        except KeyError:
            return None


        if len(gene_names) == 2:
            result = {}
            plot_df = df.filter(items=gene_names, axis=0)
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


class checkGenenames(Resource):
    def get(self):
        names = request.args.get("gene_names")
        names_validated = validate_correct_genestr(names)
        if names_validated is None:
            return {
                'outcome': 'fail',
                }
        else:
            return {
                'outcome': 'success',
                'genenames': names_validated,
                }


class markerGenes(Resource):
    def get(self):
        names = request.args.get("celltype_names")
        names_validated = validate_correct_celltypestr(names)
        if names_validated is None:
            return {
                'outcome': 'fail',
                }
        else:
            genenames = get_marker_genes(names_validated)
            return {
                'outcome': 'success',
                'genenames': genenames,
                }

