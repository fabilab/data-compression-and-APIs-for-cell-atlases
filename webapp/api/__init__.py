# Web imports
from flask import request, jsonify
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
        dataset_unified,
        get_friends,
        get_marker_genes,
        get_data_hyperoxia,
        get_celltype_abundances,
        get_data_differential,
        get_gene_ids,
        get_orthologs,
    )
from validation.genes import validate_correct_genestr
from validation.celltypes import validate_correct_celltypestr


class geneExp(Resource):
    '''API for the compressed atlas data, to be visualized as a heatmap'''
    def get(self):
        species = request.args.get("species")
        genestring = request.args.get("gene_names")
        gene_names = genestring.replace(' ', '').split(',')

        # If we are switching species, get orthologs
        new_species = request.args.get("newSpecies")
        if new_species is not None:
            gene_names = get_orthologs(
                gene_names, species, new_species,
            )
            species = new_species
            missing_genes = 'skip'
        else:
            missing_genes = 'throw'

        try:
            df = read_counts_from_file(
                    "celltype",
                    genes=gene_names,
                    species=species,
                    missing=missing_genes,
                    ).T
        except KeyError:
            return None

        # Just in case we skipped some
        gene_names = df.index.tolist()
        gene_ids = get_gene_ids(df.columns, species=species)

        dfl = np.log10(df + 0.5)

        celltypes_hierarchical = leaves_list(linkage(
            pdist(dfl.values),
            optimal_ordering=True),
        )
        if len(gene_names) <= 2:
            genes_hierarchical = list(range(len(gene_names)))
        else:
            genes_hierarchical = leaves_list(linkage(
                pdist(dfl.values.T),
                optimal_ordering=True),
            )

        result = {
            'data': df.to_dict(),
            'genes': df.columns.tolist(),
            'celltypes': df.index.tolist(),
            'celltypes_hierarchical': df.index[celltypes_hierarchical].tolist(),
            'genes_hierarchical': df.columns[genes_hierarchical].tolist(),
            'gene_ids': gene_ids,
            'species': species,
        }

        return jsonify(result)


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

        gene_id = get_gene_ids([gene])[gene]

        data = dataset_unified(gene)
        data['gene_id'] = gene_id
        return data


class geneExpHyperoxia(Resource):
    '''API for hyperoxia data'''
    def get(self):
        genestring = request.args.get("gene_names")

        # NOTE: probably do not want a full NLP-style parsing for the API,
        # but this might also be a little insufficient.
        genes = genestring.replace(' ', '').split(',')
        try:
            result = get_data_hyperoxia(genes=genes)
        except KeyError:
            return None

        for item in result:
            df = item['data']
            item['genes'] = df.columns.tolist()
            item['celltypes'] = df.index.tolist()

            df_delta = np.log2(df + 0.5) - np.log2(item['data_baseline'] + 0.5)
            new_order = leaves_list(linkage(
                pdist(df_delta.values),
                optimal_ordering=True,
                ))
            item['celltypes_hierarchical'] = df.index[new_order].tolist()
            if len(genes) > 1:
                new_order = leaves_list(linkage(
                    pdist(df_delta.values.T),
                    optimal_ordering=True,
                    ))
            else:
                new_order = [0]
            item['genes_hierarchical'] = df.columns[new_order].tolist()

            item['data'] = item['data'].to_dict()
            item['data_baseline'] = item['data_baseline'].to_dict()

        return jsonify(result)


class geneExpDifferential(Resource):
    '''API for generic differential expression'''
    def get(self):
        rqd = dict(request.args)
        conditions = [
                {'celltype': rqd['ct1'], 'timepoint': rqd['tp1'],
                 'dataset': rqd['ds1'], 'disease': rqd['dis1']},
                {'celltype': rqd['ct2'], 'timepoint': rqd['tp2'],
                 'dataset': rqd['ds2'], 'disease': rqd['dis2']},
        ]
        genes = rqd['genestr'].split(',')

        try:
            dfs = get_data_differential(conditions, genes=genes)
        except (KeyError, ValueError):
            return None

        # Get hierarchical clustering of cell types and genes
        df = dfs[0] - dfs[1]
        if len(genes) > 1:
            new_order = leaves_list(linkage(
                        pdist(df.values),
                        optimal_ordering=True,
                        ))
        else:
            new_order = [0]
        genes_hierarchical = df.index[new_order].tolist()
        new_order = leaves_list(linkage(
                    pdist(df.values.T),
                    optimal_ordering=True,
                    ))
        celltypes_hierarchical = df.columns[new_order].tolist()

        # Gene hyperlinks
        gene_ids = get_gene_ids(df.index)

        # Inject dfs into template
        heatmap_data = {
            'comparison': rqd['comparison'],
            'data': dfs[0].T.to_dict(),
            'data_baseline': dfs[1].T.to_dict(),
            'celltype': conditions[0]['celltype'],
            'celltype_baseline': conditions[1]['celltype'],
            'dataset': conditions[0]['dataset'],
            'dataset_baseline': conditions[1]['dataset'],
            'timepoint': conditions[0]['timepoint'],
            'timepoint_baseline': conditions[1]['timepoint'],
            'disease': conditions[0]['disease'],
            'disease_baseline': conditions[1]['disease'],
            'genes': dfs[0].index.tolist(),
            'celltypes': dfs[0].columns.tolist(),
            'genes_hierarchical': genes_hierarchical,
            'celltypes_hierarchical': celltypes_hierarchical,
            'gene_ids': gene_ids,
        }
        return heatmap_data


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


class celltypeAbundance(Resource):
    def get(self):
        timepoint = request.args.get("timepoint")
        kind = request.args.get("kind")

        return {
            'outcome': 'success',
            'celltypeabundance': get_celltype_abundances(timepoint, kind=kind),
            }

