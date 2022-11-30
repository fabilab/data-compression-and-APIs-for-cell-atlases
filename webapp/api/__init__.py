# Web imports
from flask import request, jsonify, abort
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
        get_data_overtime_1gene,
        get_data_overtime_1celltype,
        get_data_hyperoxia,
        get_friends,
        get_marker_genes,
        get_celltype_abundances,
        get_data_differential,
        get_data_species_comparison,
        get_gene_ids,
        get_gene_ontology_terms,
        get_genes_in_GO_term,
        get_orthologs,
    )
from validation.genes import validate_correct_genestr
from validation.celltypes import validate_correct_celltypestr


class expressionByCelltype(Resource):
    '''API for the compressed atlas data, to be visualized as a heatmap'''
    def post(self):
        '''No data is actually posted, but this is uncached and can be longer

        In particular, when many genes are requested at once, we run into the
        length limitation of GET requests. The downside of using POST is that
        there's no caching, bookmarking, and such.
        '''
        args = request.form
        return self.get(args=args)

    def get(self, args=None):
        if args is None:
            args = request.args
        species = args.get("species")
        genestring = args.get("gene_names")

        # A cap on gene names to avoid overload is reasonable
        genestring = ','.join(genestring.replace(' ', '').split(',')[:500])

        genestring = validate_correct_genestr(
                genestring, species=species, missing='skip')

        if genestring is None or genestring == '':
            return None
        gene_names = genestring.split(',')

        # If we are switching species, get orthologs
        new_species = args.get("newSpecies")
        if new_species is not None:
            gene_names = get_orthologs(
                gene_names, species, new_species,
            )[new_species]
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
                    )

            df_fractions = read_counts_from_file(
                    "celltype",
                    genes=gene_names,
                    species=species,
                    missing=missing_genes,
                    key='gene_proportion_expression',
                    )

        except KeyError:
            return None

        # Just in case we skipped some
        gene_names = df.index.tolist()
        gene_ids = get_gene_ids(df.index, species=species)

        go_terms = get_gene_ontology_terms(gene_names, species=species)

        dfl = np.log10(df + 0.5)

        idx_ct_hierarchical = leaves_list(linkage(
            pdist(dfl.T.values),
            optimal_ordering=True),
        )
        idx_ct_hierarchical = [int(x) for x in idx_ct_hierarchical]

        if len(gene_names) <= 2:
            idx_genes_hierarchical = list(range(len(gene_names)))
        else:
            idx_genes_hierarchical = leaves_list(linkage(
                pdist(dfl.values),
                optimal_ordering=True),
            )
            idx_genes_hierarchical = [int(x) for x in idx_genes_hierarchical]

        result = {
            'data': df.values.tolist(),
            'data_fractions': df_fractions.values.tolist(),
            'genes': df.index.tolist(),
            'celltypes': df.columns.tolist(),
            'celltypes_hierarchical': idx_ct_hierarchical,
            'genes_hierarchical': idx_genes_hierarchical,
            'gene_ids': gene_ids,
            'GO_terms': go_terms,
            'species': species,
        }
        return result


class expressionOvertime1Gene(Resource):
    def get(self):
        gene_name = request.args.get("gene")
        species = request.args.get("species")

        # If we are switching species, get orthologs
        new_species = request.args.get("newSpecies")
        if new_species is not None:
            gene_name = get_orthologs(
                [gene_name], species, new_species,
            )[new_species]

            if len(gene_name) == 0:
                return None
            gene_name = gene_name[0]

            species = new_species
            missing_genes = 'skip'

        else:
            missing_genes = 'throw'

        gene_id = get_gene_ids([gene_name], species=species)[gene_name]

        data = get_data_overtime_1gene(gene_name, species=species)
        data['gene_id'] = gene_id

        similar_genes = get_friends([gene_name], species=species).split(',')
        if similar_genes[0] == gene_name:
            similar_genes = similar_genes[1:]
        data['similarGenes'] = similar_genes

        return data


class expressionOvertime1Celltype(Resource):
    def get(self):
        species = request.args.get("species")

        celltype = request.args.get("celltype")
        celltype_validated = validate_correct_celltypestr(celltype)

        genestring = request.args.get("gene_names")

        # A cap on gene names to avoid overload is reasonable
        genestring = ','.join(genestring.replace(' ', '').split(',')[:500])

        genestring = validate_correct_genestr(
                genestring, species=species, missing='skip')

        if genestring is None or genestring == '':
            return None
        gene_names = genestring.split(',')

        # If we are switching species, get orthologs
        new_species = request.args.get("newSpecies")
        if new_species is not None:
            gene_names = get_orthologs(
                gene_names, species, new_species,
            )[new_species]

            if len(gene_names) == 0:
                return None

            species = new_species
            missing_genes = 'skip'

        else:
            missing_genes = 'throw'

        gene_ids = get_gene_ids(gene_names, species=species)

        data = get_data_overtime_1celltype(celltype_validated, gene_names, species=species)
        data['gene_ids'] = gene_ids
        data['celltype'] = celltype

        similar_genes = get_friends(gene_names, species=species).split(',')
        # Exclude from similar genes the ones you have already
        similar_genes = [g for g in similar_genes if g not in gene_names]

        data['similarGenes'] = similar_genes

        return data


class expressionHyperoxia(Resource):
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


class geneExpSpeciesComparison(Resource):
    '''Comparison between species'''
    def get(self):
        species = request.args.get('species')
        species_baseline = request.args.get('species_baseline')
        genes = request.args.get('genes').split(',')

        # Get the counts
        # NOTE: this function restricts to the intersection of cell types,
        # which makes the hierarchical clustering easy. In summary, both
        # genes and cell types are fully synched now
        dfs = get_data_species_comparison(species, species_baseline, genes)

        # Hierarchical clustering
        df = np.log10(dfs[0] + 0.5)

        # Get hierarchical clustering of genes
        if len(genes) > 2:
            new_order = leaves_list(linkage(
                        pdist(df.values),
                        optimal_ordering=True,
                        ))
            genes_hierarchical = df.index[new_order].tolist()
            genes_hierarchical_baseline = dfs[1].index[new_order].tolist()
        else:
            genes_hierarchical = dfs[0].index.tolist()
            genes_hierarchical_baseline = dfs[1].index.tolist()

        # Get hierarchical clustering of cell types
        # NOTE: both dfs have the same celltypes (see above note)
        new_order = leaves_list(linkage(
                    pdist(df.values.T),
                    optimal_ordering=True,
                    ))
        celltypes_hierarchical = df.columns[new_order].tolist()

        # Gene hyperlinks (they hold for both)
        gene_ids = get_gene_ids(df.index, species)

        # Inject dfs into template
        # NOTE: the whole converting DataFrame to dict of dict makes this quite
        # a bit more heavy than it should be... just use a list of lists and
        # accompanying lists of indices
        heatmap_data = {
            'data': dfs[0].T.to_dict(),
            'data_baseline': dfs[1].T.to_dict(),
            'genes': dfs[0].index.tolist(),
            'genes_baseline': dfs[1].index.tolist(),
            'celltypes': dfs[0].columns.tolist(),
            'celltypes_baseline': dfs[1].columns.tolist(),
            'genes_hierarchical': genes_hierarchical,
            'celltypes_hierarchical': celltypes_hierarchical,
            'genes_hierarchical_baseline': genes_hierarchical_baseline,
            'gene_ids': gene_ids,
            'species': species,
            'species_baseline': species_baseline,
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


class genesInGOTerm(Resource):
    def get(self):
        go_term = request.args.get("goTerm")
        species = request.args.get("species")
        genes = get_genes_in_GO_term(go_term, species)
        return ','.join(genes)


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

