# vim: fdm=indent
'''
author:     Fabio Zanini
date:       22/05/22
content:    Main flask app for compressed atlases
'''
from base64 import decode
import json
from flask import (
    Flask,
    send_from_directory,
    request,
    redirect,
    url_for,
    render_template,
    abort,
)
from flask_restful import Api
#from flask_cors import CORS
import numpy as np

from api import (
    geneExp,
    plotsForSeachGenes,
    geneFriends,
    genesInGOTerm,
    geneExpTimeUnified,
    geneExpHyperoxia,
    geneExpDifferential,
    geneExpSpeciesComparison,
    checkGenenames,
    markerGenes,
    celltypeAbundance,
)
from models import (
        get_celltype_abundances,
        get_data_differential,
        get_gene_ids,
        get_friends,
        get_data_species_comparison,
        get_orthologs,
        get_genes_in_GO_term,
        get_gsea,
        get_kegg_urls,
)
from validation.genes import validate_correct_genestr
from validation.timepoints import validate_correct_timepoint
from voice_recognition import mod as voice_control_blueprint
from text_recognition import mod as text_control_blueprint


##############################
app = Flask(__name__, static_url_path="/static", template_folder="templates")
app_api = Api(app)
# Note: this might be unsafe
#CORS(app)
with open('secret_key.txt') as f:
    app.config['SECRET_KEY'] = f.read()
##############################


##############################
# Views
##############################
@app.route("/")
def index():
    return redirect(url_for('text_control'))


# Control pages
@app.route("/text_control", methods=["GET"])
def text_control():
    """A single text bar to ask questions or post commands"""
    return render_template(
        "text_control.html",
        )

@app.route("/voice_control", methods=["GET"])
def voice_control():
    """Name says it all"""
    return render_template(
            "voice_control.html",
            )


@app.route("/heatmap_by_celltype", methods=['GET'])
def heatmap_by_celltype():
    species = request.args.get('species')
    if species is None:
        species = 'mouse'
    genestring = request.args.get("genestring")
    if genestring is None:
        pathway = request.args.get("pathway")
        if pathway is not None:
            if ' (GO' in pathway:
                pathway = pathway[:pathway.find(' (GO')]
            try:
                genes = get_genes_in_GO_term(pathway, species=species)
            except KeyError:
                return abort(404)

        else:
            genes = [
                'Col1a1,Col2a1',
                'Adh1,Col13a1,Col14a1',
                'Tgfbi,Pdgfra,Crh,Hhip,Pdgfrb',
                'Pecam1,Gja5,Vwf,Car8,Car4',
                'Ptprc,Cd19,Gzma,Cd3d,Cd68',
                'Epcam',
                ]
            if species in ('human', 'lemur'):
                genes = get_orthologs(genes, 'mouse', species)
        genestring = ','.join(genes)
    searchstring = genestring.replace(" ", "")
    return render_template(
            "heatmap_by_celltype.html",
            searchstring=searchstring,
            species=species,
            )


@app.route("/heatmap_development", methods=['GET'])
def heatmap_development():
    species = request.args.get('species')
    if species is None:
        species = 'mouse'
    genestring = request.args.get("genestring")
    if genestring is None:
        genestring = 'Car4'
    searchstring = genestring.replace(" ", "")

    similar_genes = get_friends([searchstring]).split(',')
    if similar_genes[0] == searchstring:
        similar_genes = similar_genes[1:]

    return render_template(
            "heatmap_unified.html",
            searchstring=searchstring,
            species=species,
            similarGenes=similar_genes,
            )


@app.route("/heatmap_hyperoxia", methods=["GET"])
def heatmap_hyperoxia():
    """A sort of heatmap with hyperoxia"""
    genestring = request.args.get("genestring")
    if genestring is None:
        genestring = ','.join([
            'Col1a1,Col2a1',
            'Adh1,Col13a1,Col14a1',
            'Tgfbi,Pdgfra,Crh,Hhip,Pdgfrb',
            'Pecam1,Gja5,Vwf,Car8,Car4',
            'Ptprc,Cd19,Gzma,Cd3d,Cd68',
            'Epcam',
            ])
    searchstring = genestring.replace(" ", "")
    # Default dataset/timepoints combos
    dataset_timepoints = [
        'ACZ_P7', 'Hurskainen2021_P3', 'Hurskainen2021_P7', 'Hurskainen2021_P14',
    ]
    return render_template(
            "heatmap_hyperoxia.html",
            searchstring=searchstring,
            datasetTimepoints=dataset_timepoints,
            )


# NOTE: This is where API and views break down and react would be better
@app.route("/heatmap_differential", methods=["GET"])
def heatmap_differential_genes():
    """A sort of heatmap for differential gene expression"""
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import pdist

    rqd = dict(request.args)
    conditions = [
            {'celltype': rqd['ct1'], 'timepoint': rqd['tp1'],
             'dataset': rqd['ds1'], 'disease': rqd['dis1']},
            {'celltype': rqd['ct2'], 'timepoint': rqd['tp2'],
             'dataset': rqd['ds2'], 'disease': rqd['dis2']},
    ]
    dfs = get_data_differential(
            conditions,
            kind=rqd['kind'],
            n_genes=int(rqd['n_genes']),
    )

    # Get hierarchical clustering of cell types and genes
    # FIXME: this works poorly, trying out HC on the log
    df = np.log10(dfs[0] + 0.5)
    new_order = leaves_list(linkage(
                pdist(df.values),
                optimal_ordering=True,
                ))
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

    # Set search string
    searchstring = ','.join(dfs[0].index)
    return render_template(
            "heatmap_differential.html",
            searchstring=searchstring,
            heatmapData=heatmap_data,
            )


@app.route("/list_celltypes/<timepoint>", methods=["GET"])
def list_celltypes_timepoint(timepoint):
    '''List cell types and their abundances'''
    timepoint = validate_correct_timepoint(timepoint)

    celltype_dict = get_celltype_abundances(
            timepoint,
            kind='qualitative',
            )
    return render_template(
            'list_celltypes.html',
            timepoint=timepoint,
            celltypes=celltype_dict,
            searchstring=timepoint,
            )


@app.route("/celltype_abundance/<timepoint>", methods=["GET"])
def plot_celltype_abundance(timepoint):
    '''Plot cell type abundances'''
    celltype_dict = get_celltype_abundances(
            timepoint,
            kind='quantitative',
            )
    return render_template(
            'celltype_abundance.html',
            timepoint=timepoint,
            celltypes=celltype_dict,
            searchstring=timepoint,
            )


@app.route("/barplot_gsea", methods=["GET", "POST"])
def plot_barplot_GSEA():
    '''Barplot for gene set enrichment analysis'''
    if request.method == "POST":
        args = request.form
    else:
        args = request.args

    genestring = args.get('genes')
    species = args.get('species')
    gene_set = args.get('gene_set')

    genes = validate_correct_genestr(genestring, species=species).split(',')
    if gene_set is None:
        data = get_gsea(genes, species)
    else:
        data = get_gsea(genes, species, gene_set=gene_set)

    if 'KEGG' in gene_set:
        pathway_urls = get_kegg_urls(data.index)
    else:
        pathway_urls = []

    # Cut too long results
    data = data.iloc[:15]

    return render_template(
        'barplot_gsea.html',
        species=species,
        plotData=dict(
            pathways=data.index.tolist(),
            pathways_urls=pathway_urls,
            overlap=data['Overlap'].values.tolist(),
            neglog10_p_value=(-np.log10(data['Adjusted P-value'].values)).tolist(),
            ),
        )


@app.route("/heatmap_species_comparison", methods=["GET"])
def heatmap_species_comparison():
    '''Plot heatmap of cross-species comparison'''
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import pdist

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

    # Set search string
    searchstring = ','.join(dfs[0].index)

    return render_template(
        'heatmap_species_comparison.html',
        species=species,
        heatmapData=heatmap_data,
        searchstring=searchstring,
    )


# Static assets (JS/CSS)
@app.route("/js/<path:path>")
def send_js(path):
    """JavaScript assets"""
    return send_from_directory("static/js", path)


@app.route("/css/<path:path>")
def send_css(path):
    """CSS stylesheets"""
    return send_from_directory("static/css", path)


@app.route("/favicon.ico")
def favicon():
    return send_from_directory("static", "favicon.ico")


# Blueprints
app.register_blueprint(text_control_blueprint)
app.register_blueprint(voice_control_blueprint)


# API endpoints
app_api.add_resource(geneExp, "/data/by_celltype")
app_api.add_resource(geneFriends, "/data/gene_friends")
app_api.add_resource(genesInGOTerm, "/data/genes_in_go_term")
app_api.add_resource(geneExpTimeUnified, "/data_heatmap_unified")
app_api.add_resource(geneExpHyperoxia, "/data/hyperoxia")
app_api.add_resource(geneExpDifferential, "/data/differential")
app_api.add_resource(geneExpSpeciesComparison, "/data/speciescomparison")
app_api.add_resource(checkGenenames, "/check_genenames")
app_api.add_resource(markerGenes, "/data/marker_genes")
app_api.add_resource(celltypeAbundance, "/data/celltype_abundance")
# FIXME: this should not be a separate API endpoint
app_api.add_resource(plotsForSeachGenes, "/data/by_celltype_2_genes")


# Main loop
if __name__ == "__main__":
    # app.run(debug=True)
    app.run(host="0.0.0.0", port=5000)
