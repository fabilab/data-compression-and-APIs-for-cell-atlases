from base64 import decode
import json
from flask import (
    Flask,
    send_from_directory,
    request,
    redirect,
    url_for,
    render_template,
)
from flask_restful import Api
from flask_cors import CORS

from api import (
    geneExp,
    geneExpTime,
    plotsForSeachGenes,
    geneFriends,
    geneExpTimeUnified,
    geneExpHyperoxia,
    geneExpDifferential,
    checkGenenames,
    markerGenes,
    celltypeAbundance,
)
from models import get_celltype_abundances
from validation.genes import validate_correct_genestr
from validation.timepoints import validate_correct_timepoint
from voice_recognition import mod as voice_control_blueprint
from text_recognition import mod as text_control_blueprint


##############################
app = Flask(__name__, static_url_path="/static", template_folder="templates")
app_api = Api(app)
# Note: this might be unsafe
CORS(app)
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


# Default redirects
@app.route("/heatmap_by_celltype", methods=["GET"])
def heatmap_by_celltype():
    """Landing page"""
    return redirect(url_for(
        'heatmap_by_celltype_genes',
        genestring=','.join([
            'Col1a1,Col2a1',
            'Adh1,Col13a1,Col14a1',
            'Tgfbi,Pdgfra,Crh,Hhip,Pdgfrb',
            'Pecam1,Gja5,Vwf,Car8,Car4',
            'Ptprc,Cd19,Gzma,Cd3d,Cd68',
            'Epcam',
            ])),
    )


@app.route("/heatmap_by_timepoint", methods=["GET"])
def heatmap_by_timepoint():
    """Heatmaps by timepoint (one per dataset)"""
    return redirect(url_for(
        'heatmap_by_timepoint_genes',
        genestring="Car4"),
    )


@app.route("/heatmap_unified", methods=["GET"])
def heatmap_unified():
    """One large heatmap across development, one gene"""
    return redirect(url_for(
        'heatmap_unified_genes',
        genestring='Car4'),
    )


@app.route("/heatmap_differential", methods=["GET"])
def heatmap_differential():
    """A sort of heatmap with differential expression"""
    return redirect(url_for(
        'heatmap_differential_genes',
        genestring=','.join([
            'Col1a1,Col2a1',
            'Adh1,Col13a1,Col14a1',
            'Tgfbi,Pdgfra,Crh,Hhip,Pdgfrb',
            'Pecam1,Gja5,Vwf,Car8,Car4',
            'Ptprc,Cd19,Gzma,Cd3d,Cd68',
            'Epcam',
            ])),
    )


# Generic views
@app.route("/heatmap_by_celltype/<genestring>", methods=['GET'])
def heatmap_by_celltype_genes(genestring):
    searchstring = genestring.replace(" ", "")
    return render_template(
            "heatmap_by_celltype.html",
            searchstring=searchstring,
            )


@app.route("/heatmap_by_timepoint/<genestring>", methods=['GET'])
def heatmap_by_timepoint_genes(genestring):
    searchstring = genestring.replace(" ", "")
    return render_template(
            "heatmap_by_timepoint.html",
            searchstring=searchstring,
            )


@app.route("/heatmap_unified/<genestring>", methods=['GET'])
def heatmap_unified_genes(genestring):
    searchstring = genestring.replace(" ", "")
    return render_template(
            "heatmap_unified.html",
            searchstring=searchstring,
            )


@app.route("/heatmap_hyperoxia/<genestring>", methods=["GET"])
def heatmap_hyperoxia(genestring):
    """A sort of heatmap with hyperoxia"""
    searchstring = genestring.replace(" ", "")
    return render_template(
            "heatmap_hyperoxia.html",
            searchstring=searchstring,
            )


# FIXME: this one should start to diverge
@app.route("/heatmap_differential/<genestring>", methods=["GET"])
def heatmap_differential_genes(genestring):
    """A sort of heatmap for differential gene expression"""
    searchstring = genestring.replace(" ", "")
    return render_template(
            "heatmap_differential.html",
            searchstring=searchstring,
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
app_api.add_resource(geneExpTime, "/data_timepoint")
app_api.add_resource(geneFriends, "/gene_friends")
app_api.add_resource(geneExpTimeUnified, "/data_heatmap_unified")
app_api.add_resource(geneExpHyperoxia, "/data/hyperoxia")
app_api.add_resource(geneExpDifferential, "/data/deg")
app_api.add_resource(checkGenenames, "/check_genenames")
app_api.add_resource(markerGenes, "/data/marker_genes")
app_api.add_resource(celltypeAbundance, "/data/celltype_abundance")
# FIXME: this should not be a separate API endpoint
app_api.add_resource(plotsForSeachGenes, "/data/by_celltype_2_genes")


# Main loop
if __name__ == "__main__":
    # app.run(debug=True)
    app.run(host="0.0.0.0", port=5000)
