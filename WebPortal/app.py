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

from ca_API import (
    geneExp,
    geneExpTime,
    plotsForSeachGenes,
    geneFriends,
    geneExpTimeUnified,
    checkGenenames,
)
from voice_control import mod as voice_control_blueprint


app = Flask(__name__, static_url_path="/static", template_folder="templates")
api = Api(app)
# Note: this might be unsafe
CORS(app)


# Views
@app.route("/", methods=["GET"])
def index():
    """Landing page"""
    with open("index.html") as f:
        response = f.read()
    return response


@app.route("/celltype/<genestring>", methods=['GET'])
def heatmap_by_celltype_genes(genestring):
    searchstring = genestring.replace(" ", "")
    return render_template(
            "heatmap_celltype.html",
            searchstring=searchstring,
            )


@app.route("/heatmap_by_timepoints", methods=["GET"])
def heatmap_by_timepoint():
    """Heatmaps by timepoint (one per dataset)"""
    with open("heatmap_by_timepoint.html") as f:
        response = f.read()
    return response


@app.route("/heatmap_unified", methods=["GET"])
def heatmap_unified():
    """One large heatmap across development, one gene"""
    with open("heatmap_unified.html") as f:
        response = f.read()
    return response


@app.route("/heatmap_unified/<genestring>", methods=['GET'])
def heatmap_unified_genes(genestring):
    searchstring = genestring.replace(" ", "")
    return render_template(
            "heatmap_unified.html",
            searchstring=searchstring,
            )


@app.route("/voice_control", methods=["GET"])
def voice_control():
    """Name says it all"""
    with open("voice_control.html") as f:
        response = f.read()
    return response


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


# API endpoints
api.add_resource(geneExp, "/data")
api.add_resource(plotsForSeachGenes, "/2_genes")
api.add_resource(geneExpTime, "/data_timepoint")
api.add_resource(geneFriends, "/gene_friends")
api.add_resource(geneExpTimeUnified, "/data_heatmap_unified")
api.add_resource(checkGenenames, "/check_genenames")


# Blueprints
app.register_blueprint(voice_control_blueprint)


# Main loop
if __name__ == "__main__":
    # app.run(debug=True)
    app.run(host="0.0.0.0", port=5000)
