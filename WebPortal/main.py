from base64 import decode
from flask import Flask, send_from_directory, request
from flask_restful import Api
from flask_cors import CORS

from ca_API import (
    geneExp,
    geneExpTime,
    plotsForSeachGenes,
    geneFriends,
    geneExpTimeUnified,
)


app = Flask(__name__, static_url_path='/static')
api = Api(app)
# Note: this might be unsafe
CORS(app)


# Views
@app.route('/', methods=['GET'])
def index():
    '''Landing page'''
    with open('index.html') as f:
        response = f.read()
    return response


@app.route('/heatmap_by_timepoints', methods=['GET'])
def heatmap_by_timepoint():
    '''Heatmaps by timepoint (one per dataset)'''
    with open('heatmap_by_timepoint.html') as f:
        response = f.read()
    return response


@app.route('/heatmap_unified', methods=['GET'])
def heatmap_unified():
    '''One large heatmap across development, one gene'''
    with open('heatmap_unified.html') as f:
        response = f.read()
    return response


# Static assets (JS/CSS)
@app.route('/js/<path:path>')
def send_js(path):
    '''JavaScript assets'''
    return send_from_directory('static/js', path)


@app.route('/css/<path:path>')
def send_css(path):
    '''CSS stylesheets'''
    return send_from_directory('static/css', path)


@app.route('/favicon.ico')
def favicon():
    return send_from_directory('static', 'favicon.ico')


# API endpoints
api.add_resource(geneExp, '/data')
api.add_resource(plotsForSeachGenes, '/2_genes')
api.add_resource(geneExpTime, '/data_timepoint')
api.add_resource(geneFriends, '/gene_friends')
api.add_resource(geneExpTimeUnified, '/data_heatmap_unified')


# Main loop
if __name__ == '__main__':
    #app.run(debug=True)
    app.run(host='0.0.0.0', port=5000)
