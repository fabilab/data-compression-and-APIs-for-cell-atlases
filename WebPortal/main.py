from base64 import decode
from flask import Flask, send_from_directory,request
from flask_restful import Resource, Api
from flask_cors import CORS
import pandas as pd
import h5py
import numpy as np
import json
from scipy.cluster.hierarchy import linkage,leaves_list
from scipy.spatial.distance import pdist
from helper import data_preprocessing

app = Flask(__name__, static_url_path='/static')
api = Api(app)
# Note: this might be unsafe
CORS(app)

# this is an endpoint (gives access to the file)
@app.route('/js/<path:path>')
def send_js(path):
    return send_from_directory('static/js', path)

@app.route('/css/<path:path>')
def send_css(path):
    return send_from_directory('static/css', path)

@app.route('/',methods=['GET'])
def helloworld():
    with open('index.html') as f:
        response = f.read()
    return response

# use Carsten's h5 data
class geneExpOriginal(Resource):
    def get(self):
        # # Select 5 genes of interest:
        # plot_df = df.filter(items = ['Car4','Vwf', 'Col1a1', 'Ptprc', 'Ms4a1'], axis=0)

        # get the name of genes input by the web user
        gene_names = request.args.get('gene_names')
        df = data_preprocessing(gene_names)
        if df is None:
            return None
        return json.loads(df.to_json())

class geneExpLog(Resource):
    def get(self):
        gene_names = request.args.get('gene_names')
        df = data_preprocessing(gene_names)
        if df is None:
            return None
        log_data = np.log10(0.1+df)

        return json.loads(log_data.to_json())

class geneExpHieracical(Resource):
    def get(self):
        gene_names = request.args.get('gene_names')
        df = data_preprocessing(gene_names)
        if df is None:
            return None
        
        hierachical_data = np.log10(0.1+df)
        # Hierachical Clustering
        distance = pdist(hierachical_data.values)
        Z = linkage(distance,optimal_ordering=True)

        new_order = leaves_list(Z)
        hierachical_data = hierachical_data.iloc[new_order]
        return json.loads(hierachical_data.to_json())

class geneExpHieracicalOriginal(Resource):
    def get(self):
        gene_names = request.args.get('gene_names')
        df = data_preprocessing(gene_names)
        if df is None:
            return None
        
        hierachical_data = df
        # Hierachical Clustering
        distance = pdist(hierachical_data.values)
        Z = linkage(distance,optimal_ordering=True)

        new_order = leaves_list(Z)
        hierachical_data = hierachical_data.iloc[new_order]
        return json.loads(hierachical_data.to_json())

class plotsForSeachGenes(Resource):
    def get(self):

        gene_names = request.args.get('gene_names')
        df = data_preprocessing(gene_names)
        if df is None:
            return None
        df = df.T
        a_gene_names = gene_names.split(",")
        if len(a_gene_names) == 2:
            result = {}
            plot_df = df.filter(items = a_gene_names, axis=0)
            gene1 = plot_df.index[0]
            gene2 = plot_df.index[1]

            gene1_expr = list(plot_df.loc[gene1])
            gene2_expr = list(plot_df.loc[gene2])
            # plot_data = plot_df.to_json()
            result["gene1_name"]= gene1
            result["gene2_name"] = gene2
            result["gene1_expr"] = gene1_expr
            result["gene2_expr"] = gene2_expr
            result["cell_types"] = list(plot_df.columns)

        return result

# this is an API endpoint (return data)
api.add_resource(geneExpOriginal, '/dataOrigin')
api.add_resource(geneExpLog, '/dataLog')
api.add_resource(geneExpHieracical, '/dataHierachical')
api.add_resource(geneExpHieracicalOriginal, '/dataHierachicalOriginal')
api.add_resource(plotsForSeachGenes, '/2_genes')

if __name__ == '__main__':
    app.run(debug=True)