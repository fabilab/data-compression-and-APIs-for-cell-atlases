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
from helper import data_preprocessing, dataset_by_timepoint

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

@app.route('/heatmap_by_timepoints',methods=['GET'])
def page2():
    with open('page2.html') as f:
        response = f.read()
    return response

# new end point for timepoint dataset:
class geneExpTime(Resource):
    def get(self):
        genename = request.args.get('gene')
        plot_type = request.args.get('plot_type')
        data_type = request.args.get('data_type')
        data = None
        data = dataset_by_timepoint(genename)
        if data is None:
            return None
        
        if data_type == "log10":
            data = np.log10(0.1+data)
        
        if plot_type == 'hieracical':
            distance = pdist(data.values)
            Z = linkage(distance,optimal_ordering=True)
            new_order = leaves_list(Z)
            data = data.iloc[new_order]

        return data

class geneExp(Resource):
    def get(self):
        gene_names = request.args.get('gene_names')
        plot_type = request.args.get('plot_type')
        data_type = request.args.get('data_type')
        df = None
        df = data_preprocessing(gene_names)[0]
        if df is None:
            return None
        
        if data_type == "log10":
            df = np.log10(0.1+df)
        
        if plot_type == 'hieracical':
            distance = pdist(df.values)
            Z = linkage(distance,optimal_ordering=True)
            new_order = leaves_list(Z)
            df = df.iloc[new_order]
        
        return json.loads(df.to_json())

class plotsForSeachGenes(Resource):
    def get(self):

        gene_names = request.args.get('gene_names')
        df = data_preprocessing(gene_names)[0]
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
api.add_resource(geneExp, '/data')
api.add_resource(plotsForSeachGenes, '/2_genes')
api.add_resource(geneExpTime, '/data_timepoint')

if __name__ == '__main__':
    app.run(debug=True)