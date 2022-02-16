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
class H5GeneExpression(Resource):
    def get(self):
        h5_data = h5py.File('./static/scData/condensed_lung_atlas.h5',"r")
        df = pd.DataFrame(data=np.array(h5_data['cell_type']\
            ['gene_expression_average']['block0_values']),\
            index=np.array(h5_data['cell_type']['gene_expression_average']['axis1'])\
            ,columns=np.array(h5_data['cell_type']['gene_expression_average']['axis0'])).T
        
        # current index in the dataframe is writtern as binary string.
        # We need to convert it into normal string
        new_index=[]
        for i in df.index:
            new_index.append(i.decode('utf-8'))
        # Similarly for columns name
        new_columns=[]
        for i in df.columns:
            new_columns.append(i.decode('utf-8'))

        df.index = new_index
        df.columns = new_columns

        # # Select 5 genes of interest:
        # plot_df = df.filter(items = ['Car4','Vwf', 'Col1a1', 'Ptprc', 'Ms4a1'], axis=0)
        # plot_data = plot_df.to_json()
        
        # get the name of genes input by the web user
        gene_names = request.args.get('gene_names')
        if gene_names is None:
            plot_data = df.T
        else:
            a_gene_names = gene_names.split(",")
            plot_df = df.filter(items = a_gene_names,axis=0)
            plot_data = plot_df.T
        plot_data = np.log10(0.1+plot_data)

        # Hierachical Clustering
        distance = pdist(plot_data.values)
        Z = linkage(distance,optimal_ordering=True)

        new_order = leaves_list(Z)
        plot_data = plot_data.iloc[new_order]

        return json.loads(plot_data.to_json())

class plotsForSeachGenes(Resource):
    def get(self):
        h5_data = h5py.File('./static/scData/condensed_lung_atlas.h5',"r")
        df = pd.DataFrame(data=np.array(h5_data['cell_type']\
            ['gene_expression_average']['block0_values']),\
            index=np.array(h5_data['cell_type']['gene_expression_average']['axis1'])\
            ,columns=np.array(h5_data['cell_type']['gene_expression_average']['axis0'])).T
        
        new_index=[]
        for i in df.index:
            new_index.append(i.decode('utf-8'))
        # Similarly for columns name
        new_columns=[]
        for i in df.columns:
            new_columns.append(i.decode('utf-8'))

        df.index = new_index
        df.columns = new_columns

        gene_names = request.args.get('gene_names')
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
api.add_resource(H5GeneExpression, '/data')
api.add_resource(plotsForSeachGenes, '/2_genes')

if __name__ == '__main__':
    app.run(debug=True)