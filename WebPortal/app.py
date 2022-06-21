from base64 import decode
from urllib import response
from flask import Flask, send_from_directory,request,render_template
from flask_restful import Resource, Api
from flask_cors import CORS
import pandas as pd
import h5py
import numpy as np
import json
from scipy.cluster.hierarchy import linkage,leaves_list
from scipy.spatial.distance import pdist
from ca_data_access import read_file,data_preprocessing, dataset_by_timepoint, dataset_unified, select_marker_genes
import time

app = Flask(__name__, static_url_path='/static',template_folder='templates')
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
# def helloworld():
def page1():
    return render_template('byCelltype.html',highlight='home_button',search_box='Car4,Vwf,Col1a1,Ptprc,Ms4a1')

@app.route('/heatmap_by_timepoints',methods=['GET'])
def page2():
    # with open('page2.html') as f:
    #     response = f.read()
    # return response
    return render_template('byDataset.html',highlight='page2_button',search_box='Car4')


@app.route('/unified',methods=['GET'])
def page3():
    return render_template('Unified.html',highlight='page3_button',search_box='Car4')

@app.route('/marker',methods=['GET'])
def page4():
    return render_template('markerGenes.html',highlight='page4_button',search_box='Car4')

class getAllCelltypes(Resource):
    def get(self):
        df = read_file('celltype','Car4')
        return list(df.columns)

class geneNames(Resource):
    def get(self):
        df = read_file('celltype')
        return list(df.index)

# new end point for timepoint dataset:
class geneExpTime(Resource):
    def get(self):
        genename = request.args.get('gene')
        # if plottype == 'hieracical':
        #     distance = pdist(gene_exp_df.values)
        #     # print(distance)
        #     Z = linkage(distance,optimal_ordering=True)
        #     new_order = leaves_list(Z)
        #     gene_exp_df = gene_exp_df.iloc[new_order]
        
        result = dataset_by_timepoint(genename,'celltype_dataset_timepoint')
        return result

class geneExp(Resource):
    def get(self):
        ######### 2
        start = time.time()
        print("geneEXP start")
        gene_names = request.args.get('gene_names')
        df = None
        df = data_preprocessing(gene_names,'celltype')
        if df is None:
            return None

        distance = pdist(df.values)
        Z = linkage(distance,optimal_ordering=True)
        new_order = leaves_list(Z)
        # df = df.iloc[new_order]

        response = {
            'result': df.to_dict(),
            'hierarchicalCelltypeOrder': df.index[new_order].tolist(),  # new order of the celltype
        }
        ######## 7 (store the result in a variable,print it)
        end = time.time()
        print("geneEXP end")
        print(end - start)
        return response

class plotsForSeachGenes(Resource):
    def get(self):

        gene_names = request.args.get('gene_names')
        df = data_preprocessing(gene_names,'celltype')
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

class geneExpUnified(Resource):
    def get(self):
        genename = request.args.get('gene')

        return dataset_unified(genename)

class markerGenes(Resource):
    def get(self):
        celltype = request.args.get('celltype')
        print('I am here ' + celltype)
        return select_marker_genes(celltype)

# this is an API endpoint (return data)
api.add_resource(getAllCelltypes, '/all_cell_types')
api.add_resource(geneExp, '/data')
api.add_resource(geneNames, '/all_gene_names')
api.add_resource(plotsForSeachGenes, '/2_genes')
api.add_resource(geneExpTime, '/data_timepoint')
api.add_resource(geneExpUnified, '/data_unified')
api.add_resource(markerGenes, '/markers_page')

if __name__ == '__main__':
    app.run(debug=True)