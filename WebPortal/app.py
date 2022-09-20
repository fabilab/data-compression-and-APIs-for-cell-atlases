from flask import Flask, send_from_directory,request,render_template, send_from_directory
from flask_restful import Resource, Api
from flask_cors import CORS
import os
from scipy.cluster.hierarchy import linkage,leaves_list
from scipy.spatial.distance import pdist
from ca_data_access import read_file_average_exp,read_file_proportion_exp, result_datasets, dataset_unified, result_unified_by_cell, marker_genes_expression


app = Flask(__name__, static_url_path='/static', template_folder='templates')
api = Api(app)
# Note: this might be unsafe
CORS(app)

@app.route('/favicon.ico')
def favicon():
    return send_from_directory(os.path.join(app.root_path, 'static/img'),
                          'favicon.ico',mimetype='image/x-icon')

@app.route('/logo.png')
def logo():
    return send_from_directory(os.path.join(app.root_path, 'static/img'),
                          'logo.png',mimetype='image/png')

@app.route('/about_img1.png')
def about():
    return send_from_directory(os.path.join(app.root_path, 'static/img'),
                          'about_img1.png',mimetype='image/png')

# this is an endpoint (gives access to the file)
@app.route('/js/<path:path>')
def send_js(path):
    return send_from_directory('static/js', path)

@app.route('/css/<path:path>')
def send_css(path):
    return send_from_directory('static/css', path)

@app.route('/',methods=['GET'])
def home():
    return render_template('about.html')

@app.route('/dataExplore',methods=['GET'])
def dataExplore():
    return render_template('dataExplore.html')

@app.route('/resources',methods=['GET'])
def resources():
    return render_template('resources.html')

@app.route('/userGuide',methods=['GET'])
def userGuide():
    return render_template('userGuide.html')

@app.route('/package',methods=['GET'])
def package():
    return render_template('/showPackage/index.html')

@app.route('/packageModule',methods=['GET'])
def packageModule():
    return render_template('/showPackage/modules.html')

@app.route('/heatmap_by_celltypes',methods=['GET'])
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

class getAllCellTypes(Resource):
    def get(self):
        result = {
            'Fibroblast': [
                'Adventitial fibroblast',
                'Early adventitial fibroblast',
                'Fibroblast precursor',
                'Proliferating fibroblast',
                'Alveolar fibroblast',
                'Early alveolar fibroblast'
            ],
            'ASM and MyoF':[
                'Myofibroblast and smooth muscle precursor',
                'Airway smooth muscle',
                'Early airway smooth muscle',
                'Myofibroblast',
                'Proliferating myofibroblast'
            ],
            'Mural': [
                'Vascular smooth muscle',
                'Pericyte',
                'Proliferating pericyte'
            ],
            'Immune system':[
                'B cell',
                'DC I',
                'DC II',
                'DC III',
                'IL cell',
                'Mac I',
                'Mac II',
                'Mac III',
                'Mac IV',
                'Mac V',
                'NK cell',
                'T cell',
                'basophil',
                'mast cell',
                'neutrophil'
            ],
            'Lungs':[
                'Alveolar type I',
                'Alveolar type II',
            ],
            'Vasculature':[
                'Arterial EC I',
                'Arterial EC II',
                'Car4+ capillaries',
                'Early Car4- capillaries',
                'Late Car4- capillaries',
                'Lymphatic EC',
                'Nonproliferative embryonic EC',
                'Proliferative EC',
                'Venous EC'
            ],
            'Others':[
                'Striated muscle'
            ]
        }
        return result

class getAllGeneNames(Resource):
    def get(self):
        _, df = read_file_average_exp('celltype')
        return list(df.index)


class dataDatasets (Resource):
    def get(self):
        gene = request.args.get('gene')
        result = result_datasets('celltype_dataset_timepoint',gene)
        return result

class dataGeneral(Resource):
    def get(self):
        gene_names = request.args.get('gene_names')
        df_average = None
        df_proportion = None
        # df = data_preprocessing(gene_names,'celltype')
        res, data = read_file_average_exp("celltype",gene_names)
        if not res:
            return data, 400
        df_average = data.T
        res, data = read_file_proportion_exp("celltype",gene_names)
        if not res or df_average is None or data is None:
            return None, 400
        df_proportion = data.T

        # new cell types order
        distance = pdist(df_average.values)
        Z = linkage(distance,optimal_ordering=True)
        new_order = leaves_list(Z)
        
        # new genes order
        distance_g = pdist(data.values)
        Z_g = linkage(distance_g,optimal_ordering=True)
        new_order_g = leaves_list(Z_g)

        response = {
            'result_average': df_average.to_dict(),
            'result_proportion': df_proportion.to_dict(),
            'hierarchicalGeneOrder': data.index[new_order_g].tolist(),  # new order of genes
            'hierarchicalCelltypeOrder': df_average.index[new_order].tolist(),  # new order of cell types
        }
        return response

class dataScatter(Resource):
    def get(self):
        gene_names = request.args.get('gene_names')
        # df = data_preprocessing(gene_names,'celltype')
        res, df = read_file_average_exp("celltype",gene_names)
        if not res or df is None:
            return df, 400
        a_gene_names = gene_names.split(",")
        print(a_gene_names)
        if len(a_gene_names) == 2:
            result = {}
            plot_df = df.filter(items = a_gene_names, axis=0)
            gene1 = plot_df.index[0]
            gene2 = plot_df.index[1]
            
            gene1_expr = list(plot_df.loc[gene1])
            gene2_expr = list(plot_df.loc[gene2])

            result["gene1_name"]= gene1
            result["gene2_name"] = gene2
            result["gene1_expr"] = gene1_expr
            result["gene2_expr"] = gene2_expr
            result["cell_types"] = list(plot_df.columns)

        return result

class dataUnified(Resource):
    def get(self):
        genename = request.args.get('gene')
        return dataset_unified(genename)

class UnifiedByCell(Resource):
    def get(self):
        celltype = request.args.get('celltype')
        genes = request.args.get('genes')
        res, data = result_unified_by_cell(celltype,genes)
        if not res:
            return data, 400
        return data

class dataMarkerGenes(Resource):
    def get(self):
        celltype = request.args.get('celltype')
        return marker_genes_expression(celltype)

# this is an API endpoint (return data)
api.add_resource(getAllCellTypes, '/all_cell_types')
api.add_resource(dataGeneral, '/data_general')
api.add_resource(getAllGeneNames, '/all_gene_names')
api.add_resource(dataScatter, '/data_scatter')
api.add_resource(dataDatasets, '/data_datasets')
api.add_resource(dataUnified, '/data_unified')
api.add_resource(UnifiedByCell, '/data_unified_by_cell')
api.add_resource(dataMarkerGenes, '/data_markers')

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
    # app.run(debug=True, host='127.0.0.1', port=5000)