from flask import Flask, send_from_directory
from flask_restful import Resource, Api
from flask_cors import CORS
import pandas as pd
import json

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

class GeneExpression(Resource):
    def get(self):
        df = pd.read_excel("test_data_1.xlsx")
        df = df.round(decimals=3)
        print(df)
        data_string = df.to_json()  # string type
        return json.loads(data_string)  # json object

# New data set for gene expression across cell type
class GeneExpression_2(Resource):
    def get(self):
        df = pd.read_csv("pbmc3k_subset_data.csv")
        # this data is used for heatmap and scatter plot
        # for each gene, we want the average gene expression level in each cell type
        df_avg_exp = df.groupby(['Population']).mean()
        data_string = df_avg_exp.to_json()
        return json.loads(data_string)

class GeneExpression_3(Resource):
    def get(self):
        # if we want a vilon plot for all gene expression in B cells
        df = pd.read_csv("pbmc3k_subset_data.csv")
        df_bCells = df[df['Population'] == 'B cells'].drop(['CellBarcode'],axis=1)
        df_bCells.set_index('Population',inplace=True)
        data_dic = {}
        for gene in df_bCells.columns:
            exp_data = list(df_bCells[gene])
            data_dic[gene] = exp_data
        return data_dic

# this is an API endpoint (return data)
api.add_resource(GeneExpression, '/data')
api.add_resource(GeneExpression_2, '/data2')
api.add_resource(GeneExpression_3, '/data3')

if __name__ == '__main__':
    app.run(debug=True)