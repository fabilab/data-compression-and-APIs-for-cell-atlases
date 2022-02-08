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

# this is an API endpoint (return data)
api.add_resource(GeneExpression, '/data')

if __name__ == '__main__':
    app.run(debug=True)