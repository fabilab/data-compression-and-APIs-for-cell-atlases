# import libraries
from flask import Flask, jsonify, request

# register the web app into a variable
app = Flask(__name__)

# function that return a 'hello world' text in JSON format 
@app.route('/greeting',methods=['GET'])
def helloworld():
	return 'Hello World'

if __name__ == '__main__':
	app.run(debug=True)