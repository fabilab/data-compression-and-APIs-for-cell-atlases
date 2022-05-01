# vim: fdm=indent
'''
author:     Fabio Zanini
date:       23/04/22
content:    Voice control backend module for compressed atlases
'''
from base64 import decode
import json
from flask import (
        Flask, send_from_directory, request, redirect, url_for,
        Blueprint,
        )
from .interpret_text import interpret_text
from .text_to_url import text_to_url


mod = Blueprint('text_control_blueprint', __name__)


@mod.route('/submit_text', methods=['GET'])
def text_control():
    # Extract audio Blob
    text_raw = request.args.get("text");

    # Call Google API Speech To Text
    text_int = interpret_text(text_raw)

    print('Text raw:', text_raw)
    print('Text interpreted:', text_int)

    # Redirect to the correct endpoint
    response = text_to_url(text_int)
    return response
