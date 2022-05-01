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
from .speech_to_text import convert_audio_blob_to_text
from text_recognition import text_to_response


mod = Blueprint('voice_control_blueprint', __name__)


@mod.route('/submit_audio', methods=['POST'])
def voice_control():
    # Extract audio Blob
    audio_data = request.files['audio_data']

    # Call Google API Speech To Text
    text = convert_audio_blob_to_text(audio_data)
    print('Google replied:', text)

    # Use internal algos to convert a text command to a JSON response
    response = text_to_response(text)

    return response
