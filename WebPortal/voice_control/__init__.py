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
from .interpret_text import text_to_url


mod = Blueprint('voice_control_blueprint', __name__)


@mod.route('/submit_audio', methods=['POST'])
def voice_control():
    # Extract audio Blob
    audio_data = request.files['audio_data']

    # Call Google API Speech To Text
    text = convert_audio_blob_to_text(audio_data)

    print('Google replied:', text)

    # Redirect to the correct endpoint
    url = text_to_url(text)
    return url
