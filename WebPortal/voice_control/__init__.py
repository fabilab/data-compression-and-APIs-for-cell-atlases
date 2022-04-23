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


mod = Blueprint('voice_control_blueprint', __name__)


@mod.route('/voicesocket', methods=['POST'])
def voice_control():
    audio_data = request.files['audio_data']
    
    # FIXME: debugging
    with open('/tmp/audio.wav', 'wb') as audio:
        audio_data.save(audio)
    print('file uploaded successfully')

    return ''

    ## Talk to Google via supposed Python API of STT
    #speechClient = SpeechClient.create()
    #text_stream = speechClient.magicFromGoogle(audio_data)

    ## Needed?
    #speechClient.close()

    #actual_request = interpret_text_request(text_stream)
    #if actual_request is None:
    #    return {
    #        'data': '',
    #        'error': 'Request not understood',
    #    }

    ## Redirect to the correct endpoint
    #return redirect(url_for(actual_request))
