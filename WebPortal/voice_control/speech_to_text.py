# vim: fdm=indent
'''
author:     Fabio Zanini
date:       23/04/22
content:    Convert speech to text using cloud API
'''
import os
from google.cloud import speech_v1


def check_ensure_gcloud_credentials():
    '''Ensure credentials for google cloud are found'''
    from pathlib import Path

    if os.getenv('GOOGLE_APPLICATION_CREDENTIALS') is None:
        key_path = Path(__file__).parent.parent / 'gcloud_key' / 'zaninilab_key.json'
        os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = str(key_path)


def convert_audio_blob_to_text(audio_blob):

    # convert the input audio into whatever the client requires
    config = speech_v1.RecognitionConfig(
        language_code="en-US",
        encoding='FLAC',
    )

    # FIXME: we have to convert the blob: <FileStorage: 'blob' ('audio/webm')>
    # into bytes for google cloud
    audio_bytes = audio_blob.read()
    audio = speech_v1.RecognitionAudio(
            content=audio_bytes,
    )

    request = speech_v1.RecognizeRequest(
        config=config,
        audio=audio,
    )

    print('got here, ensuring credentials')

    check_ensure_gcloud_credentials()
    client = speech_v1.SpeechClient()
    response = client.recognize(request=request)

    # TODO: extract text and confidence
    results = response.results
    if len(results) == 0:
        return ''

    for result in results:
        candidates = result.alternatives
        if len(candidates) == 0:
            continue
        cand = candidates[0]
        text = cand.transcript
        if hasattr(cand, 'confidence') and cand.confidence < 0.6:
            continue

        # Good text (sequence of words)
        return text

    return ''
