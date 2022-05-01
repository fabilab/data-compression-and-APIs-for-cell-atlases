# vim: fdm=indent
'''
author:     Fabio Zanini
date:       23/04/22
content:    Interpret text from Google TTS into a redirect URL.
'''
import numpy as np
import pandas as pd

from models import (
    gene_order,
    celltypes as celltypes_all,
    get_marker_genes,
    )
from validation import (
    validate_correct_genestr,
    validate_correct_celltypestr,
    )



# TODO
def interpret_text(text_raw):
    return text_raw
