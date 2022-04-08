#!/bin/bash

# Adapted from python-igraph, original author Tamas Nepus:
# https://github.com/igraph/python-igraph/blob/709e7023aef4f4c4c93d385f4ed11adab6f7cbae/test.sh

###############################################################################

set -e

CLEAN=0
VENV_DIR=.venv
VERBOSE=0
MAIN_FLASK_FILE=main.py

if [ x$CLEAN = x1 ]; then
    rm -rf ${VENV_DIR}
fi

PYTHON=${VENV_DIR}/bin/python
PIP=${VENV_DIR}/bin/pip
PYTEST=${VENV_DIR}/bin/pytest

if [ ! -d ${VENV_DIR} ]; then
    python -m venv ${VENV_DIR}
fi
${VENV_DIR}/bin/pip install -r requirements.txt


if [ x$VERBOSE = x1 ]; then
  echo "${VENV_DIR}/bin/${PYTHON} ${MAIN_FLASK_FILE}"
fi
${PYTHON} ${MAIN_FLASK_FILE}
