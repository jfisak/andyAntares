#!/bin/bash

BASIC_INPUT=input.dat

# input model file
INPUT_MODEL_FILE=`cat $BASIC_INPUT | grep inputmodelFile | cut -d'=' -f 2`

INPUT_COMPOSITION_FILE=`cat $BASIC_INPUT | grep compose_adata_t.dat | cut -d'=' -f 2`


