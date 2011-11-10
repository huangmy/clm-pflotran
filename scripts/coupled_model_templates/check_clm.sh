#!/bin/sh

MODEL=$1

if [ $MODEL = "clm" ]
then
    echo "This is CLM model"
    cp -f /ccs/home/gb9/COUPLED_MODEL_TEMPLATES/Srcfiles.CLM_PFLOTRAN Srcfiles
    cp -f /ccs/home/gb9/COUPLED_MODEL_TEMPLATES/Depends.CLM_PFLOTRAN Depends
fi

