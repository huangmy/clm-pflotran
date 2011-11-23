#!/bin/sh

MODEL=$1

if [ $MODEL = "clm" ]
then
    echo "This is CLM model"
    cp -f $CLM_COUPLED_MODEL/scripts/coupled_model_templates/Srcfiles.CLM_PFLOTRAN Srcfiles
    cp -f $CLM_COUPLED_MODEL/scripts/coupled_model_templates/Depends.CLM_PFLOTRAN Depends
fi

