#!/bin/bash

SANDBOX=$(dirname $(readlink -fm $0))
cd $SANDBOX

## Checking args
if [ $# -ne 3 ]
then
  echo "Usage: $0 landsat_X PPPRRR YYYY-MM-DD"
  exit 1
fi

## args
IMAGE_DATASET=$1
IMAGE_PATHROW=$2
IMAGE_DATE=$3

echo $SANDBOX

cd $SANDBOX/download
bash ./run.sh $SANDBOX/output $IMAGE_DATASET $IMAGE_PATHROW $IMAGE_DATE

# bash $SANDBOX/preprocess/run.shZ $IMAGE_DATASET $IMAGE_PATHROW $IMAGE_DATE
# docker build ./preprocess/Dockerfile --build-arg IMAGE_DATASET=$IMAGE_DATASE --build-arg IMAGE_PATHROW= $IMAGE_PATHROW --build-arg IMAGE_DATE=$IMAGE_DATE -t landsat-preproc
