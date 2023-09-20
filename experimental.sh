#!/bin/bash

FINAL_PATH=$1
METHOD=$2
INPUT_PATH=$3
OUTPUT_PATH=$4
SENSOR=$5
LANDCOVER_DATA_FILE=$6

rm -rf $FINAL_PATH/*

for i in $(seq -f "%02g" 1 30)
do
	echo Executing ${i}th experiment
	mkdir $FINAL_PATH/experimento${i}
	time sh ./execute.sh $SENSOR $INPUT_PATH $OUTPUT_PATH $METHOD $LANDCOVER_DATA_FILE
	
	sleep 105
	mv $OUTPUT_PATH/*.csv $FINAL_PATH/experimento${i}
	rm -rf $OUTPUT_PATH/*
	echo Completing ${i}th experiment
done
