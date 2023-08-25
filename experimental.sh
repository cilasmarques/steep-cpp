#!/bin/bash

FINAL_PATH=$1
METHOD=$2
INPUT_PATH=$3
OUTPUT_PATH=$4

for i in $(seq -f "%02g" 1 30)
do
	echo Executing ${i}th experiment
	mkdir $FINAL_PATH/experimento${i}
	time ./execute.sh LC08 $INPUT_PATH $OUTPUT_PATH $METHOD &&
	
	sleep 120
	mv $OUTPUT_PATH/*.csv $FINAL_PATH/experimento${i}
        mv $OUTPUT_PATH/*.txt $FINAL_PATH/experimento${i}
	rm -rf $OUTPUT_PATH/*
	echo Completing ${i}th experiment
done
