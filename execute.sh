#!/bin/bash

cd -P -- "$(dirname -- "$0")"
echo "Script executed from: ${PWD}"

SENSOR=$1
INPUT_DATA_PATH=$2
OUTPUT_DATA_PATH=$3
METHOD=$4
LANDCOVER_DATA_FILE=$5
THREADS=$6

# if [ $SENSOR == LC08 ]
#   ./src/main $INPUT_DATA_PATH/B2.tif $INPUT_DATA_PATH/B3.tif $INPUT_DATA_PATH/B4.tif $INPUT_DATA_PATH/B5.tif $INPUT_DATA_PATH/B6.tif $INPUT_DATA_PATH/B10.tif $INPUT_DATA_PATH/B7.tif $INPUT_DATA_PATH/final_tal.tif $INPUT_DATA_PATH/MTL.txt $INPUT_DATA_PATH/station.csv $LANDCOVER_DATA_FILE $OUTPUT_DATA_PATH -meth=$METHOD -nan=-3.39999995214436425e+38 > $OUTPUT_DATA_PATH/timestamp.csv &
# then
# else
# 	./src/main $INPUT_DATA_PATH/B1.tif $INPUT_DATA_PATH/B2.tif $INPUT_DATA_PATH/B3.tif $INPUT_DATA_PATH/B4.tif $INPUT_DATA_PATH/B5.tif $INPUT_DATA_PATH/B6.tif $INPUT_DATA_PATH/B7.tif $INPUT_DATA_PATH/final_tal.tif $INPUT_DATA_PATH/MTL.txt $INPUT_DATA_PATH/station.csv $LANDCOVER_DATA_FILE $OUTPUT_DATA_PATH -meth=$METHOD -nan=-3.39999995214436425e+38 > $OUTPUT_DATA_PATH/timestamp.csv &
# fi

${PWD}/src/main $INPUT_DATA_PATH/B2.tif $INPUT_DATA_PATH/B3.tif $INPUT_DATA_PATH/B4.tif $INPUT_DATA_PATH/B5.tif $INPUT_DATA_PATH/B6.tif $INPUT_DATA_PATH/B10.tif $INPUT_DATA_PATH/B7.tif $INPUT_DATA_PATH/final_tal.tif $INPUT_DATA_PATH/MTL.txt $INPUT_DATA_PATH/station.csv $LANDCOVER_DATA_FILE $OUTPUT_DATA_PATH -meth=$METHOD -nan=-3.39999995214436425e+38 -threads=$THREADS > $OUTPUT_DATA_PATH/timestamp.csv &

PID=$(pidof main)
echo "PID: $PID"
ps -aux | grep $PID
sh scripts/collect-cpu-usage.sh $PID > $OUTPUT_DATA_PATH/cpu.csv &
sh scripts/collect-memory-usage.sh $PID > $OUTPUT_DATA_PATH/mem.csv &
sh scripts/collect-disk-usage.sh $PID > $OUTPUT_DATA_PATH/disk.csv &
