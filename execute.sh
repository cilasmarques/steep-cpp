#!/bin/bash

SENSOR=$1
INPUT_PATH=$2
OUTPUT_PATH=$3
METHOD=$4

./src/main ./input/l08/B2.tif ./input/l08/B3.tif ./input/l08/B4.tif ./input/l08/B5.tif ./input/l08/B6.tif ./input/l08/B10.tif ./input/l08/B7.tif ./input/l08/final_tal.tif ./input/l08/MTL.txt ./input/l08/station.csv ./input/l08/land_cover_final.tif $OUTPUT_PATH -meth=$METHOD -nan=-3.39999995214436425e+38 > $OUTPUT_PATH/timestamp.txt &

PID=$(pidof main)
ps -aux | grep $PID
sh scripts/collect-cpu-usage.sh $PID > $OUTPUT_PATH/cpu.csv &
sh scripts/collect-memory-usage.sh $PID > $OUTPUT_PATH/mem.csv &
