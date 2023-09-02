#!/bin/bash

SENSOR=$1
INPUT_PATH=$2
OUTPUT_PATH=$3
METHOD=$4

make run INPUT_DATA_PATH=$INPUT_PATH OUTPUT_DATA_PATH=$OUTPUT_PATH METHOD=$METHOD

PID=$(pidof main)
ps -aux | grep $PID
sh scripts/collect-cpu-usage.sh $PID > $OUTPUT_PATH/cpu.csv &
sh scripts/collect-memory-usage.sh $PID > $OUTPUT_PATH/mem.csv &
