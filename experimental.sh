#!/bin/bash

cd -P -- "$(dirname -- "$0")"

FINAL_PATH=$1
METHOD=$2
INPUT_PATH=$3
OUTPUT_PATH=$4
SENSOR=$5
LANDCOVER_DATA_FILE=$6
THREADS=$7

rm -rf $FINAL_PATH/*

for i in $(seq -f "%02g" 1 3); do
  echo "Executing ${i}th experiment"

  # Start the ./src/main process
  sh ./execute.sh $SENSOR $INPUT_PATH $OUTPUT_PATH $METHOD $LANDCOVER_DATA_FILE $THREADS

  # Capture the PID of the last background process
  PID=$(pidof -s main)
  echo "PID: $PID"

  # Wait for the background process to finish
  while kill -0 $PID 2>/dev/null; do
    sleep 1
  done

  # Put the data in the final folder
  mkdir -p $FINAL_PATH/experiment${i}
  mv $OUTPUT_PATH/*.csv $FINAL_PATH/experiment${i}
  mv $OUTPUT_PATH/*.txt $FINAL_PATH/experiment${i}
  rm -rf $OUTPUT_PATH/*

  echo "Completing ${i}th experiment"
done
