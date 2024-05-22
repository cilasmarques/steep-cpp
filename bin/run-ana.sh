#!/bin/bash

current_dir=$(dirname -- "$0")
parent_dir=$(dirname -- "$current_dir")
cd -P -- "$parent_dir"

OUTPUT_DATA_PATH=./output

for i in $(seq -f "%02g" 1 3); do
  ./src/main "$@" &

  PID=$!

  # Inicia os scripts de monitoramento em background
  sh ./scripts/collect-cpu-usage.sh $PID > $OUTPUT_DATA_PATH/cpu.csv &
  sh ./scripts/collect-memory-usage.sh $PID > $OUTPUT_DATA_PATH/mem.csv &
  sh ./scripts/collect-disk-usage.sh $PID > $OUTPUT_DATA_PATH/disk.csv &
  sh ./scripts/collect-gpu-usage.sh $PID > $OUTPUT_DATA_PATH/gpu.csv &
  sh ./scripts/collect-gpu-memory-usage.sh $PID > $OUTPUT_DATA_PATH/mem-gpu.csv &

  # Aguarda o ./src/main terminar
  wait $PID

  # Para os scripts de monitoramento
  # kill $(pidof -s collect-cpu-usage.sh)
  # kill $(pidof -s collect-memory-usage.sh)
  # kill $(pidof -s collect-disk-usage.sh)
  # kill $(pidof -s collect-gpu-usage.sh)
  # kill $(pidof -s collect-gpu-memory-usage.sh)

  # Put the data in the final folder
  THREADS_NUM=`echo "$@" | grep -oP '(?<=-threads=)[0-9]+'`
  ANALYSIS_OUTPUT_PATH=$OUTPUT_DATA_PATH/analysis-$THREADS_NUM
  
  mkdir -p $ANALYSIS_OUTPUT_PATH/experiment${i}
  mv $OUTPUT_DATA_PATH/*.csv $ANALYSIS_OUTPUT_PATH/experiment${i}
  mv $OUTPUT_DATA_PATH/*.txt $ANALYSIS_OUTPUT_PATH/experiment${i}

  sleep 1
done

exit 0
