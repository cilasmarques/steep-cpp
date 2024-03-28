#!/bin/bash

# Executa ./src/main e passa todos os argumentos para ele
./src/main "$@" &

PID=$!
OUTPUT_DATA_PATH=./outputs

# Inicia os scripts de monitoramento em background
sh scripts/collect-cpu-usage.sh $PID > $OUTPUT_DATA_PATH/cpu.csv &
sh scripts/collect-memory-usage.sh $PID > $OUTPUT_DATA_PATH/mem.csv &
sh scripts/collect-disk-usage.sh $PID > $OUTPUT_DATA_PATH/disk.csv &
sh scripts/collect-gpu-usage.sh $PID > $OUTPUT_DATA_PATH/gpu.csv &
sh scripts/collect-gpu-memory-usage.sh $PID > $OUTPUT_DATA_PATH/mem-gpu.csv &

# Aguarda o ./src/main terminar
wait $PID

# Para os scripts de monitoramento
# kill $(pidof -s collect-cpu-usage.sh)
# kill $(pidof -s collect-memory-usage.sh)
# kill $(pidof -s collect-disk-usage.sh)
# kill $(pidof -s collect-gpu-usage.sh)
# kill $(pidof -s collect-gpu-memory-usage.sh)

exit 0