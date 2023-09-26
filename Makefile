GCC=g++
NVCC=nvcc
CXXFLAGS=-std=c++14 -ltiff

METHOD="$(METHOD)"
INPUT_DATA_PATH="$(INPUT_DATA_PATH)"

THREADS=256
LANDCOVER_DATA_FILE=""
OUTPUT_DATA_PATH="./output"


all:
	$(GCC) -g ./src/*.cpp -o ./src/main $(CXXFLAGS)

nvcc:
	$(NVCC) -g ./src/*.cu -o ./src/main $(CXXFLAGS)

clean:
	rm -rf ./src/main

clean-output:
	rm -rf ./output/*	

LC08:
	./src/main $(INPUT_DATA_PATH)/B2.tif $(INPUT_DATA_PATH)/B3.tif $(INPUT_DATA_PATH)/B4.tif $(INPUT_DATA_PATH)/B5.tif $(INPUT_DATA_PATH)/B6.tif $(INPUT_DATA_PATH)/B10.tif $(INPUT_DATA_PATH)/B7.tif $(INPUT_DATA_PATH)/final_tal.tif $(INPUT_DATA_PATH)/MTL.txt $(INPUT_DATA_PATH)/station.csv $(LANDCOVER_DATA_FILE) $(OUTPUT_DATA_PATH) -meth=$(METHOD) -nan=-3.39999995214436425e+38 -threads=$(THREADS) > $(OUTPUT_DATA_PATH)/timestamp.csv

TM05:
	./src/main $(INPUT_DATA_PATH)/B1.tif $(INPUT_DATA_PATH)/B2.tif $(INPUT_DATA_PATH)/B3.tif $(INPUT_DATA_PATH)/B4.tif $(INPUT_DATA_PATH)/B5.tif $(INPUT_DATA_PATH)/B6.tif $(INPUT_DATA_PATH)/B7.tif $(INPUT_DATA_PATH)/final_tal.tif $(INPUT_DATA_PATH)/MTL.txt $(INPUT_DATA_PATH)/station.csv $(LANDCOVER_DATA_FILE) $(OUTPUT_DATA_PATH) -meth=$(METHOD) -nan=-3.39999995214436425e+38 -threads=$(THREADS) > $(OUTPUT_DATA_PATH)/timestamp.csv