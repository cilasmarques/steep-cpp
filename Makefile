GCC=g++
NVCC=nvcc
CXXFLAGS=-std=c++14 -ltiff

## Docker download and preprocessing
DOCKER_OUTPUT_PATH="/home/saps/output"
DOCKER_LANDSAT="landsat_8"
DOCKER_PATHROW="215065"
DOCKER_DATE="2017-05-11"
IMAGES_OUTPUT=$(PWD)/scenes/output

## Execution
METHOD="$(METHOD)"
INPUT_DATA_PATH="$(INPUT_DATA_PATH)"
THREADS=256
BLOCKS=6504
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
	./src/main $(INPUT_DATA_PATH)/B2.tif $(INPUT_DATA_PATH)/B3.tif $(INPUT_DATA_PATH)/B4.tif $(INPUT_DATA_PATH)/B5.tif $(INPUT_DATA_PATH)/B6.tif $(INPUT_DATA_PATH)/B10.tif $(INPUT_DATA_PATH)/B7.tif $(INPUT_DATA_PATH)/final_tal.tif $(INPUT_DATA_PATH)/MTL.txt $(INPUT_DATA_PATH)/station.csv $(LANDCOVER_DATA_FILE) $(OUTPUT_DATA_PATH) -meth=$(METHOD) -nan=-3.39999995214436425e+38 -threads=$(THREADS) -blocks=$(BLOCKS) > $(OUTPUT_DATA_PATH)/timestamp.csv

TM05:
	./src/main $(INPUT_DATA_PATH)/B1.tif $(INPUT_DATA_PATH)/B2.tif $(INPUT_DATA_PATH)/B3.tif $(INPUT_DATA_PATH)/B4.tif $(INPUT_DATA_PATH)/B5.tif $(INPUT_DATA_PATH)/B6.tif $(INPUT_DATA_PATH)/B7.tif $(INPUT_DATA_PATH)/final_tal.tif $(INPUT_DATA_PATH)/MTL.txt $(INPUT_DATA_PATH)/station.csv $(LANDCOVER_DATA_FILE) $(OUTPUT_DATA_PATH) -meth=$(METHOD) -nan=-3.39999995214436425e+38 -threads=$(THREADS) -blocks=$(BLOCKS) > $(OUTPUT_DATA_PATH)/timestamp.csv

docker-run-download:
	docker run -v $(IMAGES_OUTPUT):$(DOCKER_OUTPUT_PATH) \
		-e OUTPUT_PATH=$(DOCKER_OUTPUT_PATH) \
		-e LANDSAT=$(DOCKER_LANDSAT) \
		-e PATHROW=$(DOCKER_PATHROW) \
		-e DATE=$(DOCKER_DATE) \
		cilasmarques/landsat-download

docker-run-preprocess:
	docker run -v $(IMAGES_OUTPUT):$(DOCKER_OUTPUT_PATH) \
		-e OUTPUT_PATH=$(DOCKER_OUTPUT_PATH) \
		-e LANDSAT=$(DOCKER_LANDSAT) \
		-e PATHROW=$(DOCKER_PATHROW) \
		-e DATE=$(DOCKER_DATE) \
		cilasmarques/landsat-preprocess:latest

