## Steps:

#### Install docker
```
curl -fSsL https://get.docker.com | bash
```

#### Install gdal
```
sudo add-apt-repository ppa:ubuntugis/ppa && sudo apt-get update 
sudo apt-get install gdal-bin
```

#### Install tiffs lib
```
sudo apt-get install libtiff5-dev
```

#### Make commands
* docker-landsat-download: 
    * Downloads the specified Landsat image.
* docker-landsat-preprocess: 
    * Preprocesses the Landsat image to ensure it's in the correct format for script input.
* fix-permissions: 
    * Adjusts the access permissions of Landsat images to ensure they are readable by the program.
* build-cpp:
    * Build the cpp code
* build-nvcc:
    * Build the cuda code
* exec-landsat*: 
    * Executes the script while generating an Nsight Systems (nsys) report.
* ncu-landsat*: 
    * Executes the script with an NVIDIA Compute Profiler (ncu) report.
* analysis-landsat*: 
    * Performs multiple executions of the script, each time generating an Nsight Systems (nsys) report.

