# Code for the Lap2Go Project at the ESGI140 in Portugal

## Requirements

+ cmake
+ python3
+ a reasonably up-to-date C++ compiler
+ the boost C++ libraries

## Installation instructions:

    git clone git@github.com:florianwechsung/ESGI140.git
    cd ESGI140
    mkdir build
    cd build
    cmake -DPYBIND11_PYTHON_VERSION=3.5 ..
    make
    export PYTHONPATH=`pwd`:$PYTHONPATH

Then run the sample code by doing

    cd ../python/
    python3 driver.py

This will produce a video with the race simulation:

![Sample Simulation](sample.png)
