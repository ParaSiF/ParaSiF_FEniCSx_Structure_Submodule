#!/bin/sh

# Run from this directory
cd ${0%/*} || exit 1

export PYTHONPATH=${PWD}/../../src:$PYTHONPATH
export PYTHONPATH=${PWD}../../../../../coupling_lib/MUI/wrappers/Python:$PYTHONPATH

domainFluid=${PWD}/dummyOF
domainStructure=${PWD}/structureDomain

# Ranks set to each domain
numProcsFluid=1
numProcsStructure=1

solverFluid=./PUSHER_FETCHER_1
solverStructure=structureDomainRun.py

cd ${domainFluid}

# Create build folder
mkdir build && cd build

# Check if an argument was provided
cmake -DCMAKE_PREFIX_PATH=$(pwd)/../../../../../../../coupling_lib/MUI ..

# Run make to build the executable
make 2>&1 | tee make.log && cd ..
cp build/PUSHER_FETCHER_1 ./

cd ..

# parallel run
mpirun -np ${numProcsFluid} -wdir ${domainFluid} ${solverFluid} -parallel -coupled :\
       -np ${numProcsStructure} -wdir ${domainStructure} python3 -m mpi4py ${solverStructure} 2>&1 | tee output.log

echo "Done"

# ----------------------------------------------------------------- end-of-file