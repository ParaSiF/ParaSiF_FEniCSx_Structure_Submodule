#!/bin/sh

domainFluid=${PWD}/dummyOF
domainStructure=${PWD}/structureDomain

cd ${domainFluid}
make clean
rm dispCpp.txt
cd ${domainStructure}
rm -r structureResults*
rm -r structureFSISetup/__pycache__
cd ..
rm result_compare.png
rm output.log

# ----------------------------------------------------------------- end-of-file