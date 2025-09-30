#!/bin/sh

domainFluid=${PWD}/dummyOF
domainStructure=${PWD}/structureDomain

cd ${domainFluid}
rm -f *.log
rm dispCpp.txt
rm PUSHER_FETCHER_1
cd build make clean
cd ..
rm -fr build
cd ${domainStructure}
rm -r structureResults*
rm -r structureFSISetup/__pycache__
cd ..
rm result_compare.png
rm output.log

# ----------------------------------------------------------------- end-of-file