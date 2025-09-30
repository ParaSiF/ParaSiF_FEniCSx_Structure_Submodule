#!/bin/sh

domainStructure=${PWD}/structureDomain

cd ${domainStructure}
rm -r structureResults*
rm -r structureFSISetup/__pycache__
cd ..
rm result_compare.png
rm output.log

# ----------------------------------------------------------------- end-of-file