#! bin/bash


inFile="/work/clas12/rg-m/LH2/prod1.4/dst/recon/015024/rec_clas_015024.evio.00001.hipo";
outFile="test.root"


CMD="clas12root -b \"CLAS12SkimmerTree_simple.C(\\\"${inFile}\\\",\\\"${test.root}\\\")\" "
echo "executing command: $CMD" 
eval ${CMD}
