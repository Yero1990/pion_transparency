#!/bin/sh

inFile="/work/clas12/rg-m/LH2/prod1.4/dst/recon/015024/rec_clas_015024.evio.000*.hipo";
outFile="rgm_LH2_015024.root"
target="LH2"
detected="proton"

CMD="clas12root -b \"e4nu_main.cpp(\\\"${inFile}\\\",\\\"${outFile}\\\", \\\"${target}\\\", \\\"{detected}\\\")\" "
echo "executing command: $CMD" 
eval ${CMD}
