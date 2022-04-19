#! bin/bash


inFile="/work/clas12/rg-m/LH2/prod1.4/dst/recon/015024/rec_clas_015024.evio.00001.hipo";
outFile="test.root"
target="LH2"

CMD="clas12root -b \"e4nu_main.cpp(\\\"${inFile}\\\",\\\"${outFile}\\\", \\\"${target}\\\" )\" "
echo "executing command: $CMD" 
eval ${CMD}
