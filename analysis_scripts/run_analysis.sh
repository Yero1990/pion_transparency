#!/bin/sh

# This shell script requires a target type and detected hadron
# and assumes a reaction of the type: A(e,e'h), where an electron
# interacts with target A and a electron and knocked-out/produced hadron are
# detected in final state.  eA -> e'h (missing) 

# user input  
target=$1   # target
det_had=$2  # detected hadron

# check if # of arguments is NOT 2
if [ $# -ne 2 ]; then
echo ""
echo "-------------------------"
echo "No Input Arguments Found ! "
echo "Code Usage: ./run_analysis.sh [target] [detected_hadron] "
echo "[target]: H, D"
echo "[detected_hadron]: p, pi+, pi-"
echo "Example: ./run_analysis.sh D pi+"
echo "--------------------------"
echo ""
exit 1
fi

# select characteristis run based on selected target
if [[ $target == "H" ]]; then
    inFile="/work/clas12/rg-m/LH2/prod1.4/dst/recon/015024/rec_clas_015024.evio.00001.hipo";
elif [[ $target == "D" ]]; then
    inFile="/work/clas12/rg-m/LD2/prod1.4/dst/recon/015045/rec_clas_015045.evio.00001.hipo";
fi

# select detected hadron, and set output file name
if [[ $det_had == "p" ]] || [[ $det_had == "pi+" ]] || [[ $det_had == "pi-" ]]; then
    echo "Analyzing Reaction: $target(e,e'$det_had)"
    outFile="rgm_${target}ee${det_had}_output.root"
fi

echo "Input File: $inFile"
echo "Output File: $outFile" 


CMD="clas12root -b \"e4nu_main.cpp(\\\"${inFile}\\\",\\\"${outFile}\\\", \\\"${target}\\\", \\\"${det_had}\\\")\" "


echo "executing command: $CMD" 
eval ${CMD}
