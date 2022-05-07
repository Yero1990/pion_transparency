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
echo "[target]: H, D, He4, C12, Ca40, Ca48, Ar40, Sn120"
echo "[detected_hadron]: p, pi+, pi-"
echo "Example: ./run_analysis.sh D pi+"
echo "--------------------------"
echo ""
exit 1
fi

# select characteristis run based on selected target (selected runs have beam energy 5.98636 GeV)
if [[ $target == "H" ]]; then
    inFile="/work/clas12/rg-m/LH2/prod1.4/dst/recon/015024/rec_clas_015024.evio.0000*.hipo";
elif [[ $target == "D" ]]; then
    inFile="/work/clas12/rg-m/LD2/prod1.4/dst/recon/015045/rec_clas_015045.evio.0000*.hipo";
elif [[ $target == "He4" ]]; then
    inFile="/work/clas12/rg-m/LHe/prod1.4/dst/recon/015133/rec_clas_015133.evio.0000*.hipo";
elif [[ $target == "C12" ]]; then
    inFile="/work/clas12/rg-m/Cx4/prod1.4/dst/recon/015188/rec_clas_015188.evio.0000*.hipo";
elif [[ $target == "Ca40" ]]; then
    inFile="/work/clas12/rg-m/40Ca/prod1.4/dst/recon/015392/rec_clas_015392.evio.0000*.hipo";
elif [[ $target == "Ca48" ]]; then
    inFile="/work/clas12/rg-m/48Ca/prod1.4/dst/recon/015832/rec_clas_015832.evio.0000*.hipo";
elif [[ $target == "Ar40" ]]; then
    inFile="/work/clas12/rg-m/LAr/prod1.5/dst/recon/015795/rec_clas_015795.evio.0000*.hipo";
elif [[ $target == "Sn120" ]]; then
    inFile="/work/clas12/rg-m/Snx4/prod1.4/dst/recon/015327/rec_clas_015327.evio.0000*.hipo";
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
