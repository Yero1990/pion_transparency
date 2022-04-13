/*=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=

  Author: Carlos Yero
  Date Started: April 06, 2022
  
  Pion Transparency Analysis Code - CLAS12
  (e4nu Collaboration) 


=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:*/


#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <THStack.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"
#include "clas12writer.h"
#include "HipoChain.h"

using namespace clas12;
using namespace std;

void e4nu_piTransparency(){
  
  //hipo file location: /volatile/clas12/rg-m/LAr/prod1.0/dst/recon/015672/rec_clas_015672.evio.*.hipo  (for liquid Ar, for example)
  // It seems a single run number is divided into many .hipo files which must be combined together.

  TString inFile_data = "";
  

}
