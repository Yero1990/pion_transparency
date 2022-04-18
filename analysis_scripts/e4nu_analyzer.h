#ifndef E4NU_ANALYZER_H
#define E4NU_ANALUZER_H

class e4nu_analyzer
{
  
 public:

  // constructor 
  e4nu_analyzer(TString hipo_file="", TString tgt="" );

  // destructor
  ~e4nu_analyzer();

  // main analysis method (calls relevant sub-analysis methods)
  void run_data_analysis();

  // function prototypes
  void SetHistBins();
  void CreateHist();
  void CreateTree();
  void EventLoop();
  void WriteHist();

 protected:
  
  //Set Constants
  const Double_t pi = TMath::Pi(); 
  const Double_t dtr = pi / 180.;

  // init parms
  TString inFile, target;

  // declare TFile pointers (reading/writing ROOTfiles)
  TFile *inROOT;
  TFile *outROOT;

  // particle masses variables 
  Double_t MP, MN, me;
  Double_t pip, pim, pi0;
  Double_t MD;
  
  // declare histogram binning

  // data-related user-defined variables
  TTree *data_tree;
  Long64_t nentries;

  Double_t Eb; 

  // data tree branch variables
  
};

#endif
