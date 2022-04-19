#ifndef E4NU_ANALYZER_H
#define E4NU_ANALYZER_H


class e4nu_analyzer
{
  
 public:

  // constructor 
  e4nu_analyzer(TString inHIPO_fname="", TString outROOT_fname="", TString tgt="" );

  // destructor
  ~e4nu_analyzer();

  // main analysis method (calls relevant sub-analysis methods)
  void run_data_analysis();

  // function prototypes
  void SetParticleMass(); 
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
  TString ifname, ofname, target;

  // declare data TFile / ROOTTree pointers
  TFile *outROOT;
  TTree *data_tree;
  Long64_t nentries;

  // declare particle masses
  Double_t MP;
  Double_t me;
  
  // declare particles 4-momentum 
  TLorentzVector p4_beam;
  TLorentzVector p4_target;
  TLorentzVector p4_electron;
  TLorentzVector p4_proton;
  TLorentzVector p4_q; 
  TLorentzVector p4_recoil; // recoil system 4-momenta (usually, undetected)

  // declare branch variables
  Double_t px, py, pz, chi2pid;
  Int_t pid;
  
  // declare histograms
  TH1F *H_the;    
  TH1F *H_kf;     
  TH1F *H_kfx;
  TH1F *H_kfy;
  TH1F *H_kfz;
  TH1F *H_W;      
  TH1F *H_W2;     
  TH1F *H_Q2;     
  TH1F *H_xbj;    
  TH1F *H_nu;     
  TH1F *H_q;      
  TH1F *H_qx;     
  TH1F *H_qy;     
  TH1F *H_qz;     
  TH1F *H_thq;    
  TH1F *H_phq;      	     
  TH1F *H_MM;     
  TH1F *H_MM2; 

  
};

#endif
