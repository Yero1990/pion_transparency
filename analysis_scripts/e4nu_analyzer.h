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
  Double_t me;
  Double_t MP;
  Double_t MN;
  Double_t MD;
  

  // declare generic branch variables
  Double_t br_px, br_py, br_pz, br_p;  // generic final state particle momentum components
  Double_t br_chi2pid;
  Int_t br_pid;    // particle identification (e.g., 11-> electron, 2212-> proton, etc)
  Int_t br_npart; // total number of particles per event

  // declare particles 4-momentum 
  TLorentzVector p4_beam;
  TLorentzVector p4_target;
  TLorentzVector p4_electron;
  TLorentzVector p4_hadron;
  TLorentzVector p4_q; 
  TLorentzVector p4_recoil; // recoil system 4-momenta (usually, undetected)

  // declare electron kinematics
  Double_t Eb;
  Double_t kf, kf_v2, kf_x, kf_y, kf_z;  
  Double_t q, qx, qy, qz;
  Double_t nu;
  Double_t Q2; 
  Double_t xbj;
  Double_t th_e, th_e_v2;
  Double_t th_q, ph_q;
  Double_t W, W2;
  
  // declare hadron kinematics
  Double_t pf, pf_v2, pf_x, pf_y, pf_z;
  Double_t xangle;    // opening angle between scat. e- and detected hadron
  Double_t th_h, th_h_v2;     // hadron in-palne angle
  Double_t MM, MM2;  // missing mass
  Double_t Em;       // missing energy
  Double_t Pm;       // missing momentum
  
  // --- declare histograms ---

  // electron
  TH1F *H_kf;
  TH1F *H_kf_v2;   
  TH1F *H_kfx;
  TH1F *H_kfy;
  TH1F *H_kfz;
  TH1F *H_q;      
  TH1F *H_qx;     
  TH1F *H_qy;     
  TH1F *H_qz;    
  TH1F *H_nu;
  TH1F *H_Q2;     
  TH1F *H_xbj;
  TH1F *H_the;
  TH1F *H_the_v2;
  TH1F *H_thq;    
  TH1F *H_phq;  
  TH1F *H_W;      
  TH1F *H_W2;     

  // hadron
  TH1F *H_pf;
  TH1F *H_pf_v2;
  TH1F *H_pfx;
  TH1F *H_pfy;
  TH1F *H_pfz;
  TH1F *H_thh;
  TH1F *H_thh_v2;
  TH1F *H_MM;     
  TH1F *H_MM2; 

  
};

#endif
