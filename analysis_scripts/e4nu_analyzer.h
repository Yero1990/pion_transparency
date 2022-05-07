#ifndef E4NU_ANALYZER_H
#define E4NU_ANALYZER_H

#include "./UTILS/parse_utils.h" //useful C++ string parsing utilities

class e4nu_analyzer
{
  
 public:

  // constructor 
  e4nu_analyzer(TString inHIPO_fname="", TString outROOT_fname="", TString tgt="", TString det_h="");

  // destructor
  ~e4nu_analyzer();

  // main analysis method (calls relevant sub-analysis methods)
  void run_data_analysis();

  // function prototypes
  void SetParticleMass(); 
  void SetHistBins();
  void SetCuts();
  void CreateHist();
  void CreateTree();
  void EventLoop();
  void WriteHist();

  
 protected:
  
  //Set Constants
  const Double_t pi = TMath::Pi(); 
  const Double_t dtr = pi / 180.;

  // init parms
  TString ifname, ofname, target, det_had;

  // declare data TFile / ROOTTree pointers
  TFile *outROOT;
  TTree *data_tree;
  Long64_t nentries;

  // declare detected particle masses
  Double_t me;
  Double_t MP;
  Double_t MN;
  Double_t MPip;
  Double_t MPim;
  Double_t MPi0;
  Double_t MKp;
  Double_t MKm;
  
  // declare target particle masses
  Double_t MD;
  Double_t MHe4;
  Double_t MC12;
  Double_t MCa40;
  Double_t MCa48;
  Double_t MAr40;
  Double_t MSn120;
  
  // declare generic branch variables
  Double_t br_vx, br_vy, br_vz, br_vt;  // generic final state particle vertex components, time
  Double_t br_px, br_py, br_pz, br_p;  // generic final state particle momentum components
  Double_t br_beta;
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

   //----------DEFINE CUT VARIABLES ---------

  //-------PID cuts-------

  // chi2_pid (for the hadron)
  Bool_t chi2pid_cut_flag;
  Bool_t c_chi2pid;
  Double_t c_chi2pid_min;
  Double_t c_chi2pid_max;
 
  // Z-Reaction Vertex Difference Cut
  Bool_t ztarDiff_cut_flag;
  Bool_t c_ztarDiff;
  Double_t c_ztarDiff_min;
  Double_t c_ztarDiff_max;

  //---------Kinematics Cuts----------

  //4-Momentum Transfers
  Bool_t Q2_cut_flag;
  Bool_t c_Q2;
  Double_t c_Q2_min;
  Double_t c_Q2_max;

  //Missing Energy
  Bool_t Em_cut_flag;
  Bool_t c_Em;
  Double_t c_Em_min;
  Double_t c_Em_max;

  //Invariant Mass, W
  Bool_t W_cut_flag;
  Bool_t c_W;
  Double_t c_W_min;
  Double_t c_W_max;


  // combined cuts variables
  Bool_t c_accpCuts;
  Bool_t c_pidCuts;
  Bool_t c_kinCuts;
  Bool_t c_allCuts;
  
  // declare histogram bining variables

  //-----------------------------
  // PID  Histograms Bins
  //-----------------------------
  Double_t chi2_nbins;
  Double_t chi2_xmin;
  Double_t chi2_xmax;
  
  Double_t beta_nbins;
  Double_t beta_xmin;
  Double_t beta_xmax;

  Double_t ztar_diff_nbins;
  Double_t ztar_diff_xmin;
  Double_t ztar_diff_xmax;

  Double_t zE_nbins;
  Double_t zE_xmin;
  Double_t zE_xmax;
  
  //-----------------------------
  //Kinematics Histograms Bins
  //-----------------------------
  
  // electron Kinematics
  Double_t kf_vert_nbins;   
  Double_t kf_vert_xmin;
  Double_t kf_vert_xmax;

  Double_t kf_vtime_nbins;   
  Double_t kf_vtime_xmin;
  Double_t kf_vtime_xmax;
  
  Double_t kf_nbins;   
  Double_t kf_xmin;
  Double_t kf_xmax;

  Double_t kfx_nbins;   
  Double_t kfx_xmin;
  Double_t kfx_xmax;

  Double_t kfy_nbins;   
  Double_t kfy_xmin;
  Double_t kfy_xmax;

  Double_t kfz_nbins;   
  Double_t kfz_xmin;
  Double_t kfz_xmax;

  Double_t q_nbins;
  Double_t q_xmin;
  Double_t q_xmax;

  Double_t qx_nbins;
  Double_t qx_xmin;
  Double_t qx_xmax;

  Double_t qy_nbins;
  Double_t qy_xmin;
  Double_t qy_xmax;

  Double_t qz_nbins;
  Double_t qz_xmin;
  Double_t qz_xmax;

  Double_t nu_nbins;
  Double_t nu_xmin;
  Double_t nu_xmax;

  Double_t Q2_nbins;
  Double_t Q2_xmin;
  Double_t Q2_xmax;

  Double_t xbj_nbins;
  Double_t xbj_xmin;
  Double_t xbj_xmax;
  
  Double_t the_nbins;
  Double_t the_xmin;
  Double_t the_xmax;

  Double_t phe_nbins;
  Double_t phe_xmin;
  Double_t phe_xmax;
  
  Double_t thq_nbins;
  Double_t thq_xmin;
  Double_t thq_xmax;

  Double_t phq_nbins;
  Double_t phq_xmin;
  Double_t phq_xmax;
  
  Double_t W_nbins;
  Double_t W_xmin;
  Double_t W_xmax;

  Double_t W2_nbins;
  Double_t W2_xmin;
  Double_t W2_xmax;

  // hadron
  Double_t pf_vert_nbins;   
  Double_t pf_vert_xmin;
  Double_t pf_vert_xmax;

  Double_t pf_vtime_nbins;   
  Double_t pf_vtime_xmin;
  Double_t pf_vtime_xmax;
  
  Double_t pf_nbins;
  Double_t pf_xmin;
  Double_t pf_xmax;

  Double_t pfx_nbins;
  Double_t pfx_xmin;
  Double_t pfx_xmax;

  Double_t pfy_nbins;
  Double_t pfy_xmin;
  Double_t pfy_xmax;

  Double_t pfz_nbins;
  Double_t pfz_xmin;
  Double_t pfz_xmax;
  
  Double_t thx_nbins; 
  Double_t thx_xmin;
  Double_t thx_xmax;


  Double_t MM_nbins;
  Double_t MM_xmin;
  Double_t MM_xmax;

  Double_t MM2_nbins;  
  Double_t MM2_xmin;
  Double_t MM2_xmax;

  
  Double_t Em_nbins;
  Double_t Em_xmin;
  Double_t Em_xmax;
  
  Double_t Pm_nbins;
  Double_t Pm_xmin;
  Double_t Pm_xmax;

  Double_t Pmx_lab_nbins;
  Double_t Pmx_lab_xmin;
  Double_t Pmx_lab_xmax;

  Double_t Pmy_lab_nbins;
  Double_t Pmy_lab_xmin;
  Double_t Pmy_lab_xmax;

  Double_t Pmz_lab_nbins;
  Double_t Pmz_lab_xmin;
  Double_t Pmz_lab_xmax;

  Double_t Pmx_q_nbins;
  Double_t Pmx_q_xmin;
  Double_t Pmx_q_xmax;

  Double_t Pmy_q_nbins;
  Double_t Pmy_q_xmin;
  Double_t Pmy_q_xmax;

  Double_t Pmz_q_nbins;
  Double_t Pmz_q_xmin;
  Double_t Pmz_q_xmax;

  Double_t Tx_nbins;
  Double_t Tx_xmin;
  Double_t Tx_xmax;

  Double_t Tr_nbins;
  Double_t Tr_xmin;
  Double_t Tr_xmax;
  
  Double_t thxq_nbins;
  Double_t thxq_xmin;
  Double_t thxq_xmax;

  Double_t thrq_nbins;
  Double_t thrq_xmin;
  Double_t thrq_xmax;

  Double_t phxq_nbins;
  Double_t phxq_xmin;
  Double_t phxq_xmax;

  Double_t phrq_nbins;
  Double_t phrq_xmin;
  Double_t phrq_xmax;

  // declare generic mass variables
  Double_t Mt;  // mass of target
  Double_t Mh;  // mass of detected particle X

  // declare pid variables
  Double_t e_beta;
  Double_t h_beta;
  Double_t h_chi2pid;
  Double_t e_chi2pid;
  
  // declare electron kinematic variables
  Double_t Eb;
  Double_t kf_vx, kf_vy, kf_vz, kf_vt; // e- vertex, time [cm], [ns]
  Double_t kf, kf_x, kf_y, kf_z;  
  Double_t q, qx, qy, qz;
  Double_t nu;
  Double_t Q2; 
  Double_t xbj;
  Double_t th_e;
  Double_t ph_e;
  Double_t th_q, ph_q;
  Double_t W, W2;
 
  
  // declare hadron kinematic variables
  Double_t pf_vx, pf_vy, pf_vz, pf_vt; // hadron vertex, time [cm], [ns]
  Double_t pf, pf_x, pf_y, pf_z;
  Double_t th_x;     // hadron in-palne angle
  Double_t MM, MM2;  // missing mass
  Double_t Em, Em_nuc, Em_recoil;       // missing energy, deuteron special :), production reaction
  Double_t Pm;
  Double_t Pmx_lab, Pmy_lab, Pmz_lab;  // missing momentum (lab frame)
  Double_t Pmx_q, Pmy_q, Pmz_q;  // missing momentum (lab frame)
  Double_t Tx;                    //Kinetic Energy of detected particle (proton, kaon, pion, . . .) [GeV]
  Double_t Tr;                    //Kinetic Energy of recoil system [GeV]
  Double_t th_xq;                 //In-plane angle between detected particle and q [rad]  
  Double_t th_rq;                 //In-plane angle between the recoil system and q [rad]  
  Double_t ph_xq;                 //Out-of-plane angle between detected particle and q [rad]   
  Double_t ph_rq;                 //Out-of-plane anfle between recoil system and q [rad]
 
  //user-defined variables
  Double_t ztar_diff; // difference between the z-vertex component between scattered electron and detected hadron
  Double_t zE;  // energy of the pion divided by transferred energy E_pi / nu (assuming detected hadron is a pion)

  // --- declare histograms ---
  
  //--------------------
  // Acceptance Histos
  //--------------------
  TH2F *H_the_vs_phe;
  
  //--------------------
  // Particle ID Histos
  //--------------------
  TH1F *H_chi2pid_elec;
  TH1F *H_chi2pid_had;
  TH1F *H_beta_elec;
  TH1F *H_beta_had;
  TH2F *H_beta_vs_kf;
  TH2F *H_beta_vs_pf;
  TH1F *H_ztar_diff;
  TH1F *H_zE;
  TH2F *H_pf_vs_thxq;  // detected hadron momentum (pf) vs. hadron angle relative to q-vector (th_xq)
  //--------------------
  // Kinematics Histos
  //-------------------

  
  // electron
  TH1F *H_kf_vx;
  TH1F *H_kf_vy;
  TH1F *H_kf_vz;
  TH1F *H_kf_vt;
  TH1F *H_kf;
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
  TH1F *H_phe;
  TH1F *H_thq;    
  TH1F *H_phq;  
  TH1F *H_W;      
  TH1F *H_W2;     
  
  
  // hadron
  TH1F *H_pf_vx;
  TH1F *H_pf_vy;
  TH1F *H_pf_vz;
  TH1F *H_pf_vt;
  TH1F *H_pf;
  TH1F *H_pfx;
  TH1F *H_pfy;
  TH1F *H_pfz;
  TH1F *H_thx;
  TH1F *H_MM;     
  TH1F *H_MM2; 
  TH1F *H_Em;	  
  TH1F *H_Em_nuc;
  TH1F *H_Em_recoil; 
  TH1F *H_Pm;	  
  TH1F *H_Pmx_lab;
  TH1F *H_Pmy_lab;
  TH1F *H_Pmz_lab;
  TH1F *H_Pmx_q;  
  TH1F *H_Pmy_q;  
  TH1F *H_Pmz_q;  
  TH1F *H_Tx;	  
  TH1F *H_Tr;	  
  TH1F *H_thxq;	  
  TH1F *H_thrq;	  
  TH1F *H_phxq;	  
  TH1F *H_phrq;	  
  

  // 2d kinematics
  TH2F *H_kf_vs_the;
  

  // selected kin. @ Forward Detector 
  TH1F *H_W_FD;
  TH1F *H_W_FD_sec1;
  TH1F *H_W_FD_sec2;
  TH1F *H_W_FD_sec3;
  TH1F *H_W_FD_sec4;
  TH1F *H_W_FD_sec5;
  TH1F *H_W_FD_sec6;

  
  //Create Categorical TLists to store histograms based on caterogy
  TList *accp_HList;    // acceptance histogram list
  TList *pid_HList;     // store PID histograms (i.e., Beta vs. P, etc.)
  TList *kin_HList;     // store kinematical histograms (i.e., Q2, W, th_e, etc.)

  TList *kin_HList_FD; //store kinematical histograms ONLY detected in FD (i.e., Q2, W, th_e, etc.)
  
  
};

#endif
