/*----------------------------------------------------------------------------------------------

  Author: Carlos Yero
  Date Started: April 06, 2022
  
  e4nu analyzer - CLAS12
  (e4nu Collaboration) 
  
  This code is used for sanity checks 
  and pion transparency analysis of Hall B 
  run-group M (RGM) dataset taken on 
  Fall-2021 -> Spring 2022

  The general reaction assumed is: A(e,e'h)X
  e + A -> e' + h + X
  e: incident electron
  A: target nucleus with A nucleons
  e': scattered electron
  h: knocked-out hadron(s), where:
     quasi-elastic (x ~ 1): h is a knocked-out proton (requires proton to carry most of the momentum transferred)
     inelastic, resonance, or DIS ( x < 1 ): h is a newly produced hadron
  X: recoil or residual system ("missing particle") reconstructed from momentum conservation
     of the detected particles

----------------------------------------------------------------------------------------------*/


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
// include this analyzer header file
#include "e4nu_analyzer.h"

using namespace clas12;
using namespace std;

//_______________________________________________________________________________
e4nu_analyzer::e4nu_analyzer(TString inHIPO_fname="", TString outROOT_fname="", TString tgt="", TString det_h="" )
  : ifname(inHIPO_fname), ofname(outROOT_fname), target(tgt), det_had(det_h)
{
  
  
  cout << "Start e4nu_analyzer::e4nu_analyzer ... " << endl;
  
  /* hipo file locations:
     "/volatile/clas12/rg-m/LAr/prod1.0/dst/recon/015672/rec_clas_015672.evio.*.hipo"  
     "/work/clas12/rg-m/LH2/prod1.4/dst/recon/015024/rec_clas_015024.evio.0001*.hipo"
  */
  
  /*
    //Record analysis time
    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    time_elapsed = diff.count();
  */
  
  // load CLAS12 databases / Chain multiple .hipo files
  // the user can give multiple .hipo file inputs to be chained,
  // i.e., path/to/run/*.hipo (all hipo files of a given run)
  clas12databases::SetCCDBLocalConnection("ccdb.sqlite");
  clas12databases::SetRCDBRootConnection("rcdb.root");

  SetParticleMass();
  
  //Initialize TFile Pointers
  outROOT = NULL;
  
  //Initialize TTree Pointers
  data_tree = NULL;

  //Initialize TList Pointers
  accp_HList = NULL;
  pid_HList = NULL;
  kin_HList = NULL;
  kin_HList_FD = NULL;
  
  // --- initialize histogram pointers ----

  //--------------------
  // Acceptance Histos
  //--------------------
  H_the_vs_phe = NULL;
  
  //--------------------
  // Particle ID Histos
  //--------------------
  H_chi2pid_elec = NULL;
  H_chi2pid_had = NULL;
  H_beta_elec = NULL;
  H_beta_had   = NULL;
  H_beta_vs_kf = NULL;
  H_beta_vs_pf = NULL;
  H_ztar_diff = NULL;
  H_zE = NULL;
  H_pf_vs_thxq = NULL;
  
  //--------------------
  // Kinematics Histos
  //-------------------
  
  // electron kinematics
  H_kf_vx = NULL;
  H_kf_vy = NULL;
  H_kf_vz = NULL;
  H_kf_vt = NULL;
  H_kf   =  NULL;
  H_kfx  =  NULL; 
  H_kfy  =  NULL; 
  H_kfz  =  NULL;
  H_q    =  NULL;
  H_qx 	 =  NULL; 
  H_qy 	 =  NULL; 
  H_qz 	 =  NULL; 
  H_nu   =  NULL;
  H_Q2   =  NULL; 
  H_xbj  =  NULL;
  H_the  =  NULL;
  H_phe  =  NULL;
  H_thq	 =  NULL; 
  H_phq	 =  NULL;  
  H_W    =  NULL; 
  H_W2   =  NULL; 
  
  
  // hadron kinematics
  H_pf_vx = NULL;
  H_pf_vy = NULL;
  H_pf_vz = NULL;
  H_pf_vt = NULL;
  H_pf   =  NULL;
  H_pfx  =  NULL;
  H_pfy  =  NULL;
  H_pfz  =  NULL;
  H_thx  =  NULL;
  H_MM   =  NULL; 
  H_MM2  =  NULL; 
  H_Em   =  NULL;
  H_Em_nuc   =  NULL;
  H_Em_recoil   =  NULL;
  H_Pm       =  NULL;
  H_Pmx_lab  =  NULL;
  H_Pmy_lab  =  NULL;
  H_Pmz_lab  =  NULL;
  H_Pmx_q    =  NULL;  
  H_Pmy_q    =  NULL;  
  H_Pmz_q    =  NULL;  
  H_Tx       =  NULL;	  
  H_Tr       =  NULL;	  
  H_thxq       = NULL;	
  H_thrq       = NULL;	
  H_phxq       = NULL;	
  H_phrq       = NULL;
 
  // 2D kinematics
  H_kf_vs_the  = NULL;
 

  // selected kin. @ Forward Detector 
  H_W_FD        = NULL;
  H_W_FD_sec1   = NULL;
  H_W_FD_sec2   = NULL;
  H_W_FD_sec3   = NULL;
  H_W_FD_sec4   = NULL;
  H_W_FD_sec5   = NULL;
  H_W_FD_sec6   = NULL; 
  
}

//_______________________________________________________________________________
e4nu_analyzer::~e4nu_analyzer()
{

    cout << "Start e4nu_analyzer::~e4nu_analyzer ... " << endl;

    // delete TFile Pointers
    delete outROOT;   outROOT = NULL;
    
    
    // delete TTree Pointers
    delete data_tree; data_tree = NULL;

    //Delete TList Pointers
    delete accp_HList; accp_HList = NULL;
    delete pid_HList;  pid_HList  = NULL;
    delete kin_HList;  kin_HList  = NULL;

    delete kin_HList_FD;  kin_HList_FD  = NULL;

    // --- delete histogram pointers ---

    //--------------------
    // Acceptance Histos
    //--------------------
    delete H_the_vs_phe; H_the_vs_phe = NULL;
    
    //--------------------
    // Particle ID Histos
    //--------------------
    delete H_chi2pid_elec; H_chi2pid_elec = NULL;
    delete H_chi2pid_had; H_chi2pid_had = NULL;
    delete H_beta_elec; H_beta_elec = NULL;
    delete H_beta_had; H_beta_had = NULL;
    delete H_beta_vs_kf; H_beta_vs_kf = NULL;
    delete H_beta_vs_pf; H_beta_vs_pf = NULL;
    delete H_ztar_diff; H_ztar_diff = NULL;
    delete H_zE; H_zE = NULL;
    delete H_pf_vs_thxq; H_pf_vs_thxq = NULL;
    
    //--------------------
    // Kinematics Histos
    //-------------------
    // electron
    delete H_kf_vx ; H_kf_vx = NULL;
    delete H_kf_vy ; H_kf_vy = NULL;
    delete H_kf_vz ; H_kf_vz = NULL;
    delete H_kf_vt ; H_kf_vt = NULL;
    
    delete H_kf  ;  H_kf   =  NULL;
    delete H_kfx ;  H_kfx   =  NULL; 
    delete H_kfy ;  H_kfy   =  NULL; 
    delete H_kfz ;  H_kfz   =  NULL; 
    delete H_q   ;  H_q    =  NULL; 
    delete H_qx  ;  H_qx   =  NULL; 
    delete H_qy  ;  H_qy   =  NULL; 
    delete H_qz  ;  H_qz   =  NULL; 
    delete H_nu  ;  H_nu   =  NULL;
    delete H_Q2  ;  H_Q2   =  NULL; 
    delete H_xbj ;  H_xbj  =  NULL;
    delete H_the ;  H_the  =  NULL;
    delete H_phe ;  H_phe  =  NULL;
    delete H_thq ;  H_thq  =  NULL; 
    delete H_phq ;  H_phq  =  NULL; 
    delete H_W   ;  H_W    =  NULL; 
    delete H_W2  ;  H_W2   =  NULL;     
    
    // hadron
    delete H_pf_vx ; H_pf_vx = NULL;
    delete H_pf_vy ; H_pf_vy = NULL;
    delete H_pf_vz ; H_pf_vz = NULL;
    delete H_pf_vt ; H_pf_vt = NULL;
    delete H_pf  ;  H_pf   =  NULL;
    delete H_pfx  ;  H_pfx   =  NULL;
    delete H_pfy  ;  H_pfy   =  NULL;
    delete H_pfz  ;  H_pfz   =  NULL;
    delete H_thx  ;  H_thx   = NULL;
    delete H_MM  ;  H_MM   =  NULL; 
    delete H_MM2 ;  H_MM2  =  NULL; 
    delete H_Em      ;	  H_Em       =  NULL;	
    delete H_Em_nuc  ;	  H_Em_nuc   =  NULL;
    delete H_Em_recoil  ;	  H_Em_recoil   =  NULL;	
    delete H_Pm      ;	  H_Pm       =  NULL;	
    delete H_Pmx_lab ;	  H_Pmx_lab  =  NULL;	
    delete H_Pmy_lab ;	  H_Pmy_lab  =  NULL;	
    delete H_Pmz_lab ;	  H_Pmz_lab  =  NULL;	
    delete H_Pmx_q   ;    H_Pmx_q    =  NULL; 
    delete H_Pmy_q   ;    H_Pmy_q    =  NULL; 
    delete H_Pmz_q   ;    H_Pmz_q    =  NULL; 
    delete H_Tx      ;	  H_Tx       =  NULL;	
    delete H_Tr      ;	  H_Tr       =  NULL;	
    delete H_thxq;	   H_thxq       = NULL;	
    delete H_thrq;	   H_thrq       = NULL;	
    delete H_phxq;	   H_phxq       = NULL;	
    delete H_phrq;	   H_phrq       = NULL;	
    

    // 2D kinematics
   
    delete H_kf_vs_the;  H_kf_vs_the  = NULL;
   

    // selected kin. @ Forward Detector 
    delete H_W_FD;      H_W_FD       = NULL;
    delete H_W_FD_sec1; H_W_FD_sec1  = NULL;
    delete H_W_FD_sec2; H_W_FD_sec2  = NULL;
    delete H_W_FD_sec3; H_W_FD_sec3  = NULL;
    delete H_W_FD_sec4; H_W_FD_sec4  = NULL;
    delete H_W_FD_sec5; H_W_FD_sec5  = NULL;
    delete H_W_FD_sec6; H_W_FD_sec6  = NULL;
    

  
}

//_______________________________________________________________________________
void e4nu_analyzer::SetParticleMass()
{
  
  cout << "Start e4nu_analyzer::SetParticleMass() ... " << endl;

  // load particle database from PDG to get particle mass (GeV)
  auto db=TDatabasePDG::Instance();
  Double_t amu2GeV = 0.931494;  // GeV
    
  // detected leptons / hadrons (some of these might not be used at all)
  me   = db->GetParticle(11)  ->Mass();  // electron
  MP   = db->GetParticle(2212)->Mass();  // proton
  MN   = db->GetParticle(2112)->Mass();  // neutron
  MPip = db->GetParticle(211) ->Mass();  // pi+
  MPim = db->GetParticle(-211)->Mass();  // pi-
  MPi0 = db->GetParticle(111) ->Mass();  // pi0
  MKp  = db->GetParticle(321) ->Mass();  // K+
  MKm  = db->GetParticle(-321)->Mass();  // K-

  // target mass (amu -> GeV) (obtained from NIST: https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses)
  MD     = 2.01410177812 * amu2GeV ;
  MHe4   = 4.00260325413 * amu2GeV ;
  MC12   = 12.0000000    * amu2GeV ;
  MCa40  = 39.962590863  * amu2GeV ;
  MCa48  = 47.95252276   * amu2GeV ;
  MAr40  = 39.9623831237 * amu2GeV ; 
  MSn120 = 119.90220163  * amu2GeV ;

  // set target mass
  if      (target=="H")     { Mt = MP     ;} 
  else if (target=="D")     { Mt = MD     ;}
  else if (target=="He4")   { Mt = MHe4   ;}
  else if (target=="C12")   { Mt = MC12   ;}
  else if (target=="Ca40")  { Mt = MCa40  ;}
  else if (target=="Ca48")  { Mt = MCa48  ;}
  else if (target=="Ar40")  { Mt = MAr40  ;}
  else if (target=="Sn120") { Mt = MSn120 ;}
  else {
    cout << "Please enter one of the targets: \n "
      "H, D, He4, C12, Ca40, Ca48, Ar40, Sn120 " << endl;
    exit(0);
  }
  
  // set detected (X) primary hadron mass A(e,e'X)r
  if(det_had=="p") { Mh = MP; }
  else if(det_had=="pi+")    { Mh = MPip; }
  else if(det_had=="pi-")    { Mh = MPim; }
  else if(det_had=="pi0")    { Mh = MPi0; }
  else {
    cout << "Please enter detected hadron: \n "
      "p, pi+, pi-, pi0 " << endl;
    exit(0);
  }

  cout << Form("Analyzing Reaction: %s(e,e'%s)", target.Data(), det_had.Data()) << endl;
  
}

//_______________________________________________________________________________
void e4nu_analyzer::SetHistBins()
{
  cout << "Start e4nu_analyzer::SetHistBins() ... " << endl;
  
  TString input_HBinFileName = "histo_bins.inp";
  
  //-------------------------
  // PID Histograms Binning
  //-------------------------
  
  chi2_nbins    = stod(split(FindString("chi2_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  chi2_xmin     = stod(split(FindString("chi2_xmin",  input_HBinFileName.Data())[0], '=')[1]);  
  chi2_xmax     = stod(split(FindString("chi2_xmax",  input_HBinFileName.Data())[0], '=')[1]);  

  beta_nbins    = stod(split(FindString("beta_nbins",  input_HBinFileName.Data())[0], '=')[1]);                                                                                                 
  beta_xmin     = stod(split(FindString("beta_xmin",  input_HBinFileName.Data())[0], '=')[1]);                                                                                                     
  beta_xmax     = stod(split(FindString("beta_xmax",  input_HBinFileName.Data())[0], '=')[1]);

  ztar_diff_nbins	= stod(split(FindString("ztar_diff_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  ztar_diff_xmin 	= stod(split(FindString("ztar_diff_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  ztar_diff_xmax 	= stod(split(FindString("ztar_diff_xmax",  input_HBinFileName.Data())[0], '=')[1]);
  
  zE_nbins	= stod(split(FindString("zE_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  zE_xmin 	= stod(split(FindString("zE_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  zE_xmax 	= stod(split(FindString("zE_xmax",  input_HBinFileName.Data())[0], '=')[1]);

  //---------------------------------
  // Kinematics Histograms Binning
  //---------------------------------

  // detected electron
  kf_vert_nbins	= stod(split(FindString("kf_vert_nbins",  	input_HBinFileName.Data())[0], '=')[1]);  
  kf_vert_xmin	= stod(split(FindString("kf_vert_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  kf_vert_xmax	= stod(split(FindString("kf_vert_xmax",  	input_HBinFileName.Data())[0], '=')[1]);

  kf_vtime_nbins = stod(split(FindString("kf_vtime_nbins",  	input_HBinFileName.Data())[0], '=')[1]);  
  kf_vtime_xmin	 = stod(split(FindString("kf_vtime_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  kf_vtime_xmax	 = stod(split(FindString("kf_vtime_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  
  kf_nbins	= stod(split(FindString("kf_nbins",  	input_HBinFileName.Data())[0], '=')[1]);  
  kf_xmin	= stod(split(FindString("kf_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  kf_xmax	= stod(split(FindString("kf_xmax",  	input_HBinFileName.Data())[0], '=')[1]);

  kfx_nbins	= stod(split(FindString("kfx_nbins",  	input_HBinFileName.Data())[0], '=')[1]);  
  kfx_xmin	= stod(split(FindString("kfx_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  kfx_xmax	= stod(split(FindString("kfx_xmax",  	input_HBinFileName.Data())[0], '=')[1]);

  kfy_nbins	= stod(split(FindString("kfy_nbins",  	input_HBinFileName.Data())[0], '=')[1]);  
  kfy_xmin	= stod(split(FindString("kfy_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  kfy_xmax	= stod(split(FindString("kfy_xmax",  	input_HBinFileName.Data())[0], '=')[1]);

  kfz_nbins	= stod(split(FindString("kfz_nbins",  	input_HBinFileName.Data())[0], '=')[1]);  
  kfz_xmin	= stod(split(FindString("kfz_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  kfz_xmax	= stod(split(FindString("kfz_xmax",  	input_HBinFileName.Data())[0], '=')[1]);

  q_nbins      	= stod(split(FindString("q_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  q_xmin       	= stod(split(FindString("q_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  q_xmax       	= stod(split(FindString("q_xmax",  input_HBinFileName.Data())[0], '=')[1]);
               				          
  qx_nbins     	= stod(split(FindString("qx_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  qx_xmin      	= stod(split(FindString("qx_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  qx_xmax      	= stod(split(FindString("qx_xmax",  input_HBinFileName.Data())[0], '=')[1]);
               				          
  qy_nbins     	= stod(split(FindString("qy_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  qy_xmin      	= stod(split(FindString("qy_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  qy_xmax      	= stod(split(FindString("qy_xmax",  input_HBinFileName.Data())[0], '=')[1]);
               				          
  qz_nbins     	= stod(split(FindString("qz_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  qz_xmin      	= stod(split(FindString("qz_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  qz_xmax      	= stod(split(FindString("qz_xmax",  input_HBinFileName.Data())[0], '=')[1]);
  
  nu_nbins     	= stod(split(FindString("nu_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  nu_xmin      	= stod(split(FindString("nu_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  nu_xmax      	= stod(split(FindString("nu_xmax",  input_HBinFileName.Data())[0], '=')[1]);
  
  Q2_nbins     	= stod(split(FindString("Q2_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  Q2_xmin      	= stod(split(FindString("Q2_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  Q2_xmax      	= stod(split(FindString("Q2_xmax",  input_HBinFileName.Data())[0], '=')[1]);
               				          
  xbj_nbins     = stod(split(FindString("xbj_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  xbj_xmin      = stod(split(FindString("xbj_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  xbj_xmax      = stod(split(FindString("xbj_xmax",  input_HBinFileName.Data())[0], '=')[1]);
  
  the_nbins    	= stod(split(FindString("the_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  the_xmin     	= stod(split(FindString("the_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  the_xmax     	= stod(split(FindString("the_xmax",  input_HBinFileName.Data())[0], '=')[1]);

  phe_nbins    	= stod(split(FindString("phe_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  phe_xmin     	= stod(split(FindString("phe_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  phe_xmax     	= stod(split(FindString("phe_xmax",  input_HBinFileName.Data())[0], '=')[1]);

  thq_nbins    	= stod(split(FindString("thq_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  thq_xmin     	= stod(split(FindString("thq_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  thq_xmax     	= stod(split(FindString("thq_xmax",  input_HBinFileName.Data())[0], '=')[1]);
               				          
  phq_nbins    	= stod(split(FindString("phq_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  phq_xmin     	= stod(split(FindString("phq_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  phq_xmax     	= stod(split(FindString("phq_xmax",  input_HBinFileName.Data())[0], '=')[1]);
     
  W_nbins      	= stod(split(FindString("W_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  W_xmin       	= stod(split(FindString("W_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  W_xmax       	= stod(split(FindString("W_xmax",  input_HBinFileName.Data())[0], '=')[1]);
               				          
  W2_nbins     	= stod(split(FindString("W2_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  W2_xmin      	= stod(split(FindString("W2_xmin",  input_HBinFileName.Data())[0], '=')[1]);
  W2_xmax      	= stod(split(FindString("W2_xmax",  input_HBinFileName.Data())[0], '=')[1]);


  // detected hadron kinematics (and also recoil system)
  pf_vert_nbins	= stod(split(FindString("pf_vert_nbins",  	input_HBinFileName.Data())[0], '=')[1]);  
  pf_vert_xmin	= stod(split(FindString("pf_vert_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  pf_vert_xmax	= stod(split(FindString("pf_vert_xmax",  	input_HBinFileName.Data())[0], '=')[1]);

  pf_vtime_nbins = stod(split(FindString("pf_vtime_nbins",  	input_HBinFileName.Data())[0], '=')[1]);  
  pf_vtime_xmin	 = stod(split(FindString("pf_vtime_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  pf_vtime_xmax	 = stod(split(FindString("pf_vtime_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  
  pf_nbins	= stod(split(FindString("pf_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  pf_xmin	= stod(split(FindString("pf_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  pf_xmax	= stod(split(FindString("pf_xmax",  	input_HBinFileName.Data())[0], '=')[1]);

  pfx_nbins	= stod(split(FindString("pfx_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  pfx_xmin	= stod(split(FindString("pfx_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  pfx_xmax	= stod(split(FindString("pfx_xmax",  	input_HBinFileName.Data())[0], '=')[1]);

  pfy_nbins	= stod(split(FindString("pfy_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  pfy_xmin	= stod(split(FindString("pfy_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  pfy_xmax	= stod(split(FindString("pfy_xmax",  	input_HBinFileName.Data())[0], '=')[1]);

  pfz_nbins	= stod(split(FindString("pfz_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  pfz_xmin	= stod(split(FindString("pfz_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  pfz_xmax	= stod(split(FindString("pfz_xmax",  	input_HBinFileName.Data())[0], '=')[1]);

  thx_nbins	= stod(split(FindString("thx_nbins",  	input_HBinFileName.Data())[0], '=')[1]);  
  thx_xmin 	= stod(split(FindString("thx_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  thx_xmax 	= stod(split(FindString("thx_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
     

  MM_nbins     	 = stod(split(FindString("MM_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  MM_xmin      	 = stod(split(FindString("MM_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  MM_xmax      	 = stod(split(FindString("MM_xmax",  	input_HBinFileName.Data())[0], '=')[1]);

  MM2_nbins	= stod(split(FindString("MM2_nbins",  	input_HBinFileName.Data())[0], '=')[1]);  
  MM2_xmin 	= stod(split(FindString("MM2_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  MM2_xmax 	= stod(split(FindString("MM2_xmax",  	input_HBinFileName.Data())[0], '=')[1]);

  Em_nbins     	 = stod(split(FindString("Em_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  Em_xmin      	 = stod(split(FindString("Em_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Em_xmax      	 = stod(split(FindString("Em_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  Pm_nbins     	 = stod(split(FindString("Pm_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  Pm_xmin      	 = stod(split(FindString("Pm_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Pm_xmax      	 = stod(split(FindString("Pm_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
               				               
  Pmx_lab_nbins	 = stod(split(FindString("Pmx_lab_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  Pmx_lab_xmin 	 = stod(split(FindString("Pmx_lab_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Pmx_lab_xmax 	 = stod(split(FindString("Pmx_lab_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  Pmy_lab_nbins	 = stod(split(FindString("Pmy_lab_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  Pmy_lab_xmin 	 = stod(split(FindString("Pmy_lab_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Pmy_lab_xmax 	 = stod(split(FindString("Pmy_lab_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  Pmz_lab_nbins	 = stod(split(FindString("Pmz_lab_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  Pmz_lab_xmin 	 = stod(split(FindString("Pmz_lab_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Pmz_lab_xmax 	 = stod(split(FindString("Pmz_lab_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  Pmx_q_nbins  	 = stod(split(FindString("Pmx_q_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  Pmx_q_xmin   	 = stod(split(FindString("Pmx_q_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Pmx_q_xmax   	 = stod(split(FindString("Pmx_q_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  Pmy_q_nbins  	 = stod(split(FindString("Pmy_q_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  Pmy_q_xmin   	 = stod(split(FindString("Pmy_q_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Pmy_q_xmax   	 = stod(split(FindString("Pmy_q_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  Pmz_q_nbins  	 = stod(split(FindString("Pmz_q_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  Pmz_q_xmin   	 = stod(split(FindString("Pmz_q_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Pmz_q_xmax   	 = stod(split(FindString("Pmz_q_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  Tx_nbins     	 = stod(split(FindString("Tx_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  Tx_xmin      	 = stod(split(FindString("Tx_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Tx_xmax      	 = stod(split(FindString("Tx_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  Tr_nbins     	 = stod(split(FindString("Tr_nbins",  	input_HBinFileName.Data())[0], '=')[1]);
  Tr_xmin      	 = stod(split(FindString("Tr_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  Tr_xmax      	 = stod(split(FindString("Tr_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
	       				  	       
  thxq_nbins   	 = stod(split(FindString("thxq_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  thxq_xmin    	 = stod(split(FindString("thxq_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  thxq_xmax    	 = stod(split(FindString("thxq_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  thrq_nbins   	 = stod(split(FindString("thrq_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  thrq_xmin    	 = stod(split(FindString("thrq_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  thrq_xmax    	 = stod(split(FindString("thrq_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
               				               
  phxq_nbins   	 = stod(split(FindString("phxq_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  phxq_xmin    	 = stod(split(FindString("phxq_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  phxq_xmax    	 = stod(split(FindString("phxq_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  	       				  	       
  phrq_nbins   	 = stod(split(FindString("phrq_nbins",  input_HBinFileName.Data())[0], '=')[1]);
  phrq_xmin    	 = stod(split(FindString("phrq_xmin",  	input_HBinFileName.Data())[0], '=')[1]);
  phrq_xmax    	 = stod(split(FindString("phrq_xmax",  	input_HBinFileName.Data())[0], '=')[1]);
  
    
}

//_______________________________________________________________________________
void e4nu_analyzer::SetCuts()
{

  TString input_CutFileName = "histo_cuts.inp";

  cout << "Start e4nu_analyzer::SetCuts() ... " << endl;

  //------PID Cuts-----

  // chi2pid Cut
  chi2pid_cut_flag = stoi(split(FindString("chi2pid_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_chi2pid_min = stod(split(FindString("c_chi2pid_min", input_CutFileName.Data())[0], '=')[1]);
  c_chi2pid_max = stod(split(FindString("c_chi2pid_max", input_CutFileName.Data())[0], '=')[1]);
  
  // Z-Reaction Vertex Difference Cut
  ztarDiff_cut_flag = stoi(split(FindString("ztarDiff_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_ztarDiff_min = stod(split(FindString("c_ztarDiff_min", input_CutFileName.Data())[0], '=')[1]);
  c_ztarDiff_max = stod(split(FindString("c_ztarDiff_max", input_CutFileName.Data())[0], '=')[1]);

  // z = E_had/nu (fraction of energy transferred to hadron) Cut
  zE_cut_flag = stoi(split(FindString("zE_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_zE_min = stod(split(FindString("c_zE_min", input_CutFileName.Data())[0], '=')[1]);
  c_zE_max = stod(split(FindString("c_zE_max", input_CutFileName.Data())[0], '=')[1]);
  
  //-----Kinematics Cuts------
  
  //4-Momentum Transfers [GeV^2]
  Q2_cut_flag = stoi(split(FindString("Q2_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_Q2_min = stod(split(FindString("c_Q2_min", input_CutFileName.Data())[0], '=')[1]);
  c_Q2_max = stod(split(FindString("c_Q2_max", input_CutFileName.Data())[0], '=')[1]);

  //Missing Energy [GeV]
  Em_cut_flag = stoi(split(FindString("Em_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_Em_min = stod(split(FindString("c_Em_min", input_CutFileName.Data())[0], '=')[1]);
  c_Em_max = stod(split(FindString("c_Em_max", input_CutFileName.Data())[0], '=')[1]);

  //Invariant Mass, W [GeV]
  W_cut_flag = stoi(split(FindString("W_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_W_min = stod(split(FindString("c_W_min", input_CutFileName.Data())[0], '=')[1]);
  c_W_max = stod(split(FindString("c_W_max", input_CutFileName.Data())[0], '=')[1]);

  //Missing Mass, MM [GeV]
  MM_cut_flag = stoi(split(FindString("MM_cut_flag", input_CutFileName.Data())[0], '=')[1]);
  c_MM_min = stod(split(FindString("c_MM_min", input_CutFileName.Data())[0], '=')[1]);
  c_MM_max = stod(split(FindString("c_MM_max", input_CutFileName.Data())[0], '=')[1]);
}

//_______________________________________________________________________________
void e4nu_analyzer::CreateHist()
{
  
  cout << "Start e4nu_analyzer::CreateHist() ... " << endl;

  SetHistBins();

  //Create TLists to store categorical histograms
  accp_HList  = new TList();
  pid_HList  = new TList();
  kin_HList  = new TList();
  kin_HList_FD  = new TList();

  //--------------------
  // Acceptance Histos
  //--------------------
  H_the_vs_phe = new TH2F("H_the_vs_phe", "e^{-} #theta_{e} vs. # phi_{e}; #phi_{e} [deg]; #theta_{e} [deg]", phe_nbins, phe_xmin, phe_xmax, the_nbins, the_xmin, the_xmax);      

  //--------------------
  // Particle ID Histos
  //--------------------
  H_chi2pid_elec = new TH1F("H_chi2pid_elec", "e^{-} #chi^{2} PID", chi2_nbins, chi2_xmin, chi2_xmax);
  H_chi2pid_had = new TH1F("H_chi2pid_had", Form("%s #chi^{2} PID", det_had.Data()), chi2_nbins, chi2_xmin, chi2_xmax);

  H_beta_elec = new TH1F("H_beta_elec",  "e^{-} #beta",                                                     beta_nbins, beta_xmin, beta_xmax);  
  H_beta_had  = new TH1F("H_beta_had",  Form("%s #beta", det_had.Data()),                                   beta_nbins, beta_xmin, beta_xmax); 			     	           
  H_beta_vs_kf = new TH2F("H_beta_vs_kf", " #beta_{e} vs. e^{-} Momentum; k_{f} [GeV/c]; #beta_{e}", kf_nbins, kf_xmin, kf_xmax, beta_nbins, beta_xmin, beta_xmax);      
  H_beta_vs_pf = new TH2F("H_beta_vs_pf", " #beta_{h} vs. Hadron Momentum; p_{f} [GeV/c]; #beta_{h}", pf_nbins, pf_xmin, pf_xmax, beta_nbins, beta_xmin, beta_xmax);

  //difference in reaction vertex z (user-defined)
  H_ztar_diff = new TH1F("H_ztar_diff", "Ztar Difference; z-Target Difference [cm]; Counts", ztar_diff_nbins, ztar_diff_xmin, ztar_diff_xmax);

  // fraction of the total energy transferred to the pion Epi / nu
  H_zE = new TH1F("H_zE", Form("Energy Transferred to %s, E_{h}/#nu; E_{h}/#nu; Counts", det_had.Data()), zE_nbins, zE_xmin, zE_xmax);
  H_pf_vs_thxq = new TH2F("H_pf_vs_thxq", Form("%s p_{f} vs. #theta_{xq}", det_had.Data()), thxq_nbins, thxq_xmin, thxq_xmax,  pf_nbins, pf_xmin, pf_xmax);

  //--------------------
  // Kinematics Histos
  //-------------------
  
  // electron
  H_kf_vx   = new TH1F("H_kf_vx",       "e^{-} x-vertex ",                                    kf_vert_nbins,  kf_vert_xmin,  kf_vert_xmax  );
  H_kf_vy   = new TH1F("H_kf_vy",       "e^{-} y-vertex ",                                    kf_vert_nbins,  kf_vert_xmin,  kf_vert_xmax  );
  H_kf_vz   = new TH1F("H_kf_vz",       "e^{-} z-vertex ",                                    kf_vert_nbins,  kf_vert_xmin,  kf_vert_xmax  );
  H_kf_vt   = new TH1F("H_kf_vt",       "e^{-} time @ Vertex ",                               kf_vtime_nbins, kf_vtime_xmin, kf_vtime_xmax );
  H_kf      = new TH1F("H_kf",          "e^{-} momentum, |#vec{k}_{f}|",                      kf_nbins,  kf_xmin,  kf_xmax  );
  H_kfx     = new TH1F("H_kfx",         "|#vec{k}_{f,x}|",                                    kfx_nbins, kfx_xmin, kfx_xmax );
  H_kfy     = new TH1F("H_kfy",         "|#vec{k}_{f,y}|",                                    kfy_nbins, kfy_xmin, kfy_xmax );			    
  H_kfz     = new TH1F("H_kfz",         "|#vec{k}_{f,z}|",                                    kfz_nbins, kfz_xmin, kfz_xmax );
  H_q       = new TH1F("H_q",           "3-momentum transfer, |#vec{q}|",                     q_nbins,   q_xmin,   q_xmax   );		    
  H_qx      = new TH1F("H_qx",          "|#vec{q}_{x}|",                                      qx_nbins,  qx_xmin,  qx_xmax  );				    
  H_qy      = new TH1F("H_qy",          "|#vec{q}_{y}|",                                      qy_nbins,  qy_xmin,  qy_xmax  );				    
  H_qz      = new TH1F("H_qz",          "|#vec{q}_{z}|",                                      qz_nbins,  qz_xmin,  qz_xmax  );
  H_nu      = new TH1F("H_nu",          "energy transfer, #nu",                               nu_nbins,  nu_xmin,  nu_xmax  );
  H_Q2      = new TH1F("H_Q2",          "4-momentum transfer, Q^{2}",                         Q2_nbins,  Q2_xmin,  Q2_xmax  ); 		    
  H_xbj     = new TH1F("H_xbj",         "x-Bjorken",                                          xbj_nbins, xbj_xmin, xbj_xmax );
  H_the     = new TH1F("H_the",         "e^{-} in-plane scat. angle, #theta_{e}",             the_nbins, the_xmin, the_xmax );
  H_phe     = new TH1F("H_phe",         "e^{-} out-of-plane scat. angle, #phi_{e}",           phe_nbins, phe_xmin, phe_xmax );
  H_thq     = new TH1F("H_thq",         "#vec{q} in-plane angle w.r.t +z(lab), #theta_{q}",   thq_nbins, thq_xmin, thq_xmax );     
  H_phq     = new TH1F("H_phq",         "#vec{q} out-of-plane angle w.r.t +z(lab), #phi_{q}", phq_nbins, phq_xmin, phq_xmax );
  H_W       = new TH1F("H_W",           "invariant mass, W",                                  W_nbins,   W_xmin,   W_xmax   ); 				    
  H_W2      = new TH1F("H_W2",          "invariant mass, W^{2}",                              W2_nbins,  W2_xmin,  W2_xmax  );   
  
  // hadron
  H_pf_vx     = new TH1F("H_pf_vx",     Form(" %s x-vertex ",      det_had.Data()),                         pf_vert_nbins,  pf_vert_xmin,  pf_vert_xmax );
  H_pf_vy     = new TH1F("H_pf_vy",     Form(" %s y-vertex ",      det_had.Data()),                         pf_vert_nbins,  pf_vert_xmin,  pf_vert_xmax );
  H_pf_vz     = new TH1F("H_pf_vz",     Form(" %s z-vertex ",      det_had.Data()),                         pf_vert_nbins,  pf_vert_xmin,  pf_vert_xmax );
  H_pf_vt     = new TH1F("H_pf_vt",     Form(" %s time @ vertex ", det_had.Data()),                         pf_vtime_nbins, pf_vtime_xmin, pf_vtime_xmax );
  H_pf        = new TH1F("H_pf",        Form(" %s momentum (detected), p_{f}", det_had.Data()),             pf_nbins, pf_xmin, pf_xmax  );
  H_pfx       = new TH1F("H_pfx",       Form(" %s momentum, X-comp. p_{fx}",   det_had.Data()),             pfx_nbins, pfx_xmin, pfx_xmax );
  H_pfy       = new TH1F("H_pfy",       Form(" %s momentum, Y-comp. p_{fy}",   det_had.Data()),             pfy_nbins, pfy_xmin, pfy_xmax);
  H_pfz       = new TH1F("H_pfz",       Form(" %s momentum, Z-comp. p_{fz}",   det_had.Data()),             pfz_nbins, pfz_xmin, pfz_xmax  );
  H_thx       = new TH1F("H_thx",       Form(" %s scat. angle (detected), #theta_{x}",  det_had.Data() ),   thx_nbins, thx_xmin, thx_xmax);
  H_MM        = new TH1F("H_MM",        "missing mass, M_{miss}",                                           MM_nbins, MM_xmin, MM_xmax);        		    
  H_MM2       = new TH1F("H_MM2",       "missing mass Squared, M^{2}_{miss}",                               MM2_nbins, MM2_xmin, MM2_xmax); 	    
  H_Em        = new TH1F("H_Emiss",     "missing energy",                                                   Em_nbins, Em_xmin, Em_xmax);   
  H_Em_nuc    = new TH1F("H_Em_nuc",    "nuclear missing energy",                                           Em_nbins, Em_xmin, Em_xmax);
  H_Em_recoil = new TH1F("H_Em_recoil", "recoil missing energy",                                            Em_nbins, Em_xmin, Em_xmax); 
  H_Pm        = new TH1F("H_Pm",        "missing momentum, P_{m}",                                          Pm_nbins, Pm_xmin, Pm_xmax); 
  H_Pmx_lab   = new TH1F("H_Pmx_Lab",   "P_{m, x} (Lab)",                                                   Pmx_lab_nbins, Pmx_lab_xmin, Pmx_lab_xmax);         
  H_Pmy_lab   = new TH1F("H_Pmy_Lab",   "P_{m, y} (Lab)",                                                   Pmy_lab_nbins, Pmy_lab_xmin, Pmy_lab_xmax);    
  H_Pmz_lab   = new TH1F("H_Pmz_Lab",   "P_{m, z} (Lab)",                                                   Pmz_lab_nbins, Pmz_lab_xmin, Pmz_lab_xmax);  
  H_Pmx_q     = new TH1F("H_Pmx_q",     "P_{m, xq} (w.r.t #vec{q}) ",                                       Pmx_q_nbins, Pmx_q_xmin, Pmx_q_xmax);   
  H_Pmy_q     = new TH1F("H_Pmy_q",     "P_{m, yq} (w.r.t #vec{q}) ",                                       Pmy_q_nbins, Pmy_q_xmin, Pmy_q_xmax); 
  H_Pmz_q     = new TH1F("H_Pmz_q",     "P_{m, zq} (along #vec{q}) ",                                       Pmz_q_nbins, Pmz_q_xmin, Pmz_q_xmax); 
  H_Tx        = new TH1F("H_Tx",        Form(" %s kinetic energy, T_{x} (detected)", det_had.Data()),       Tx_nbins, Tx_xmin, Tx_xmax);     
  H_Tr        = new TH1F("H_Tr",        "recoil system kinetic energy, T_{r} (recoil)",                     Tr_nbins, Tr_xmin, Tr_xmax);  
  H_thxq      = new TH1F("H_thxq",      Form(" %s in-plane angle, #theta_{xq}", det_had.Data()),            thxq_nbins, thxq_xmin, thxq_xmax);
  H_thrq      = new TH1F("H_thrq",      "recoil system in-plane angle, #theta_{rq}",                        thrq_nbins, thrq_xmin, thrq_xmax);
  H_phxq      = new TH1F("H_phxq",      Form(" %s out-of-plane angle, #phi_{xq}", det_had.Data()),          phxq_nbins, phxq_xmin, phxq_xmax);
  H_phrq      = new TH1F("H_phrq",      "recoil system out-of-plane angle, #phi_{rq}",                      phrq_nbins, phrq_xmin, phrq_xmax);

  // 2d kinematics
  H_kf_vs_the  = new TH2F("H_kf_vs_the", "e^{-} Momentum vs. #theta_{e}; #theta_{e} [deg]; k_{f} [GeV/c]", the_nbins, the_xmin, the_xmax, kf_nbins, kf_xmin, kf_xmax);      

  
  // selected kin. @ Forward Detector 
  H_W_FD        =  new TH1F("H_W_FD",        "FD: Invariant Mass, W",       W_nbins,   W_xmin,   W_xmax);
  H_W_FD_sec1   =  new TH1F("H_W_FD_sec1",   "FD sec1: Invariant Mass, W",  W_nbins,   W_xmin,   W_xmax);
  H_W_FD_sec2   =  new TH1F("H_W_FD_sec2",   "FD sec2: Invariant Mass, W",  W_nbins,   W_xmin,   W_xmax);
  H_W_FD_sec3   =  new TH1F("H_W_FD_sec3",   "FD sec3: Invariant Mass, W",  W_nbins,   W_xmin,   W_xmax);
  H_W_FD_sec4   =  new TH1F("H_W_FD_sec4",   "FD sec4: Invariant Mass, W",  W_nbins,   W_xmin,   W_xmax);
  H_W_FD_sec5   =  new TH1F("H_W_FD_sec5",   "FD sec5: Invariant Mass, W",  W_nbins,   W_xmin,   W_xmax);
  H_W_FD_sec6   =  new TH1F("H_W_FD_sec6",   "FD sec6: Invariant Mass, W",  W_nbins,   W_xmin,   W_xmax);


  
  // Add Histos to TList

  //--------------------
  // Acceptance Histos
  //--------------------
  accp_HList->Add( H_the_vs_phe );
  
  //--------------------
  // Particle ID Histos
  //--------------------
  pid_HList->Add( H_chi2pid_elec);
  pid_HList->Add( H_chi2pid_had);  
  pid_HList->Add( H_beta_elec);
  pid_HList->Add( H_beta_had );
  pid_HList->Add( H_beta_vs_kf );
  pid_HList->Add( H_beta_vs_pf );
  pid_HList->Add( H_ztar_diff);
  pid_HList->Add( H_zE);
  pid_HList->Add( H_pf_vs_thxq );
  
  //--------------------
  // Kinematics Histos
  //-------------------
  
  // electron kinematic histos
  kin_HList->Add( H_kf_vx);
  kin_HList->Add( H_kf_vy);
  kin_HList->Add( H_kf_vz);
  kin_HList->Add( H_kf_vt);
  kin_HList->Add( H_kf   );
  kin_HList->Add( H_kfx  );
  kin_HList->Add( H_kfy  );
  kin_HList->Add( H_kfz  );
  kin_HList->Add( H_q    );
  kin_HList->Add( H_qx   );
  kin_HList->Add( H_qy   );
  kin_HList->Add( H_qz   );
  kin_HList->Add( H_nu   );
  kin_HList->Add( H_Q2   ); 
  kin_HList->Add( H_xbj  );
  kin_HList->Add( H_the  );
  kin_HList->Add( H_phe  );
  kin_HList->Add( H_thq  );
  kin_HList->Add( H_phq  );
  kin_HList->Add( H_W    );
  kin_HList->Add( H_W2   );
  
  
  // hadron kinematic histos
  kin_HList->Add(  H_pf_vx    );
  kin_HList->Add(  H_pf_vy    );
  kin_HList->Add(  H_pf_vz    );
  kin_HList->Add(  H_pf_vt    );
  kin_HList->Add(  H_pf       );
  kin_HList->Add(  H_pfx      );
  kin_HList->Add(  H_pfy      );
  kin_HList->Add(  H_pfz      );
  kin_HList->Add(  H_thx      );
  kin_HList->Add(  H_MM       );
  kin_HList->Add(  H_MM2      );
  kin_HList->Add(  H_Em       );
  kin_HList->Add(  H_Em_nuc   );
  kin_HList->Add(  H_Em_recoil);
  kin_HList->Add(  H_Pm       );
  kin_HList->Add(  H_Pmx_lab  );
  kin_HList->Add(  H_Pmy_lab  );
  kin_HList->Add(  H_Pmz_lab  );
  kin_HList->Add(  H_Pmx_q    );
  kin_HList->Add(  H_Pmy_q    );
  kin_HList->Add(  H_Pmz_q    );
  kin_HList->Add(  H_Tx       );
  kin_HList->Add(  H_Tr       );
  kin_HList->Add(  H_thxq     );
  kin_HList->Add(  H_thrq     );
  kin_HList->Add(  H_phxq     );
  kin_HList->Add(  H_phrq     );
 
  // 2d kinematics

  kin_HList->Add( H_kf_vs_the  );
 

  // selected kin. @ Forward Detector 
  kin_HList_FD->Add( H_W_FD );
  kin_HList_FD->Add( H_W_FD_sec1 );
  kin_HList_FD->Add( H_W_FD_sec2 );
  kin_HList_FD->Add( H_W_FD_sec3 );
  kin_HList_FD->Add( H_W_FD_sec4 );
  kin_HList_FD->Add( H_W_FD_sec5 );
  kin_HList_FD->Add( H_W_FD_sec6 );

  
}

//_______________________________________________________________________________
void e4nu_analyzer::CreateTree()
{

  cout << "Start e4nu_analyzer::CreateTree() ... " << endl;

  // Define TTree and variables to write to it
  data_tree = new TTree("data_tree","");

  // Define variables to put in TTree here
  // (remember, we want all particle types, and be able to make
  // cuts on particle types, eaach sector, and either the FD or CND)
  
  // define data TTree branches
  data_tree->Branch("vx",     &br_vx,    "vx/D");
  data_tree->Branch("vy",     &br_vy,    "vx/D");
  data_tree->Branch("vz",     &br_vz,    "vx/D");
  data_tree->Branch("vt",     &br_vt,    "vt/D");  
  data_tree->Branch("px",     &br_px,    "px/D");
  data_tree->Branch("py",     &br_py,    "px/D");
  data_tree->Branch("pz",     &br_pz,    "px/D");
  data_tree->Branch("beta",   &br_beta,  "beta/D");
  data_tree->Branch("pid",    &br_pid,   "pid/I");
  data_tree->Branch("npart",  &br_npart, "npart/I");
  data_tree->Branch("chi2pid",&br_chi2pid,"chi2pid/D");

}

//_______________________________________________________________________________
void e4nu_analyzer::EventLoop()
{

  cout << "Start e4nu_analyzer::EventLoop() ... " << endl;

  // add many hipo files together
  clas12root::HipoChain chain;
  chain.Add(ifname);
  chain.SetReaderTags({0});
  
  auto config_c12 = chain.GetC12Reader();
  chain.db()->turnOffQADB();

  auto& c12 = chain.C12ref();
  auto& rcdbData= config_c12->rcdb()->current();//struct with all relevent rcdb values        

   // get beam energy from database (not working, need to ask)
  Double_t Eb = 5.98636; // GeV  | this command does NOT work (it gives 0) : rcdbData.beam_energy/1000;
  
  //Define event counter
  int evt_cnt=0;
  
 // Loop over all events
  while(chain.Next())
    {
      
      //cout << "===========" << endl;
      //cout << "event #: " << evt_cnt << endl;
      //cout << "===========" << endl;

      // get std::vector array of all final state particles per event
      auto particles = c12->getDetParticles();
      br_npart = particles.size();
      //cout << "# of particles: " << particles.size() << endl;
      for(auto& p : particles)
	{
	  // get all the info for every particle type here, and fill the tree outside
	  br_vx = p->par()->getVx();
	  br_vy = p->par()->getVy();
	  br_vz = p->par()->getVz();
	  br_vt = p->par()->getVt();
	  br_px = p->par()->getPx();
	  br_py = p->par()->getPy();
	  br_pz = p->par()->getPz();
	  br_p  = p->par()->getP();
	  br_beta = p->par()->getBeta();
	  br_pid = p->par()->getPid();	  
	  br_chi2pid = p->par()->getChi2Pid();

	  //cout << "pid: " << br_pid << endl;
	  //cout << " region: " << p->getRegion() << endl;
	  //cout << "hadron region: " << p->getRegion() << endl;
	  
	  //cout << " sector: " << p->getSector() << endl;
	  //cout << "hadron sector: " << p->getSector() << endl;
	
	}
      
      //Fill Tree Here !
      data_tree->Fill();
      
      // get std::vector array of final state particles for the reaction of interest
      // if there are multiple of the same particle, e.g. 3 protons, but only want 1 proton,
      // then require protons.size()>=1 later on.
      auto electrons = c12->getByID(11);
      auto protons   = c12->getByID(2212);
      auto neutrons  = c12->getByID(2112);
      auto pip       = c12->getByID(211);  // pion +
      auto pim       = c12->getByID(-211); // pion -
      auto pi0       = c12->getByID(111);  // pion 0


      // define boolean flag for reaction of interest
      Bool_t Aeep   = false;   // single proton knockout, (e,e'p)
      Bool_t Aeepip = false;  // pion+ electro-production (e,e'pi+) : ep -> e'pi+ (X) : virtual photon strikes quark in nucleon, producing a pi+  + (continuum of "missing" particles)
      Bool_t Aeepim = false;  // pion- electro-production (e,e'pi-) : ep -> e'pi- (X) : virtual photon strikes quark in neutron, producing a pi-  +  (continuum of "missing" particles)
      Bool_t Aeepi0 = false;  // pion- electro-production (e,e'pi0) : ep -> e'pi0 (X) : virtual photon strikes quark in neutron, producing a pi0  +  (continuum of "missing" particles)

      //======================================================
      //
      // select final state particles of interest to analyze
      //
      //======================================================
      
      // single e-proton knockout, A(e,e'p), require exactly 1 electron, at least one proton and at least 2 particles in final state 
      Aeep   = electrons.size()==1 && protons.size()>=1 && particles.size()>=2;
      
      // pion+ electro-production A(e,e'pi+) : eN -> e'pi+ (X) :  virtual photon strikes quark in nucleon, producing a pi+ + (continuum of "missing" particles)
      Aeepip = electrons.size()==1 && pip.size()==1 && particles.size()>=2;

      // pion- electro-production A(e,e'pi-) : eN -> e'pi+ (X) :  virtual photon strikes quark in nucleon, producing a pi- + (continuum of "missing" particles)
      Aeepim = electrons.size()==1 && pim.size()>=1 && particles.size()>=2;

      // pion- electro-production A(e,e'pi0) : eN -> e'pi0 (X) :  virtual photon strikes quark in nucleon, producing a pi0 + (continuum of "missing" particles)
      Aeepi0 = electrons.size()==1 && pi0.size()>=1 && particles.size()>=2;

      Bool_t reaction_type = false;

      if(det_had=="p"){

	reaction_type = Aeep;

	if(reaction_type){
	  // knocked-out protons vertex (i.e., interaction point location)
	  pf_vx = protons[0]->par()->getVx();
	  pf_vy = protons[0]->par()->getVy();
	  pf_vz = protons[0]->par()->getVz();
	  pf_vt = protons[0]->par()->getVt();
	  
	  // knocked-out protons 3-momentum
	  pf_x = protons[0]->par()->getPx(); 
	  pf_y = protons[0]->par()->getPy();
	  pf_z = protons[0]->par()->getPz(); 
	  pf   = protons[0]->par()->getP();
	  
	  // knocked-out protons beta
	  h_beta = protons[0]->par()->getBeta();
	  h_chi2pid = protons[0]->par()->getChi2Pid();
	  
	  // detected particle in-plane angle
	  th_x = protons[0]->getTheta()*TMath::RadToDeg();
	}
      }
      
      else if(det_had=="pi+"){
	
	reaction_type = Aeepip;

	if(reaction_type){
	  // electro-produced pions+ vertex (i.e., interaction point location)
	  pf_vx = pip[0]->par()->getVx();
	  pf_vy = pip[0]->par()->getVy();
	  pf_vz = pip[0]->par()->getVz();
	  pf_vt = pip[0]->par()->getVt();
	  
	  // electro-produced pions+ 3-momentum
	  pf_x = pip[0]->par()->getPx(); 
	  pf_y = pip[0]->par()->getPy();
	  pf_z = pip[0]->par()->getPz(); 
	  pf   = pip[0]->par()->getP();
	  
	  // knocked-out pions+ beta
	  h_beta = pip[0]->par()->getBeta();
	  h_chi2pid = pip[0]->par()->getChi2Pid();
	   
	  // detected particle in-plane angle
	  th_x = pip[0]->getTheta()*TMath::RadToDeg();
	}
      }
      
      else if(det_had=="pi-"){
	
	reaction_type = Aeepim;

	if(reaction_type){
	  // electro-produced pions- vertex (i.e., interaction point location)
	  pf_vx = pim[0]->par()->getVx();
	  pf_vy = pim[0]->par()->getVy();
	  pf_vz = pim[0]->par()->getVz();
	  pf_vt = pim[0]->par()->getVt();
	  
	  // electro-produced pions- 3-momentum
	  pf_x = pim[0]->par()->getPx(); 
	  pf_y = pim[0]->par()->getPy();
	  pf_z = pim[0]->par()->getPz(); 
	  pf   = pim[0]->par()->getP();
	  
	  // knocked-out pions- beta
	  h_beta = pim[0]->par()->getBeta();
	  h_chi2pid = pim[0]->par()->getChi2Pid();
	   
	  // detected particle in-plane angle
	  th_x = pim[0]->getTheta()*TMath::RadToDeg();
	}
      }

      else if(det_had=="pi0"){
	
	reaction_type = Aeepim;
	
	if(reaction_type){
	  // electro-produced pions- vertex (i.e., interaction point location)
	  pf_vx = pi0[0]->par()->getVx();
	  pf_vy = pi0[0]->par()->getVy();
	  pf_vz = pi0[0]->par()->getVz();
	  pf_vt = pi0[0]->par()->getVt();
	  
	  // electro-produced pions- 3-momentum
	  pf_x = pi0[0]->par()->getPx(); 
	  pf_y = pi0[0]->par()->getPy();
	  pf_z = pi0[0]->par()->getPz(); 
	  pf   = pi0[0]->par()->getP();
	  
	  // knocked-out pions- beta
	  h_beta = pi0[0]->par()->getBeta();
	  h_chi2pid = pi0[0]->par()->getChi2Pid();
	   
	  // detected particle in-plane angle
	  th_x = pi0[0]->getTheta()*TMath::RadToDeg();
	}
      }
      
      
      if ( reaction_type ){
	

	//NOTE: if there are multiple protons, the proton with highest momentum is the most
	// likely to have been directly hit by the virtual photon, so in order to get the maximum, one can do:
	//double leading_proton = *max_element(protons_momentum.begin(), protons_momentum.end());, where protons_momentum [] is an array
	// NOTE2:  I found that actually, for an array of protons, the array is organized from higher to lower momentum, so protons[0]->par()->getP() will be the leading
	// and protons[1]->par()->getP() will be the second, ...

	/*
	std::vector<double> pvec;
	for (int i=0; i<protons.size(); i++){
	  pvec.push_back( protons[i]->par()->getP() );
	}
	
	// get maximum momentum value
	double pvec_max =  *max_element( pvec.begin(), pvec.end() );
	double idx = max_element(pvec.begin(),pvec.end()) - pvec.begin();
	 
	std::cout << "The vector elements are : ";
	
	for(int i=0; i < pvec.size(); i++){
	  std::cout << pvec.at(i) << ' ';
	}
	
	cout << "pvec_max = " << pvec_max << endl;
	cout << "pvec_max idx = " << idx << endl;
	*/

	
	// --- get vertex/momentum components of scattered electron ---

	// scattered electrons vertex (i.e., interaction point location)
	kf_vx = electrons[0]->par()->getVx();
	kf_vy = electrons[0]->par()->getVy();
	kf_vz = electrons[0]->par()->getVz();
	kf_vt = electrons[0]->par()->getVt();
	
	// scattered electrons 3-momentum
	kf_x = electrons[0]->par()->getPx(); 
	kf_y = electrons[0]->par()->getPy();
	kf_z = electrons[0]->par()->getPz(); 
	kf   = electrons[0]->par()->getP();

	// scattered electron beta
	e_beta = electrons[0]->par()->getBeta();	
	e_chi2pid = electrons[0]->par()->getChi2Pid();
	 
	// set 4-momenta of beam, target, scattered electron and primary hadron detected 
	p4_beam.SetXYZM(0., 0., Eb, me);
	p4_target.SetXYZM(0.,0.,0., Mt);
	p4_electron.SetXYZM(kf_x, kf_y, kf_z, me);
	p4_hadron.SetXYZM(pf_x, pf_y, pf_z, Mh);		
	
	// 4-momentum transferred
	p4_q = p4_beam - p4_electron;

	// calculate electron kinematics
	Q2 = -p4_q.M2(); // 4-momentum squared (Q^2) or "virtuality"
	nu = p4_q.E();   // energy part ( nu ) of 4-vector p4_q
	xbj = Q2 / (2.*MP*nu); 
	q = p4_q.P();    // 3-momentum part |q| of 4-vector p4_q
	qx = p4_q.Px();
	qy = p4_q.Py();
	qz = p4_q.Pz();
	th_q = p4_q.Theta()*TMath::RadToDeg();
	ph_q = p4_q.Phi()*TMath::RadToDeg();
	th_e = electrons[0]->getTheta()*TMath::RadToDeg();
	ph_e = electrons[0]->getPhi()*TMath::RadToDeg();

	//invariant mass squared 
	W2 = (p4_q + p4_target).M2();
	W =  (p4_q + p4_target).M();
	//if(W2>0){
	//  W = sqrt(W2); // INVARIANT MASS
	//}

	
	// calculate hadron (and "missing" particles) kinematics
	// This we want to make a function and call for different hadrons
	// 4-momentum of undetected recoil system
	p4_recoil = p4_q + p4_target - p4_hadron;  

	// missing momentum components in the lab coordinate system
	Pmx_lab = p4_recoil.X();
	Pmy_lab = p4_recoil.Y();
	Pmz_lab = p4_recoil.Z();
	
	// Calculate angles of detected and recoil system wrt q-vector 
	// xq and bq are the 3-momentum vectors of "detected" and "recoil" expressed in
	// the coordinate system where q is the z-axis and the x-axis
	// lies in the scattering plane (defined by q and e') and points
	// in the direction of e', so the out-of-plane angle lies within
	// -90<phi_xq<90deg if the hadron is detected on the downstream/forward side of q.
	TRotation rot_to_q;
	rot_to_q.SetZAxis( p4_q.Vect(), p4_electron.Vect()).Invert();
	TVector3 xq = p4_hadron.Vect();
	TVector3 bq = p4_recoil.Vect();
	xq *= rot_to_q;  
	bq *= rot_to_q;

	th_xq     = xq.Theta()*TMath::RadToDeg();   //"theta_xq" (in-plane angle of detected particle relative to q-vector)
	ph_xq     = xq.Phi()*TMath::RadToDeg();     //"out-of-plane angle", "phi" 
	th_rq     = bq.Theta()*TMath::RadToDeg();
	ph_rq     = bq.Phi()*TMath::RadToDeg();

	// Missing momentum and components wrt q-vector
	// The definition of p_miss as the negative of the undetected recoil
	// momentum is the standard nuclear physics convention.
	TVector3 p_miss = -bq;
	Pm   = p_miss.Mag();  

	// The missing momentum components in the q coordinate system.
	Pmx_q = p_miss.X();
	Pmy_q = p_miss.Y();
	Pmz_q = p_miss.Z();


	// invariant mass of the recoil system a.k.a. "missing" mass 
	MM2 = p4_recoil.M2();  
	MM = p4_recoil.M();
	
	//if(MM2>0){
	//  MM = sqrt(MM2); // INVARIANT MASS
	//}

	// calculate kinetic energues of detected (x) and recoil (r) system
	Tx =  p4_hadron.E() - p4_hadron.M();
	Tr =  p4_recoil.E() - p4_recoil.M();
	
	// ------- missing energy calculations -------

	// does not apply for hydrogen elastics (definition for electro-breakup of nucleons)
	// --> mass diff. between free and bounf nucleons
	Em_nuc    = nu - Tx - Tr;
	
	// definition for electro-production reactions (virtual photon hits quark in nucleon leading to electro-production of other particles)
	Em        = nu + p4_target.M() - p4_hadron.E(); 
	Em_recoil = p4_recoil.E();

	//reaction vertex z difference
	ztar_diff = kf_vz - pf_vz;  

	//energy transferred to the hadron
	zE = p4_hadron.E() / nu;
	
	//cout << "MM = " << MM << endl;
	//double MM_v2 = sqrt(Em*Em - Pm*Pm);
	//cout << "MM_v2 = " << MM_v2 << endl;
	
	//--------------DEFINE CUTS--------------------
	
	//----PID Cuts----

	//chi2pid on hadrons
	if(chi2pid_cut_flag) {c_chi2pid = h_chi2pid >=c_chi2pid_min &&   h_chi2pid <=c_chi2pid_max;}
	else{c_chi2pid=1;}

	// z-reaction vertex difference between electron and detected hadron
	if(ztarDiff_cut_flag) {c_ztarDiff = ztar_diff>=c_ztarDiff_min && ztar_diff<=c_ztarDiff_max;}
	else{c_ztarDiff=1;}

	// z = E_had/nu (fraction of energy transferred to hadron) Cut
	if(zE_cut_flag) {c_zE = zE>=c_zE_min && zE<=c_zE_max;}
	else{c_zE=1;}
	
	c_pidCuts = c_chi2pid && c_ztarDiff && c_zE;
	
	//----Kinematic Cuts----
	//Q2
	if(Q2_cut_flag){c_Q2 = Q2>=c_Q2_min && Q2<=c_Q2_max;}
	else{c_Q2=1;}

	//Missing Energy, Em (assuming electro-production or H(e,e'p) )
	
	if(Em_cut_flag){
	  if(det_had=="p" && (target!="H") ){c_Em = Em_nuc>=c_Em_min && Em_nuc<=c_Em_max;} // breakup reactions A(e,e'p)
	  else {c_Em = Em>=c_Em_min && Em<=c_Em_max;} 
	} 
	else{c_Em=1;}
		
	//Invariant Mass, W
	if(W_cut_flag){c_W = W>=c_W_min && W<=c_W_max;}
	else{c_W=1;}

	//Missing Mass, MM
	if(MM_cut_flag){c_MM = MM>=c_MM_min && MM<=c_MM_max;}
	else{c_MM=1;}
	
	c_kinCuts = c_Q2 && c_Em && c_W && c_MM;
	
	// combine all cuts
	c_allCuts = c_pidCuts && c_kinCuts;

       	
	//Apply cuts
	if( c_allCuts ) {

	  // Fill Histograms
	  
	  //--------------------
	  // Acceptance Histos
	  //--------------------
	  H_the_vs_phe ->Fill(ph_e, th_e);
	  
	  //--------------------
	  // Particle ID Histos
	  //--------------------
	  H_chi2pid_elec->Fill(e_chi2pid);
	  H_chi2pid_had->Fill(h_chi2pid);
	  H_beta_elec->Fill(e_beta);
	  H_beta_had->Fill(h_beta);
	  H_beta_vs_kf ->Fill(kf, e_beta);
	  H_beta_vs_pf ->Fill(pf, h_beta);
	  H_ztar_diff ->Fill(ztar_diff);
	  H_zE -> Fill(zE);
	  H_pf_vs_thxq ->Fill(th_xq, pf);
	    
	  //--------------------
	  // Kinematics Histos
	  //-------------------
	  
	  // electron kinematics
	  H_kf_vx  ->Fill(kf_vx);
	  H_kf_vy  ->Fill(kf_vy);
	  H_kf_vz  ->Fill(kf_vz);
	  H_kf_vt  ->Fill(kf_vt); 
	  H_kf  ->Fill(kf);
	  H_kfx ->Fill(kf_x); 
	  H_kfy ->Fill(kf_y); 
	  H_kfz ->Fill(kf_z);
	  H_q   ->Fill(q);
	  H_qx  ->Fill(qx);
	  H_qy  ->Fill(qy);
	  H_qz  ->Fill(qz);  
	  H_nu  ->Fill(nu);  
	  H_Q2  ->Fill(Q2);
	  H_xbj ->Fill(xbj); 	
	  H_the ->Fill(th_e);
	  H_phe ->Fill(ph_e);		
	  H_thq ->Fill(th_q);
	  H_phq ->Fill(ph_q);	
	  H_W   ->Fill(W);   
	  H_W2  ->Fill(W2);  	
	  
	  
	  // hadron kinematics
	  H_pf_vx  ->Fill(pf_vx);
	  H_pf_vy  ->Fill(pf_vy);
	  H_pf_vz  ->Fill(pf_vz);
	  H_pf_vt  ->Fill(pf_vt); 
	  H_pf   ->Fill(pf);
	  H_pfx  ->Fill(pf_x);
	  H_pfy  ->Fill(pf_y);
	  H_pfz  ->Fill(pf_z);
	  H_thx  ->Fill(th_x);
	  H_MM   ->Fill(MM);  
	  H_MM2  ->Fill(MM2); 
	  H_Em   ->Fill(Em);	  
	  H_Em_nuc    ->Fill(Em_nuc);
	  H_Em_recoil ->Fill(Em_recoil); 
	  H_Pm        ->Fill(Pm);	  
	  H_Pmx_lab   ->Fill(Pmx_lab);
	  H_Pmy_lab   ->Fill(Pmy_lab);
	  H_Pmz_lab   ->Fill(Pmz_lab);
	  H_Pmx_q     ->Fill(Pmx_q);  
	  H_Pmy_q     ->Fill(Pmy_q);  
	  H_Pmz_q     ->Fill(Pmz_q);  
	  H_Tx        ->Fill(Tx);	  
	  H_Tr        ->Fill(Tr);	  
	  H_thxq      ->Fill(th_xq);	  
	  H_thrq      ->Fill(th_rq);	  
	  H_phxq      ->Fill(ph_xq);	  
	  H_phrq      ->Fill(ph_rq);
	  
	  
	  // 2d kinematics
	  H_kf_vs_the ->Fill(th_e, kf);
	  
	  // Fill certain kin. variables per region (either detected in Forward or Central Detector, FD - 2000, CD - 4000)
	  
	  //cout << "electrons[0]->getRegion() --> " << electrons[0]->getRegion() << endl;
	  //Forward Detector
	  if(electrons[0]->getRegion()==2000){
	    //cout <<  "electrons[0]->getSector() --> " <<  electrons[0]->getSector() << endl;
	    H_W_FD->Fill(W);
	    if( electrons[0]->getSector()==1 ) {H_W_FD_sec1->Fill(W);}
	    if( electrons[0]->getSector()==2 ) {H_W_FD_sec2->Fill(W);}
	    if( electrons[0]->getSector()==3 ) {H_W_FD_sec3->Fill(W);}
	    if( electrons[0]->getSector()==4 ) {H_W_FD_sec4->Fill(W);}
	    if( electrons[0]->getSector()==5 ) {H_W_FD_sec5->Fill(W);}
	    if( electrons[0]->getSector()==6 ) {H_W_FD_sec6->Fill(W);}
	    
	    
	  }
	  
	} // end analysis cuts
	
      } // end final state particle requirement
      
      
      
     
      /* NOTE: particles detected are separated into central detector (CD) and forward detector (FD) ( < 35 deg in-plane)
	 From Justin Esteeves, currently the CD has much worse resolution and so one should look at events reconstructed from
	 these detectors separately in order to compare how well reconstruction is done.
      */

	 

      if (evt_cnt % 1000 == 0){  
	cout << "Events Analyzed: " << std::setprecision(2) << double(evt_cnt) << std::flush << "\r";
      }

      // increment event number
      evt_cnt++;

      
    }
  
}

//_______________________________________________________________________________
void e4nu_analyzer::WriteHist()
{
  cout << "Start e4nu_analyzer::WriteHist() ... " << endl;

  //Create Output ROOTfile
  outROOT = new TFile(ofname.Data(), "RECREATE");

  //Make directories to store histograms based on category
  outROOT->mkdir("accp_plots");
  outROOT->mkdir("pid_plots");
  outROOT->mkdir("kin_plots");
  outROOT->mkdir("kin_plots_FD");
      
  //Write Kinematics histos to kin_plots directory
  outROOT->cd("accp_plots");
  accp_HList->Write();

  outROOT->cd("pid_plots");
  pid_HList->Write();
  
  outROOT->cd("kin_plots");
  kin_HList->Write();

  outROOT->cd("kin_plots_FD");
  kin_HList_FD->Write();
  
  outROOT->Close();
}

//_______________________________________________________________________________
void e4nu_analyzer::run_data_analysis()
{
  
 
  CreateHist();
  SetCuts();
  CreateTree();
  EventLoop();
  WriteHist();
  
}
