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
  : ifname(inHIPO_fname), ofname(outROOT_fname), target(tgt), detected_hadron(det_h)
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
  
  // --- initialize histogram pointers ----

  // electron kinematics
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
  H_the_vs_phe = NULL;
  H_kf_vs_the  = NULL;

}

//_______________________________________________________________________________
e4nu_analyzer::~e4nu_analyzer()
{

    cout << "Start e4nu_analyzer::~e4nu_analyzer ... " << endl;

    // delete TFile Pointers
    delete outROOT;   outROOT = NULL;
    
    
    // delete TTree Pointers
    delete data_tree; data_tree = NULL;
    
    // --- delete histogram pointers ---

    // electron
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
    delete H_the_vs_phe; H_the_vs_phe = NULL;
    delete H_kf_vs_the;  H_kf_vs_the  = NULL;
    
}

//_______________________________________________________________________________
void e4nu_analyzer::SetParticleMass()
{
  
  cout << "Start e4nu_analyzer::SetParticleMass() ... " << endl;

  // load particle database from PDG to get particle mass
  auto db=TDatabasePDG::Instance();
  
  MP = db->GetParticle(2212)->Mass();  
  me = db->GetParticle(11)->Mass();

  // set target mass
  if(target=="LH2"){ Mt = MP; }

  // set detected primary hadron mass A(e,e'p)X
  if(detected_hadron=="proton") { Mh = MP; }
  
}

//_______________________________________________________________________________
void e4nu_analyzer::SetHistBins()
{
  cout << "Start e4nu_analyzer::SetHistBins() ... " << endl;
  
  //---------------------------------
  // Kinematics Histograms Binning
  //---------------------------------

  TString input_HBinFileName = "histo_bins.inp";
  // detected electron
  
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
void e4nu_analyzer::CreateHist()
{
  
  cout << "Start e4nu_analyzer::CreateHist() ... " << endl;

  SetHistBins();

  // electron
  H_kf      = new TH1F("H_kf",  "Final e^{-} Momentum",                       kf_nbins,  kf_xmin,  kf_xmax );
  H_kfx     = new TH1F("H_kfx", "Final e^{-} Momentum (x-comp)",              kfx_nbins, kfx_xmin, kfx_xmax );
  H_kfy     = new TH1F("H_kfy", "Final e^{-} Momentum (y-comp)",              kfy_nbins, kfy_xmin, kfy_xmax );			    
  H_kfz     = new TH1F("H_kfz", "Final e^{-} Momentum (z-comp)",              kfz_nbins, kfz_xmin, kfz_xmax );
  H_q       = new TH1F("H_q",   "3-Momentum Transfer, |#vec{q}|",             q_nbins,   q_xmin,   q_xmax);		    
  H_qx      = new TH1F("H_qx",  "|#vec{q}_{x}|",                              qx_nbins,  qx_xmin,  qx_xmax);				    
  H_qy      = new TH1F("H_qy",  "|#vec{q}_{y}|",                              qy_nbins,  qy_xmin,  qy_xmax);				    
  H_qz      = new TH1F("H_qz",  "|#vec{q}_{z}|",                              qz_nbins,  qz_xmin,  qz_xmax);
  H_nu      = new TH1F("H_nu",  "Energy Transfer, #nu",                       nu_nbins,  nu_xmin,  nu_xmax);
  H_Q2      = new TH1F("H_Q2",  "4-Momentum Transfer, Q^{2}",                 Q2_nbins,  Q2_xmin,  Q2_xmax); 		    
  H_xbj     = new TH1F("H_xbj", "x-Bjorken",                                  xbj_nbins, xbj_xmin, xbj_xmax);
  H_the     = new TH1F("H_the", "Electron In-Plane Scattering Angle, #theta_{e}",      the_nbins, the_xmin, the_xmax);
  H_phe     = new TH1F("H_phe", "Electron Out-Of-Plane Scattering Angle, #phi_{e}",      phe_nbins, phe_xmin, phe_xmax);
  H_thq     = new TH1F("H_thq", "In-Plane Angle w.r.t +z(lab), #theta_{q}",   thq_nbins, thq_xmin, thq_xmax);     
  H_phq     = new TH1F("H_phq", "Out-of-Plane Angle w.r.t +z(lab), #phi_{q}", phq_nbins, phq_xmin, phq_xmax);
  H_W       = new TH1F("H_W",   "Invariant Mass, W",                          W_nbins,   W_xmin,   W_xmax); 				    
  H_W2      = new TH1F("H_W2",  "Invariant Mass, W^{2}",                      W2_nbins,  W2_xmin,  W2_xmax); 			     		 				    

  // hadron
  H_pf      = new TH1F("H_pf",    "Final Hadron Momentum (detected), p_{f}",              pf_nbins, pf_xmin, pf_xmax  );
  H_pfx      = new TH1F("H_pfx",  "Final Hadron Momentum, X-comp. p_{fx}",                    pfx_nbins, pfx_xmin, pfx_xmax );
  H_pfy      = new TH1F("H_pfy",  "Final Hadron Momentum, Y-comp. p_{fy}",                    pfy_nbins, pfy_xmin, pfy_xmax);
  H_pfz      = new TH1F("H_pfz",  "Final Hadron Momentum, Z-comp p_{fz}",                     pfz_nbins, pfz_xmin, pfz_xmax  );
  H_thx      = new TH1F("H_thx",  "Final Hadron Scatteting Angle (detected), #theta_{x}",              thx_nbins, thx_xmin, thx_xmax);
  H_MM      = new TH1F("H_MM",  "Missing Mass, M_{miss}",                       MM_nbins, MM_xmin, MM_xmax);        		    
  H_MM2     = new TH1F("H_MM2", "Missing Mass Squared, M^{2}_{miss}",          MM2_nbins, MM2_xmin, MM2_xmax); 	    
  H_Em      = new TH1F("H_Emiss","Missing Energy",                            Em_nbins, Em_xmin, Em_xmax);   
  H_Em_nuc    = new TH1F("H_Em_nuc","Nuclear Missing Energy",                    Em_nbins, Em_xmin, Em_xmax);
  H_Em_recoil  = new TH1F("H_Em_recoil","Recoil Missing Energy",                Em_nbins, Em_xmin, Em_xmax); 
  H_Pm      = new TH1F("H_Pm","Missing Momentum, P_{miss}",                    Pm_nbins, Pm_xmin, Pm_xmax); 
  H_Pmx_lab = new TH1F("H_Pmx_Lab","P_{miss, x} (Lab)",                    Pmx_lab_nbins, Pmx_lab_xmin, Pmx_lab_xmax);         
  H_Pmy_lab = new TH1F("H_Pmy_Lab","P_{miss, y} (Lab)",                    Pmy_lab_nbins, Pmy_lab_xmin, Pmy_lab_xmax);    
  H_Pmz_lab = new TH1F("H_Pmz_Lab","P_{miss, z} (Lab)",                    Pmz_lab_nbins, Pmz_lab_xmin, Pmz_lab_xmax);  
  H_Pmx_q   = new TH1F("H_Pmx_q","P_{miss, xq} (w.r.t #vec{q}) ",                    Pmx_q_nbins, Pmx_q_xmin, Pmx_q_xmax);   
  H_Pmy_q   = new TH1F("H_Pmy_q","P_{miss, yq} (w.r.t #vec{q}) ",                    Pmy_q_nbins, Pmy_q_xmin, Pmy_q_xmax); 
  H_Pmz_q   = new TH1F("H_Pmz_q","P_{miss, zq} (along #vec{q}) ",                    Pmz_q_nbins, Pmz_q_xmin, Pmz_q_xmax); 
  H_Tx      = new TH1F("H_Tx", "Kinetic Energy, T_{x} (detected)",                    Tx_nbins, Tx_xmin, Tx_xmax);     
  H_Tr      = new TH1F("H_Tr", "Kinetic Energy, T_{r} (recoil)",                      Tr_nbins, Tr_xmin, Tr_xmax);  
  H_thxq    = new TH1F("H_thxq", "In-Plane Angle, #theta_{xq}",                    thxq_nbins, thxq_xmin, thxq_xmax);
  H_thrq    = new TH1F("H_thrq", "In-Plane Angle, #theta_{rq}",                    thrq_nbins, thrq_xmin, thrq_xmax);
  H_phxq    = new TH1F("H_phxq", "Out-of-Plane Angle, #phi_{xq}",                    phxq_nbins, phxq_xmin, phxq_xmax);
  H_phrq    = new TH1F("H_phrq", "Out-of-Plane Angle, #phi_{rq}",                    phrq_nbins, phrq_xmin, phrq_xmax);

  // 2d kinematics
  H_the_vs_phe = new TH2F("H_the_vs_phe", "e^{-} #theta_{e} vs. # phi_{e}; #theta_{e} [deg]; #phi_{e} [deg]", phe_nbins, phe_xmin, phe_xmax, the_nbins, the_xmin, the_xmax);      
  H_kf_vs_the = new TH2F("H_kf_vs_the", "e^{-} Momentum vs. #theta_{e}; k_{f}; #theta_{e} [deg];", the_nbins, the_xmin, the_xmax, kf_nbins, kf_xmin, kf_xmax);      
  
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
  data_tree->Branch("px",     &br_px,    "px/D");
  data_tree->Branch("py",     &br_py,    "px/D");
  data_tree->Branch("pz",     &br_pz,    "px/D");
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
      cout << "===========" << endl;
      cout << "event #: " << evt_cnt << endl;
      cout << "===========" << endl;

      // get std::vector array of all final state particles per event
      auto particles = c12->getDetParticles();
      br_npart = particles.size();
      cout << "# of particles: " << particles.size() << endl;
      for(auto& p : particles)
	{
	  // get all the info for every particle type here, and fill the tree outside
	  
	  br_px = p->par()->getPx();
	  br_py = p->par()->getPy();
	  br_pz = p->par()->getPz();
	  br_p  = p->par()->getP();
	  br_pid = p->par()->getPid();	  
	  br_chi2pid = p->par()->getChi2Pid();

	}
      
      //Fill Tree Here !
      data_tree->Fill();
      
      // get std::vector array of final state particles for the reaction of interest
      // if there are multiple of the same particle, e.g. 3 protons, but only want 1 proton,
      // then require protons.size()==1 later on.
      auto electrons = c12->getByID(11);
      auto protons   = c12->getByID(2212);


      //======================================================
      //
      // select final state particle of interest to analyze
      //
      //======================================================
      
      // (e,e'p)
      if ( electrons.size()==1 && protons.size()==1 && particles.size()>=2 ){
	
	
	// --- get momentum components of final state particles ---

	// scattered electrons 3-momentum
	kf_x = electrons[0]->par()->getPx(); 
	kf_y = electrons[0]->par()->getPy();
	kf_z = electrons[0]->par()->getPz(); 
	kf   = electrons[0]->par()->getP();
	
	// knocked-out protons 3-momentum
	pf_x = protons[0]->par()->getPx(); 
	pf_y = protons[0]->par()->getPy();
	pf_z = protons[0]->par()->getPz(); 
	pf   = protons[0]->par()->getP();

	
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

	// detected particle in-plane angle
	th_x = protons[0]->getTheta()*TMath::RadToDeg();

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

	// calculate kinetic energues of detected and recoil system
	Tx =  p4_hadron.E() - p4_hadron.M();
	Tr =  p4_recoil.E() - p4_recoil.M();
	
	// missing energy calculations
	Em_nuc    = nu - Tx - Tr;
	Em        = nu + p4_target.M() - p4_hadron.E();
	Em_recoil = p4_recoil.E();
	
	
	
	// Fill Histograms 

	// electron kinematics
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
	H_the_vs_phe ->Fill(th_e, ph_e);
	H_kf_vs_the ->Fill(kf, th_e);
	  
      } // end final state particle requirement

      
      
     
      /* NOTE: particles detected are separated into central detector (CD) and forward detector (FD) ( < 35 deg in-plane)
	 From Justin Esteeves, currently the CD has much worse resolution and so one should look at events reconstructed from
	 these detectors separately in order to compare how well reconstruction is done.
      */
      // -- example from Justin ---
      /*
      if(p->getRegion()==CD)
	{
	  p_chi2_cd->Fill(p->par()->getChi2Pid());
	  p_chi2_cd_v = p->par()->getChi2Pid();
	}
      if(p->getRegion()==FD)
	{
	  p_chi2_fd->Fill(p->par()->getChi2Pid());
	  p_chi2_fd_v = p->par()->getChi2Pid();
	}
      */
	 
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
      
  outROOT->cd();
  
  // write leaf branches 
  data_tree->Write();
  
  // write histograms objects

  // electron kinematics
  H_kf  ->Write();
  H_kfx ->Write(); 
  H_kfy ->Write(); 
  H_kfz ->Write(); 
  H_q   ->Write();
  H_qx  ->Write();
  H_qy  ->Write();
  H_qz  ->Write();
  H_nu  ->Write();  
  H_Q2  ->Write();  
  H_xbj ->Write();
  H_the ->Write();
  H_thq ->Write();
  H_phq ->Write();
  H_W   ->Write();   
  H_W2  ->Write();

  // hadron (detected and "missing") kinematics
  H_pf  ->Write();
  H_pfx  ->Write();
  H_pfy  ->Write();
  H_pfz  ->Write();
  H_thx  ->Write();
  H_MM  ->Write();    
  H_MM2 ->Write();
  H_Em  ->Write();
  H_Em_nuc ->Write();
  H_Em_recoil ->Write();
  H_Pm      ->Write(); 
  H_Pmx_lab ->Write(); 
  H_Pmy_lab ->Write(); 
  H_Pmz_lab ->Write(); 
  H_Pmx_q   ->Write(); 
  H_Pmy_q   ->Write(); 
  H_Pmz_q   ->Write(); 
  H_Tx      ->Write(); 
  H_Tr      ->Write(); 
  H_thxq    ->Write(); 
  H_thrq    ->Write(); 
  H_phxq    ->Write(); 
  H_phrq    ->Write(); 
  
  outROOT->Close();
}

//_______________________________________________________________________________
void e4nu_analyzer::run_data_analysis()
{
  
 
  CreateHist();
  CreateTree();
  EventLoop();
  WriteHist();
  
}
