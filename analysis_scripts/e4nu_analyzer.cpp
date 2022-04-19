/*=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=

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
// include this analyzer header file
#include "e4nu_analyzer.h"

using namespace clas12;
using namespace std;

//_______________________________________________________________________________
e4nu_analyzer::e4nu_analyzer(TString inHIPO_fname="", TString outROOT_fname="", TString tgt="" )
: ifname(inHIPO_fname), ofname(outROOT_fname), target(tgt)
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
  
  // initialize histogram pointers
  H_the  =  NULL; 
  H_kf   =  NULL; 
  H_W    =  NULL; 
  H_W2   =  NULL; 
  H_Q2   =  NULL; 
  H_xbj  =  NULL; 
  H_nu   =  NULL; 
  H_q    =  NULL; 
  H_qx 	 =  NULL; 
  H_qy 	 =  NULL; 
  H_qz 	 =  NULL; 
  H_thq	 =  NULL; 
  H_phq	 =  NULL;  
  H_MM   =  NULL; 
  H_MM2  =  NULL; 





}

//_______________________________________________________________________________
e4nu_analyzer::~e4nu_analyzer()
{

    cout << "Start e4nu_analyzer::~e4nu_analyzer ... " << endl;

    // delete TFile Pointers
    delete outROOT;   outROOT = NULL;
    
    
    // delete TTree Pointers
    delete data_tree; data_tree = NULL;
    
    // delete histogram pointers
    delete H_the ;  H_the  =  NULL; 
    delete H_kf  ;  H_kf   =  NULL; 
    delete H_W   ;  H_W    =  NULL; 
    delete H_W2  ;  H_W2   =  NULL; 
    delete H_Q2  ;  H_Q2   =  NULL; 
    delete H_xbj ;  H_xbj  =  NULL; 
    delete H_nu  ;  H_nu   =  NULL; 
    delete H_q   ;  H_q    =  NULL; 
    delete H_qx  ;  H_qx   =  NULL; 
    delete H_qy  ;  H_qy   =  NULL; 
    delete H_qz  ;  H_qz   =  NULL; 
    delete H_thq ;  H_thq  =  NULL; 
    delete H_phq ;  H_phq  =  NULL; 
    delete H_MM  ;  H_MM   =  NULL; 
    delete H_MM2 ;  H_MM2  =  NULL; 
    

    
    

}

//_______________________________________________________________________________
void e4nu_analyzer::SetParticleMass()
{
  
  cout << "Start e4nu_analyzer::SetParticleMass() ... " << endl;

  // load particle database from PDG to get particle mass
  auto db=TDatabasePDG::Instance();
  
  MP = db->GetParticle(2212)->Mass();  
  me = db->GetParticle(11)->Mass();
  
}

//_______________________________________________________________________________
void e4nu_analyzer::SetHistBins()
{
    cout << "Start e4nu_analyzer::SetHistBins() ... " << endl;


}

//_______________________________________________________________________________
void e4nu_analyzer::CreateHist()
{
  
  cout << "Start e4nu_analyzer::CreateHist() ... " << endl;

  SetHistBins();

  H_the     = new TH1F("H_the", "Electron Scattering Angle, #theta_{e}", 200, 0.5, 180);	    
  H_kf      = new TH1F("H_kf",  "Final e^{-} Momentum", 100, 0.5, 6);			    
  H_W       = new TH1F("H_W",   "Invariant Mass, W", 100, 0.1, 5); 				    
  H_W2      = new TH1F("H_W2",  "Invariant Mass, W^{2}", 100, -5, 5); 			    
  H_Q2      = new TH1F("H_Q2",  "4-Momentum Transfer, Q^{2}", 100, 0.1, 5); 		    
  H_xbj     = new TH1F("H_xbj", "x-Bjorken", 100, 0.2, 2);  				    
  H_nu      = new TH1F("H_nu",  "Energy Transfer, #nu", 100, 0.2, 6); 			    
  H_q       = new TH1F("H_q",   "3-Momentum Transfer, |#vec{q}|", 100, 0.2, 6);		    
  H_qx      = new TH1F("H_qx",  "|#vec{q}_{x}|", 100, 0.2, 6);				    
  H_qy      = new TH1F("H_qy",  "|#vec{q}_{y}|", 100, 0.2, 6);				    
  H_qz      = new TH1F("H_qz",  "|#vec{q}_{z}|",                              100,  0.2, 6     );				    
  H_thq     = new TH1F("H_thq", "In-Plane Angle w.r.t +z(lab), #theta_{q}",   100,  0,   180   );     
  H_phq     = new TH1F("H_phq", "Out-of-Plane Angle w.r.t +z(lab), #phi_{q}", 100, -180, 180   );       											    
  H_MM      = new TH1F("H_MM",  "Missing Mass, M_{miss}",                     100, -5,   5     );        		    
  H_MM2     = new TH1F("H_MM2", "Missing Mass Squared, M^{2}_{miss}",         100, -5,   5     ); 	    



  
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
  data_tree->Branch("px",&px,"px/D");
  data_tree->Branch("py",&py,"px/D");
  data_tree->Branch("pz",&pz,"px/D");
  data_tree->Branch("pid",&pid,"pid/I");
  data_tree->Branch("chi2pid",&chi2pid,"chi2pid/D");

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

      
      // TESTING / DEBUGGING
      //cout << "beamE = " << beam_energy << endl;

      // get std::vector array of all final state particles per event
      auto particles = c12->getDetParticles();

      for(auto& p : particles)
	{
	  // get all the info for every particle type here, and fill the tree outside

	  px = p->par()->getPx();
	  py = p->par()->getPy();
	  pz = p->par()->getPz();
	  pid = p->par()->getPid();
	  chi2pid = p->par()->getChi2Pid();
	  
	}
      
      //Fill Tree Here !
      data_tree->Fill();
      
      // get std::vector array of final state particles for the reaction of interest
      // if there are multiple of the same particle, e.g. 3 protons, but only want 1 proton,
      // then require protons.size()==1 later on.
      auto electrons = c12->getByID(11);
      auto protons   = c12->getByID(2212);


      //select final state particle of interest to analyze
      if (electrons.size()==1 && protons.size()==1){


	// --- get momentum components of final state particles ---

	// scattered electrons
	double px_e = electrons[0]->par()->getPx(); 
	double py_e = electrons[0]->par()->getPy();
	double pz_e = electrons[0]->par()->getPz(); 
	double p_e = electrons[0]->par()->getP();
	// knocked-out protons
	double px_p = protons[0]->par()->getPx(); 
	double py_p = protons[0]->par()->getPy();
	double pz_p = protons[0]->par()->getPz(); 

	// set 4-momenta of particle of interest
	p4_beam.SetXYZM(0., 0., Eb, me);
	p4_target.SetXYZM(0.,0.,0., MP);
	p4_electron.SetXYZM(px_e, py_e, pz_e, me);
	p4_proton.SetXYZM(px_p, py_p, pz_p, MP);		
	
	// 4-momentum trasnferred
	p4_q = p4_beam - p4_electron;

	// define electron kinematics
	double Q2 = -p4_q.M2(); // 4-momentum squared (Q^2) or "virtuality"
	double q = p4_q.P();    // momentum part |q| of 4-vector p4_q
	double th_q = p4_q.Theta()*TMath::RadToDeg();
	double ph_q = p4_q.Phi()*TMath::RadToDeg();
	double nu = p4_q.E();   // energy part ( nu ) of 4-vector p4_q
	double xbj = Q2 / (2.*MP*nu); 
	double th_e = p4_beam.Angle( p4_electron.Vect() ) *TMath::RadToDeg();
	
	// define recoil (usually undetected) system kinematics
	p4_recoil = p4_q + p4_target - p4_proton;  

	double W = -1000.;
	double MM = -1000.;

	double MM2 = p4_recoil.M2();  //recoil ("MISSING") mass squared of the system
	if(MM2>0){
	  MM = sqrt(MM2); // INVARIANT MASS
	}
	
	double W2 = (p4_q + p4_target).M2(); //invariant mass squared 
	if(W2>0){
	  W = sqrt(W2); // INVARIANT MASS
	}
	
	
	// Fill Histograms 
	H_the ->Fill(th_e);
	H_kf  ->Fill(p_e); 
	H_W   ->Fill(W);   
	H_W2  ->Fill(W2);  
	H_Q2  ->Fill(Q2);  
	H_xbj ->Fill(xbj); 
	H_nu  ->Fill(nu);  
	H_q   ->Fill(q);      
	H_thq ->Fill(th_q);
	H_phq ->Fill(ph_q);
	     		   
	H_MM  ->Fill(MM);  
	H_MM2 ->Fill(MM2); 

	
      }

      
      
     
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
  H_the ->Write();
  H_kf  ->Write(); 
  H_W   ->Write();   
  H_W2  ->Write();  
  H_Q2  ->Write();  
  H_xbj ->Write(); 
  H_nu  ->Write();  
  H_q   ->Write();   
  H_thq ->Write();
  H_phq ->Write();
  
  H_MM  ->Write();    
  H_MM2 ->Write();

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
