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
  
  // --- initialize histogram pointers ----

  // electron
  H_kf   =  NULL;
  H_kf_v2   =  NULL; 
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
  H_the_v2  =  NULL;
  H_thq	 =  NULL; 
  H_phq	 =  NULL;  
  H_W    =  NULL; 
  H_W2   =  NULL; 
 
  // hadron
  H_pf   =  NULL;
  H_pf_v2 = NULL;
  H_pfx  =  NULL;
  H_pfy  =  NULL;
  H_pfz  =  NULL;
  H_thh  =  NULL;
  H_thh_v2  =  NULL;
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
    
    // --- delete histogram pointers ---

    // electron
    delete H_kf  ;  H_kf   =  NULL;
    delete H_kf_v2  ;  H_kf_v2   =  NULL; 
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
    delete H_the_v2 ;  H_the_v2  =  NULL;
    delete H_thq ;  H_thq  =  NULL; 
    delete H_phq ;  H_phq  =  NULL; 
    delete H_W   ;  H_W    =  NULL; 
    delete H_W2  ;  H_W2   =  NULL;     

    // hadron
    delete H_pf  ;  H_pf   =  NULL;
    delete H_pf_v2  ;  H_pf_v2   =  NULL;
    delete H_pfx  ;  H_pfx   =  NULL;
    delete H_pfy  ;  H_pfy   =  NULL;
    delete H_pfz  ;  H_pfz   =  NULL;
    delete H_thh  ;  H_thh   = NULL;
    delete H_thh_v2  ;  H_thh_v2   = NULL;
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


  // electron
  H_kf      = new TH1F("H_kf",  "Final e^{-} Momentum",                       100,  0.5, 6  );
  H_kf_v2      = new TH1F("H_kf_v2",  "Final e^{-} Momentum",                 100,  0.5, 6  );
  H_kfx     = new TH1F("H_kfx","Final e^{-} Momentum (x-comp)",               100, -50,  50 );
  H_kfy     = new TH1F("H_kfy","Final e^{-} Momentum (y-comp)",               100, -50,  50 );			    
  H_kfz     = new TH1F("H_kfz","Final e^{-} Momentum (z-comp)",               100, -50,  50 );
  H_q       = new TH1F("H_q",   "3-Momentum Transfer, |#vec{q}|",             100,  0.2, 6  );		    
  H_qx      = new TH1F("H_qx",  "|#vec{q}_{x}|",                              100,  0.2, 6  );				    
  H_qy      = new TH1F("H_qy",  "|#vec{q}_{y}|",                              100,  0.2, 6  );				    
  H_qz      = new TH1F("H_qz",  "|#vec{q}_{z}|",                              100,  0.2, 6  );
  H_nu      = new TH1F("H_nu",  "Energy Transfer, #nu",                       100,  0.2, 6  );
  H_Q2      = new TH1F("H_Q2",  "4-Momentum Transfer, Q^{2}",                 100,  0.1, 5  ); 		    
  H_xbj     = new TH1F("H_xbj", "x-Bjorken",                                  100,  0.2, 2  );
  H_the     = new TH1F("H_the", "Electron Scattering Angle, #theta_{e}",      200,  0.5, 180);
  H_the_v2  = new TH1F("H_the_v2", "Electron Scattering Angle, #theta_{e}",   200,  0.5, 180);
  H_thq     = new TH1F("H_thq", "In-Plane Angle w.r.t +z(lab), #theta_{q}",   100,  0,   180);     
  H_phq     = new TH1F("H_phq", "Out-of-Plane Angle w.r.t +z(lab), #phi_{q}", 100, -180, 180);
  H_W       = new TH1F("H_W",   "Invariant Mass, W",                          100,  0.1, 5  ); 				    
  H_W2      = new TH1F("H_W2",  "Invariant Mass, W^{2}",                      100, -5,   5  ); 			     		 				    

  // hadron
  H_pf      = new TH1F("H_pf",  "Final Hadron Momentum, p_{f}",               100, 0,    5  );
  H_pf_v2      = new TH1F("H_pf_v2",  "Final Hadron Momentum, p_{f}",               100, 0,    5  );
  H_pfx      = new TH1F("H_pfx",  "Final Hadron Momentum, p_{fx}",               100, -5,    5  );
  H_pfy      = new TH1F("H_pfy",  "Final Hadron Momentum, p_{fy}",               100, -5,    5  );
  H_pfz      = new TH1F("H_pfz",  "Final Hadron Momentum, p_{fz}",               100, -5,    5  );
  H_thh      = new TH1F("H_thh",  "Final Hadron Angle, #theta_{h}",              100,  0,    180);
  H_thh_v2      = new TH1F("H_thh_v2",  "Final Hadron Angle, #theta_{h}",              100,  0,    180);
  H_MM      = new TH1F("H_MM",  "Missing Mass, M_{miss}",                     100, -5,   5  );        		    
  H_MM2     = new TH1F("H_MM2", "Missing Mass Squared, M^{2}_{miss}",         100, -5,   5  ); 	    



  
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
  data_tree->Branch("px",     &px,    "px/D");
  data_tree->Branch("py",     &py,    "px/D");
  data_tree->Branch("pz",     &pz,    "px/D");
  data_tree->Branch("pid",    &pid,   "pid/I");
  data_tree->Branch("npart",  &npart, "npart/I");
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

      // get std::vector array of all final state particles per event
      auto particles = c12->getDetParticles();
      npart = particles.size();
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


      //select final state particle of interest to analyze

      // (e,e'p) 
      if ( electrons.size()==1 && protons.size()==1 && particles.size()==2 ){
	
	
	// --- get momentum components of final state particles ---

	// scattered electrons 3-momentum
	kf_x = electrons[0]->par()->getPx(); 
	kf_y = electrons[0]->par()->getPy();
	kf_z = electrons[0]->par()->getPz(); 
	kf   = electrons[0]->par()->getP();
	kf_v2 = sqrt(kf_x*kf_x + kf_y*kf_y + kf_z*kf_z);
	
	// knocked-out protons 3-momentum
	pf_x = protons[0]->par()->getPx(); 
	pf_y = protons[0]->par()->getPy();
	pf_z = protons[0]->par()->getPz(); 
	pf   = protons[0]->par()->getP();
	pf_v2 = sqrt(pf_x*pf_x + pf_y*pf_y + pf_z*pf_z);
	
	// set 4-momenta of particle of interest
	p4_beam.SetXYZM(0., 0., Eb, me);
	p4_target.SetXYZM(0.,0.,0., MP);
	p4_electron.SetXYZM(kf_x, kf_y, kf_z, me);
	p4_hadron.SetXYZM(pf_x, pf_y, pf_z, MP);		
	
	// 4-momentum trasnferred
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
	th_e = electrons[0]->par()->getTheta()*TMath::RadToDeg();
	th_e_v2 = p4_beam.Angle( p4_electron.Vect() ) *TMath::RadToDeg();
	
	// define recoil (usually undetected) system kinematics
	p4_recoil = p4_q + p4_target - p4_hadron;  

	xangle = p4_hadron.Angle( p4_electron.Vect() );
	th_h = protons[0]->par()->getTheta()*TMath::RadToDeg();
	th_h_v2 = xangle - th_e_v2;
	
	W = -1000.;
	MM = -1000.;

	MM2 = p4_recoil.M2();  //recoil ("MISSING") mass squared of the system

	if(MM2>0){
	  MM = sqrt(MM2); // INVARIANT MASS
	}
	
	W2 = (p4_q + p4_target).M2(); //invariant mass squared 
	if(W2>0){
	  W = sqrt(W2); // INVARIANT MASS
	}
	
	
	// Fill Histograms 

	// electron kinematics
	H_kf  ->Fill(kf);
	H_kf_v2  ->Fill(kf_v2); 
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
	H_the_v2 ->Fill(th_e_v2);
	H_thq ->Fill(th_q);
	H_phq ->Fill(ph_q);	
	H_W   ->Fill(W);   
	H_W2  ->Fill(W2);  	

	// hadron kinematics
	H_pf  ->Fill(pf);
	H_pf_v2  ->Fill(pf_v2);
	H_pfx  ->Fill(pfx);
	H_pfy  ->Fill(pfy);
	H_pfz  ->Fill(pfz);
	H_thh  ->Fill(th_h);
	H_thh_v2  ->Fill(th_h_v2);
	H_MM  ->Fill(MM);  
	H_MM2 ->Fill(MM2); 

	
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
  
  H_kf  ->Write();
  H_kf_v2  ->Write(); 
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
  H_the_v2 ->Write();
  H_thq ->Write();
  H_phq ->Write();
  H_W   ->Write();   
  H_W2  ->Write();

  H_pf  ->Write();
  H_pf_v2 -> Write();
  H_pfx  ->Write();
  H_pfy  ->Write();
  H_pfz  ->Write();
  H_thh  ->Write();
  H_thh_v2 ->Write();
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
