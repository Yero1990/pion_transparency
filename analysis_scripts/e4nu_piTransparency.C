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

  // Record start time
  auto start = std::chrono::high_resolution_clock::now();

  // output file (to write leaf variables)
  TFile *outROOT = new TFile("e4nu_test_cyero.root","RECREATE");

  // load databases (check with Justin where are these located)
  clas12databases::SetCCDBLocalConnection("ccdb.sqlite");
  clas12databases::SetRCDBRootConnection("rcdb.root");

  //LH2 target (we can use to construct W and check.
  TString inFile = "/work/clas12/rg-m/LH2/prod1.4/dst/recon/015024/rec_clas_015024.evio.000*.hipo";

  // Chain multiple .hipo files
  clas12root::HipoChain chain;
  chain.Add(inFile);
  chain.SetReaderTags({0});
  auto config_c12 = chain.GetC12Reader();
  chain.db()->turnOffQADB();

  auto& c12 = chain.C12ref();
  auto& rcdbData= config_c12->rcdb()->current();//struct with all relevent rcdb values        
  auto Eb = 5.98636; // GeV  | this command does NOT work (it gives 0) : rcdbData.beam_energy/1000;
  
  cout << "Beam energy is " << Eb << endl;
  
  // Load particle database from PDG to get mass
  auto db=TDatabasePDG::Instance();
  double MP = db->GetParticle(2212)->Mass();  
  double me = db->GetParticle(11)->Mass();

  //Define 4-momentum particles
  TLorentzVector p4_beam;
  TLorentzVector p4_target;
  TLorentzVector p4_electron;
  TLorentzVector p4_proton;
  TLorentzVector p4_q; 
  TLorentzVector p4_recoil; // recoil system 4-momenta (usually, undetected)


  //Define Histograms
  TH1F *H_the     = new TH1F("H_the", "Electron Scattering Angle, #theta_{e}", 200, 0.5, 180);
  TH1F *H_kf      = new TH1F("H_kf", "Final e^{-} Momentum", 100, 0.5, 6);
  TH1F *H_W       = new TH1F("H_W", "Invariant Mass, W", 100, 0.1, 5); 
  TH1F *H_W2      = new TH1F("H_W2", "Invariant Mass, W^{2}", 100, -5, 5); 
  TH1F *H_Q2      = new TH1F("H_Q2","4-Momentum Transfer, Q^{2}", 100, 0.1, 5); 
  TH1F *H_xbj     = new TH1F("H_xbj", "x-Bjorken", 100, 0.2, 2);  
  TH1F *H_nu      = new TH1F("H_nu","Energy Transfer, #nu", 100, 0.2, 6); 
  TH1F *H_q       = new TH1F("H_q", "3-Momentum Transfer, |#vec{q}|", 100, 0.2, 6);
  TH1F *H_qx      = new TH1F("H_qx", "|#vec{q}_{x}|", 100, 0.2, 6);
  TH1F *H_qy      = new TH1F("H_qy", "|#vec{q}_{y}|", 100, 0.2, 6);
  TH1F *H_qz      = new TH1F("H_qz", "|#vec{q}_{z}|", 100, 0.2, 6);
  TH1F *H_thq     = new TH1F("H_thq", "In-Plane Angle w.r.t +z(lab), #theta_{q}", 100, 0, 180); 
  TH1F *H_phq     = new TH1F("H_phq", "Out-of-Plane Angle w.r.t +z(lab), #phi_{q}", 100, -180, 180);
  	     
  TH1F *H_MM      = new TH1F("H_MM","Missing Mass, M_{miss}", 100, -5, 5 );        
  TH1F *H_MM2     = new TH1F("H_MM2","Missing Mass Squared, M^{2}_{miss}", 100, -5, 5); 


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

      
      /*
      double th_e = electrons[0]->getTheta();  //radians
      
      double Ef = sqrt( pow(me,2) + pow(electrons[0]->par()->getP(), 2) );
      double Q2 = 4.*Eb*Ef*pow(sin(th_e/2.), 2);
      double nu = Eb - Ef;
      double W = pow(MP,2) + 2*MP*nu - Q2;

      if (electrons.size()==1 && protons.size()==1) {
	cout << "Ef = " << Ef << endl;
	cout << "Q2 = " << Q2 << endl;
	cout << "nu = " << nu << endl;
	cout << "W = " << W << endl;
      }       
      */
      
      /*
      cout << "particles size: " << particles.size() << endl;
      cout << "electrons size: " << electrons.size() << endl;
      cout << "protons size: "   << protons.size() << endl;

      cout << "proton mass: " << db->GetParticle(2212)->Mass() << endl;
      cout << "electron mass: " << db->GetParticle(11)->Mass() << endl;  
      */


      // select ONLY 1 electron and 1 proton in final state: H(e,e'p) -> but about about everythin else?
      // is it:  1 e- + 1 p + whatever else is detected?
      /*
      if (electrons.size()==1 && protons.size()==1){
	
	double th_e = electrons[0]->getTheta();  //radians                                                                                                             
                                                                                                                                                                     
	double Ef = sqrt( pow(me,2) + pow(electrons[0]->par()->getP(), 2) );                                                                                           
	double Q2 = 4.*Eb*Ef*pow(sin(th_e/2.), 2);                                                                                                                     
	double nu = Eb - Ef;                                                                                                                                           
	double W = pow(MP,2) + 2*MP*nu - Q2;   
	
	cout << "Ef = " << Ef << endl;                                                                                                                               
        cout << "Q2 = " << Q2 << endl;                                                                                                                               
        cout << "nu = " << nu << endl;                                                                                                                               
        cout << "W = " << W << endl;   

	cout << "ELECTRONS == 1 ! ! ! ! " << endl;
	// loop over pth element in particles array
	for(int i=0;i<particles.size();i++){
	  
	  cout << "particle pid --> " << particles[i]->par()->getPid() << endl;	 
	  cout << "particle time --> " << particles[i]->getTime() << endl;
	  cout << "particle path --> " <<particles[i]->getPath() << endl;
	  cout << "particle sector -->"	<< particles[i]->getSector() << endl;
	  cout << "particle region --> " <<particles[i]->getRegion() << endl;
	  cout << "particle theta --> " << particles[i]->getTheta() << endl;
	  cout << "electron theta -->" << electrons[0]->getTheta() << endl;
	  cout << "particle phi --> " << particles[i]->getPhi() << endl;

	  cout << "px = " << particles[i]->par()->getPx() << endl;
	  cout << "px_elec = " << electrons[0]->par()->getPx() << endl; 
	  cout << "px = " << particles[i]->par()->getPy() << endl; 
	  cout << "px = " << particles[i]->par()->getPz() << endl; 
	  cout << "p = " << particles[i]->par()->getP() << endl; 

	}
      }
      
      */
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


  outROOT->cd();
  
  // Write Histograms 
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
