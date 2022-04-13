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
  TFile *tree_file = new TFile("e4nu_test_cyero.root","RECREATE");

  // load databases (check with Justin where are these located)
  clas12databases::SetCCDBLocalConnection("ccdb.sqlite");
  clas12databases::SetRCDBRootConnection("rcdb.root");

  //LH2 target (we can use to construct W and check.
  TString inFile = "/work/clas12/rg-m/LH2/prod1.4/dst/recon/015024/rec_clas_015024.evio.00001.hipo";

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


  //Define some useful variables
  int evt_cnt=0;
  

  // Loop over all events
  while(chain.Next())
    {
      cout << "===========" << endl;
      cout << "event #: " << evt_cnt << endl;
      cout << "===========" << endl;

      
      // TESTING / DEBUGGING
      //cout << "beamE = " << beam_energy << endl;

      auto particles = c12->getDetParticles();
      
      // get particles by type (returns vector array of particles, for example, an event might have 3 protons)
      auto electrons = c12->getByID(11);
      auto protons   = c12->getByID(2212);

      
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
