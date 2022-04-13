#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
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


void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
	      rp->par()->getPz(),p4.M());

}

void CLAS12SkimmerTree_simple(TString inFile = "", TString outputFile = "", double beamE = 0){
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();


  Double_t mD = 1.8756;


  Double_t q2 = 0;
  Double_t w2 = 0;
  Double_t px = 0;
  Double_t py = 0;
  Double_t pz = 0;
  Double_t vx = 0;
  Double_t vy = 0;
  Double_t vz = 0;
  Double_t chi2pid = 0;

  /////////////////////////////////////
  TFile *tree_file = new TFile(outputFile+".root","RECREATE");

  TTree *tree = new TTree("tree","Low Energy Data");
  tree->Branch("px",&px,"px/D");
  tree->Branch("py",&py,"px/D");
  tree->Branch("pz",&pz,"px/D");
  tree->Branch("chi2pid",&chi2pid,"chi2pid/D");


  //initialising clas12writer with path to output file
  clas12databases::SetCCDBLocalConnection("ccdb.sqlite");
  //clas12databases::SetRCDBRootConnection("rcdb.root");

  clas12root::HipoChain chain;
  chain.Add(inFile);
  chain.SetReaderTags({0});
  auto config_c12=chain.GetC12Reader();
  chain.db()->turnOffQADB();

  auto& c12=chain.C12ref();
  
  // How do we actually get the true beam energy from Hall B?
  // Where are the results from the Hall B Energy Measurement?
  auto beam_energy = 6.0;
  
  auto db=TDatabasePDG::Instance();
  double mass_p = db->GetParticle(2212)->Mass();

  //some particles
  TLorentzVector beam(0,0,beam_energy, beam_energy);
  TLorentzVector target(0,0,0,mD);
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector pr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector pip(0,0,0,db->GetParticle(211)->Mass());
  TLorentzVector pim(0,0,0,db->GetParticle(-211)->Mass());
   
  gBenchmark->Start("timer");
 
  int counter = 0;
  while(chain.Next())
    {

      auto particles = c12->getDetParticles(); //particles is now a std::vector of particles for this event

      beam.SetE(beam_energy);
      beam.SetPz(beam_energy);

      if( beamE > 1e-4)
	{
	  if( counter == 0)
	    cout<<"Beam energy set manually: "<<beamE<<endl;
	  beam.SetE(beamE);
	  beam.SetPz(beamE);
	}

    q2 = 0;
    w2 = 0.;
    px = 0.;
    py = 0.;
    pz = 0.;
    vx = 0.;
    vy = 0.;
    vz = 0.;
    chi2pid = 0.;

    for(auto& p : particles)
      {

	 px = p->par()->getPx();
	 py = p->par()->getPy();
	 pz = p->par()->getPz();
	 vx = p->par()->getVx();
	 vy = p->par()->getVx();
	 vz = p->par()->getVx();
	 chi2pid = p->par()->getChi2Pid();

	 /*
	p->par()->getVt();
	p->par()->getBeta();
	p->par()->getChi2Pid();
	p->par()->getCharge();
	p->par()->getStatus();
	p->par()->getP(); //magnitude of momentum
	 */

      }


    tree->Fill();
    }
  
  tree_file->cd();
  tree->Write();
  tree_file->Close();
  
  gBenchmark->Stop("timer");
  gBenchmark->Print("timer");
  
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count()<< " read events = "<<counter<<endl;
  }
