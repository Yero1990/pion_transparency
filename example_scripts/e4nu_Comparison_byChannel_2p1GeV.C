/*////////////////////////////////////////////////////
Original code by J. Estee, MIT
Code updated for comparisons by J. L. Barrow, MIT/TAU
////////////////////////////////////////////////////*/

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

void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp)
{
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
	      rp->par()->getPz(),p4.M());

}

//MC comparison added by J. L. Barrow (JLB)
void e4nu_Comparison_byChannel_2p1GeV(TString inFile_data = "", TString inFile_MC = "", TString outputFile = "", double beamE = 0., double omega_threshold = 0., double vz_min = -10., double vz_max = 10.)
{
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();

  ofstream outfile_byChannel("outfile_byChannel.txt");
  ofstream outfile_byProcess("outfile_byProcesswTopology.txt");

  //initialising clas12writer with path to output file
  //clas12databases::SetCCDBLocalConnection("ccdb.sqlite");
  //clas12databases::SetRCDBRootConnection("rcdb.root");
  //Changed by JLB
  clas12root::HipoChain chain_data; clas12root::HipoChain chain_MC;
  //JLB
  chain_data.Add(inFile_data); chain_MC.Add(inFile_MC);
  //JLB
  chain_data.SetReaderTags({0}); chain_MC.SetReaderTags({0});
  auto config_c12_data=chain_data.GetC12Reader(); auto config_c12_MC=chain_MC.GetC12Reader();
  auto& c12_data=chain_data.C12ref(); auto& c12_MC=chain_MC.C12ref();
  //struct with all relevent rcdb values
  //auto& rcdbData= config_c12_data->rcdb()->current(); auto& rcdbMC= config_c12_MC->rcdb()->current();
  auto beam_energy_data = beamE, beam_energy_MC = beamE;
  auto omega_threshold_data = omega_threshold, omega_threshold_MC = omega_threshold;
  cout << "Beam energy in data is: " << beam_energy_data << "GeV, and the beam energy in MC is: " << beam_energy_MC << "GeV" << endl;
  
  auto db=TDatabasePDG::Instance();
  chain_data.db()->turnOffQADB(); chain_MC.db()->turnOffQADB();
  double mass_p = db->GetParticle(2212)->Mass();
  double mass_n = db->GetParticle(2112)->Mass();
  double mass_pip = db->GetParticle(211)->Mass();
  double mass_pim = db->GetParticle(-211)->Mass();

  //Make some particles
  TLorentzVector beam_data(0,0,beam_energy_data, beam_energy_data);
  TLorentzVector beam_MC(0,0,beam_energy_MC, beam_energy_MC);
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector pr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector nr(0,0,0,db->GetParticle(2112)->Mass());
  TLorentzVector pip(0,0,0,db->GetParticle(211)->Mass());
  TLorentzVector pim(0,0,0,db->GetParticle(-211)->Mass());
   
  gBenchmark->Start("timer");
 
  //JLB: Stacked versions of 1D histograms (show with "nostack" option)
  //These stacked histograms will store the 1D histograms as they are made
  THStack *mult_p_stacked = new THStack("mult_p_stacked","Proton Multiplicity Comparison;Multiplicity of p;Counts");
  THStack *mult_n_stacked = new THStack("mult_n_stacked","Neutron Multiplicity Comparison;Multiplicity of n;Counts");
  THStack *mult_pip_stacked = new THStack("mult_pip_stacked","#pi^{+} Multiplicity Comparison;Multiplicity of #pi^{+};Counts");
  THStack *mult_pim_stacked = new THStack("mult_pim_stacked","#pi^{-} Multiplicity Comparison;Multiplicity of #pi^{-};Counts");
  //THStack *q2_mc_h_stacked = new THStack("q2_mc_h_stacked","Q^{2} Comparison;Q^{2};Counts");
  THStack *q2_stacked = new THStack("q2_stacked","Q^{2} Data to GENIE Comparison;Q^{2};Counts");
  THStack *p_chi2_stacked = new THStack("p_chi2_stacked","Proton #chi^{2} Data to GENIE Comparison;Proton #chi^{2};Counts");
  THStack *n_chi2_stacked = new THStack("n_chi2_stacked","Neutron #chi^{2} Data to GENIE Comparison;Neutron #chi^{2};Counts");
  THStack *pip_chi2_stacked = new THStack("pip_chi2_stacked","#pi^{+} #chi^{2} Data to GENIE Comparison;#pi^{+} #chi^{2};Counts");
  THStack *pim_chi2_stacked = new THStack("pim_chi2_stacked","#pi^{-} #chi^{2} Data to GENIE Comparison;#pi^{-} #chi^{2};Counts");
  THStack *energy_transfer_0_stacked = new THStack("energy_transfer_0_stacked","Energy transfer  5 < #theta < 10 Data to GENIE Comparison; #omega (GeV); Counts");
  THStack *energy_transfer_1_stacked = new THStack("energy_transfer_1_stacked","Energy transfer  10 < #theta < 15 Data to GENIE Comparison; #omega (GeV); Counts");
  THStack *energy_transfer_2_stacked = new THStack("energy_transfer_2_stacked","Energy transfer  15 < #theta < 25 Data to GENIE Comparison; #omega (GeV); Counts");
  THStack *energy_transfer_3_stacked = new THStack("energy_transfer_3_stacked","Energy transfer  25 < #theta < 35 Data to GENIE Comparison; #omega (GeV); Counts");
  THStack *invariant_mass_stacked = new THStack("invariant_mass_stacked","Invariant Mass Data to GENIE Comparison of p+#pi^{+} System;Mass (GeV/c^{2});Counts");
  THStack *el_vz_h_stacked = new THStack("el_vz_h_stacked","Z-Vertex Data to GENIE Comparison;Z-vertex (cm);Counts");
  THStack *W_stacked = new THStack("W_stacked","Invariant Mass, W, Data to GENIE Comparison;W (GeV/c^{2});Counts");
  THStack *theta_stacked = new THStack("theta_stacked","Angular Data to GENIE Comparison;#theta ({}^{#circ});Counts");
  THStack *omega_stacked = new THStack("omega_stacked","Energy Transfer Data to GENIE Comparison;#omega (GeV);Counts");
  THStack *p_proton_stacked = new THStack("p_proton_stacked","Proton Momentum Data to GENIE Comparison;Momentum (GeV/c);Counts");
  THStack *p_miss_stacked = new THStack("p_miss_stacked","Missing Momentum Data to GENIE Comparison;Missing Momentum (GeV/c);Counts");
  THStack *p_miss_theta_stacked = new THStack("p_miss_theta_stacked","Missing Momentum Angle #theta Comparison;#theta;Counts");
  THStack *xb_stacked = new THStack("xb_stacked","Bjorken-x Data to GENIE Comparison;x_{B};Counts");
  THStack *y_stacked = new THStack("y_stacked","Bjorken-y Scaling Variable Data to GENIE Comparison;y = -q + #sqrt{#omega^{2} + 2 #omega M};Counts");

  //Stacked 1D histograms for data comparisons by channel
  THStack *q2_byChannel_stacked = new THStack("q2_byChannel_stacked","Data: Q^{2} Distributions by Reconstructed Channel;Q^{2};Counts");
  THStack *totrecoe_byChannel_stacked = new THStack("totrecoe_byChannel_stacked","Data: Calorimetric Energy Distributions by Reconstructed Channel;Calorimetric Energy (GeV);Counts");
  THStack *W_byChannel_stacked = new THStack("W_byChannel_stacked","Data: Invariant Mass, W, Distributions by Reconstructed Channel;W (GeV/c^{2]);Counts");

  //Histograms for stacking by DC sector
  //W, Invariant mass
  THStack *W_bySector_data_stacked = new THStack("W_bySector_data_stacked","Data: Invariant Mass, W, Distributions by DC Sector;W (GeV/c^{2]);Counts");
  TH1D *W_bySector_data = new TH1D("W_bySector_data","Data, All Sectors;W (GeV/c^{2};Counts",400,0.5,2.2);
  TH1D *W_S1_data = new TH1D("W_S1_data","Datax6, Sector #1;W (GeV/c^{2};Counts",400,0.5,2.2);
  TH1D *W_S2_data = new TH1D("W_S2_data","Datax6, Sector #2;W (GeV/c^{2};Counts",400,0.5,2.2);
  TH1D *W_S3_data = new TH1D("W_S3_data","Datax6, Sector #3;W (GeV/c^{2};Counts",400,0.5,2.2);
  TH1D *W_S4_data = new TH1D("W_S4_data","Datax6, Sector #4;W (GeV/c^{2};Counts",400,0.5,2.2);
  TH1D *W_S5_data = new TH1D("W_S5_data","Datax6, Sector #5;W (GeV/c^{2};Counts",400,0.5,2.2);
  TH1D *W_S6_data = new TH1D("W_S6_data","Datax6, Sector #6;W (GeV/c^{2};Counts",400,0.5,2.2);
  //Q^2
  THStack *q2_bySector_data_stacked = new THStack("q2_bySector_data_stacked","Data: Q^{2} Distributions by DC Sector;Q^{2};Counts");
  TH1D *q2_bySector_data = new TH1D("q2_bySector_data","Data, All Sectors;Q^{2};Counts",400,0.,1.2);
  TH1D *q2_S1_data = new TH1D("q2_S1_data","Datax6, Sector #1;Q^{2};Counts",400,0.,1.2);
  TH1D *q2_S2_data = new TH1D("q2_S2_data","Datax6, Sector #2;Q^{2};Counts",400,0.,1.2);
  TH1D *q2_S3_data = new TH1D("q2_S3_data","Datax6, Sector #3;Q^{2};Counts",400,0.,1.2);
  TH1D *q2_S4_data = new TH1D("q2_S4_data","Datax6, Sector #4;Q^{2};Counts",400,0.,1.2);
  TH1D *q2_S5_data = new TH1D("q2_S5_data","Datax6, Sector #5;Q^{2};Counts",400,0.,1.2);
  TH1D *q2_S6_data = new TH1D("q2_S6_data","Datax6, Sector #6;Q^{2};Counts",400,0.,1.2);
  //Calorimetric Energy, 1pXn0pi
  THStack *totrecoe_1pXn0pi_bySector_data_stacked = new THStack("totrecoe_1pXn0pi_bySector_data_stacked","Data: Total Calorimetric Energy of 1pXn0#pi^{#pm} by DC Sector;Calorimetric Energy (GeV);Counts");
  TH1D *totrecoe_1pXn0pi_bySector_data = new TH1D("totrecoe_bySector_data","Data, All Sectors;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn0pi_S1_data = new TH1D("totrecoe_S1_data","Datax6, Sector #1;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn0pi_S2_data = new TH1D("totrecoe_S2_data","Datax6, Sector #2;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn0pi_S3_data = new TH1D("totrecoe_S3_data","Datax6, Sector #3;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn0pi_S4_data = new TH1D("totrecoe_S4_data","Datax6, Sector #4;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn0pi_S5_data = new TH1D("totrecoe_S5_data","Datax6, Sector #5;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn0pi_S6_data = new TH1D("totrecoe_S6_data","Datax6, Sector #6;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  //Calorimetric Energy, 1pXn1pim
  THStack *totrecoe_1pXn1pim_bySector_data_stacked = new THStack("totrecoe_1pXn1pim_bySector_data_stacked","Data: Total Calorimetric Energy of 1pXn1#pi^{-} by DC Sector;Calorimetric Energy (GeV);Counts");
  TH1D *totrecoe_1pXn1pim_bySector_data = new TH1D("totrecoe_bySector_data","Data, All Sectors;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn1pim_S1_data = new TH1D("totrecoe_S1_data","Datax6, Sector #1;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn1pim_S2_data = new TH1D("totrecoe_S2_data","Datax6, Sector #2;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn1pim_S3_data = new TH1D("totrecoe_S3_data","Datax6, Sector #3;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn1pim_S4_data = new TH1D("totrecoe_S4_data","Datax6, Sector #4;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn1pim_S5_data = new TH1D("totrecoe_S5_data","Datax6, Sector #5;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn1pim_S6_data = new TH1D("totrecoe_S6_data","Datax6, Sector #6;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  //Calorimetric Energy, 1pXn2pi
  THStack *totrecoe_1pXn2pi_bySector_data_stacked = new THStack("totrecoe_1pXn2pi_bySector_data_stacked","Data: Total Calorimetric Energy of 1pXn2#pi^{#pm} by DC Sector;Calorimetric Energy (GeV);Counts");
  TH1D *totrecoe_1pXn2pi_bySector_data = new TH1D("totrecoe_bySector_data","Data, All Sectors;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn2pi_S1_data = new TH1D("totrecoe_S1_data","Datax6, Sector #1;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn2pi_S2_data = new TH1D("totrecoe_S2_data","Datax6, Sector #2;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn2pi_S3_data = new TH1D("totrecoe_S3_data","Datax6, Sector #3;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn2pi_S4_data = new TH1D("totrecoe_S4_data","Datax6, Sector #4;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn2pi_S5_data = new TH1D("totrecoe_S5_data","Datax6, Sector #5;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn2pi_S6_data = new TH1D("totrecoe_S6_data","Datax6, Sector #6;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  //Calorimetric Energy, 0p1n2pi
  THStack *totrecoe_0p1n2pi_bySector_data_stacked = new THStack("totrecoe_0p1n2pi_bySector_data_stacked","Data: Total Calorimetric Energy of 0p1n2#pi^{#pm} by DC Sector;Calorimetric Energy (GeV);Counts");
  TH1D *totrecoe_0p1n2pi_bySector_data = new TH1D("totrecoe_bySector_data","Data, All Sectors;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_0p1n2pi_S1_data = new TH1D("totrecoe_S1_data","Datax6, Sector #1;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_0p1n2pi_S2_data = new TH1D("totrecoe_S2_data","Datax6, Sector #2;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_0p1n2pi_S3_data = new TH1D("totrecoe_S3_data","Datax6, Sector #3;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_0p1n2pi_S4_data = new TH1D("totrecoe_S4_data","Datax6, Sector #4;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_0p1n2pi_S5_data = new TH1D("totrecoe_S5_data","Datax6, Sector #5;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_0p1n2pi_S6_data = new TH1D("totrecoe_S6_data","Datax6, Sector #6;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  //W vs. Theta
  TH2D *W_theta_S1_data = new TH2D("W_theta_S1_data","Datax6, Sector #1: W vs #theta;W (GeV/c^{2});#theta",400,0.,2.,400,0.,20.);
  TH2D *W_theta_S2_data = new TH2D("W_theta_S2_data","Datax6, Sector #2: W vs #theta;W (GeV/c^{2});#theta",400,0.,2.,400,0.,20.);
  TH2D *W_theta_S3_data = new TH2D("W_theta_S3_data","Datax6, Sector #3: W vs #theta;W (GeV/c^{2});#theta",400,0.,2.,400,0.,20.);
  TH2D *W_theta_S4_data = new TH2D("W_theta_S4_data","Datax6, Sector #4: W vs #theta;W (GeV/c^{2});#theta",400,0.,2.,400,0.,20.);
  TH2D *W_theta_S5_data = new TH2D("W_theta_S5_data","Datax6, Sector #5: W vs #theta;W (GeV/c^{2});#theta",400,0.,2.,400,0.,20.);
  TH2D *W_theta_S6_data = new TH2D("W_theta_S6_data","Datax6, Sector #6: W vs #theta;W (GeV/c^{2});#theta",400,0.,2.,400,0.,20.);
  //Now for GENIE simulation
  //W, Invariant mass
  THStack *W_bySector_MC_stacked = new THStack("W_bySector_MC_stacked","GENIE: Invariant Mass, W, Distributions by DC Sector;W (GeV/c^{2]);Counts");
  TH1D *W_bySector_MC = new TH1D("W_bySector_MC","GENIE, All Sectors;W (GeV/c^{2};Counts",400,0.5,2.2);
  TH1D *W_S1_MC = new TH1D("W_S1_MC","GENIEx6, Sector #1;W (GeV/c^{2};Counts",400,0.5,2.2);
  TH1D *W_S2_MC = new TH1D("W_S2_MC","GENIEx6, Sector #2;W (GeV/c^{2};Counts",400,0.5,2.2);
  TH1D *W_S3_MC = new TH1D("W_S3_MC","GENIEx6, Sector #3;W (GeV/c^{2};Counts",400,0.5,2.2);
  TH1D *W_S4_MC = new TH1D("W_S4_MC","GENIEx6, Sector #4;W (GeV/c^{2};Counts",400,0.5,2.2);
  TH1D *W_S5_MC = new TH1D("W_S5_MC","GENIEx6, Sector #5;W (GeV/c^{2};Counts",400,0.5,2.2);
  TH1D *W_S6_MC = new TH1D("W_S6_MC","GENIEx6, Sector #6;W (GeV/c^{2};Counts",400,0.5,2.2);
  //Q^2
  THStack *q2_bySector_MC_stacked = new THStack("q2_bySector_MC_stacked","GENIE: Q^{2} Distributions by DC Sector;Q^{2};Counts");
  TH1D *q2_bySector_MC = new TH1D("q2_bySector_MC","GENIE, All Sectors;Q^{2};Counts",400,0.,1.2);
  TH1D *q2_S1_MC = new TH1D("q2_S1_MC","GENIEx6, Sector #1;Q^{2};Counts",400,0.,1.2);
  TH1D *q2_S2_MC = new TH1D("q2_S2_MC","GENIEx6, Sector #2;Q^{2};Counts",400,0.,1.2);
  TH1D *q2_S3_MC = new TH1D("q2_S3_MC","GENIEx6, Sector #3;Q^{2};Counts",400,0.,1.2);
  TH1D *q2_S4_MC = new TH1D("q2_S4_MC","GENIEx6, Sector #4;Q^{2};Counts",400,0.,1.2);
  TH1D *q2_S5_MC = new TH1D("q2_S5_MC","GENIEx6, Sector #5;Q^{2};Counts",400,0.,1.2);
  TH1D *q2_S6_MC = new TH1D("q2_S6_MC","GENIEx6, Sector #6;Q^{2};Counts",400,0.,1.2);
  //Calorimetric Energy, 1pXn0pi
  THStack *totrecoe_1pXn0pi_bySector_MC_stacked = new THStack("totrecoe_1pXn0pi_bySector_MC_stacked","GENIE: Total Calorimetric Energy of 1pXn0#pi^{#pm} by DC Sector;Calorimetric Energy (GeV);Counts");
  TH1D *totrecoe_1pXn0pi_bySector_MC = new TH1D("totrecoe_bySector_MC","GENIE, All Sectors;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn0pi_S1_MC = new TH1D("totrecoe_S1_MC","GENIEx6, Sector #1;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn0pi_S2_MC = new TH1D("totrecoe_S2_MC","GENIEx6, Sector #2;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn0pi_S3_MC = new TH1D("totrecoe_S3_MC","GENIEx6, Sector #3;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn0pi_S4_MC = new TH1D("totrecoe_S4_MC","GENIEx6, Sector #4;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn0pi_S5_MC = new TH1D("totrecoe_S5_MC","GENIEx6, Sector #5;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn0pi_S6_MC = new TH1D("totrecoe_S6_MC","GENIEx6, Sector #6;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  //Calorimetric Energy, 1pXn1pim
  THStack *totrecoe_1pXn1pim_bySector_MC_stacked = new THStack("totrecoe_1pXn1pim_bySector_MC_stacked","GENIE: Total Calorimetric Energy of 1pXn1#pi^{-} by DC Sector;Calorimetric Energy (GeV);Counts");
  TH1D *totrecoe_1pXn1pim_bySector_MC = new TH1D("totrecoe_bySector_MC","GENIE, All Sectors;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn1pim_S1_MC = new TH1D("totrecoe_S1_MC","GENIEx6, Sector #1;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn1pim_S2_MC = new TH1D("totrecoe_S2_MC","GENIEx6, Sector #2;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn1pim_S3_MC = new TH1D("totrecoe_S3_MC","GENIEx6, Sector #3;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn1pim_S4_MC = new TH1D("totrecoe_S4_MC","GENIEx6, Sector #4;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn1pim_S5_MC = new TH1D("totrecoe_S5_MC","GENIEx6, Sector #5;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn1pim_S6_MC = new TH1D("totrecoe_S6_MC","GENIEx6, Sector #6;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  //Calorimetric Energy, 1pXn2pi
  THStack *totrecoe_1pXn2pi_bySector_MC_stacked = new THStack("totrecoe_1pXn2pi_bySector_MC_stacked","GENIE: Total Calorimetric Energy of 1pXn2#pi^{#pm} by DC Sector;Calorimetric Energy (GeV);Counts");
  TH1D *totrecoe_1pXn2pi_bySector_MC = new TH1D("totrecoe_bySector_MC","GENIE, All Sectors;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn2pi_S1_MC = new TH1D("totrecoe_S1_MC","GENIEx6, Sector #1;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn2pi_S2_MC = new TH1D("totrecoe_S2_MC","GENIEx6, Sector #2;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn2pi_S3_MC = new TH1D("totrecoe_S3_MC","GENIEx6, Sector #3;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn2pi_S4_MC = new TH1D("totrecoe_S4_MC","GENIEx6, Sector #4;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn2pi_S5_MC = new TH1D("totrecoe_S5_MC","GENIEx6, Sector #5;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_1pXn2pi_S6_MC = new TH1D("totrecoe_S6_MC","GENIEx6, Sector #6;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  //Calorimetric Energy, 0p1n2pi
  THStack *totrecoe_0p1n2pi_bySector_MC_stacked = new THStack("totrecoe_0p1n2pi_bySector_MC_stacked","GENIE: Total Calorimetric Energy of 0p1n2#pi^{#pm} by DC Sector;Calorimetric Energy (GeV);Counts");
  TH1D *totrecoe_0p1n2pi_bySector_MC = new TH1D("totrecoe_bySector_MC","GENIE, All Sectors;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_0p1n2pi_S1_MC = new TH1D("totrecoe_S1_MC","GENIEx6, Sector #1;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_0p1n2pi_S2_MC = new TH1D("totrecoe_S2_MC","GENIEx6, Sector #2;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_0p1n2pi_S3_MC = new TH1D("totrecoe_S3_MC","GENIEx6, Sector #3;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_0p1n2pi_S4_MC = new TH1D("totrecoe_S4_MC","GENIEx6, Sector #4;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_0p1n2pi_S5_MC = new TH1D("totrecoe_S5_MC","GENIEx6, Sector #5;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  TH1D *totrecoe_0p1n2pi_S6_MC = new TH1D("totrecoe_S6_MC","GENIEx6, Sector #6;Calorimetric Energy (GeV);Counts",400,0.65,4.5);
  //W vs. Theta
  TH2D *W_theta_S1_MC = new TH2D("W_theta_S1_MC","GENIEx6, Sector #1: W vs #theta;W (GeV/c^{2});#theta",400,0.,2.,400,0.,20.);
  TH2D *W_theta_S2_MC = new TH2D("W_theta_S2_MC","GENIEx6, Sector #2: W vs #theta;W (GeV/c^{2});#theta",400,0.,2.,400,0.,20.);
  TH2D *W_theta_S3_MC = new TH2D("W_theta_S3_MC","GENIEx6, Sector #3: W vs #theta;W (GeV/c^{2});#theta",400,0.,2.,400,0.,20.);
  TH2D *W_theta_S4_MC = new TH2D("W_theta_S4_MC","GENIEx6, Sector #4: W vs #theta;W (GeV/c^{2});#theta",400,0.,2.,400,0.,20.);
  TH2D *W_theta_S5_MC = new TH2D("W_theta_S5_MC","GENIEx6, Sector #5: W vs #theta;W (GeV/c^{2});#theta",400,0.,2.,400,0.,20.);
  TH2D *W_theta_S6_MC = new TH2D("W_theta_S6_MC","GENIEx6, Sector #6: W vs #theta;W (GeV/c^{2});#theta",400,0.,2.,400,0.,20.);

  //1D histograms for stacking with data-only comparisons by reconstructed channel
  //Q^2 distributions
  TH1D *q2_byChannel_data = new TH1D("q2_byChannel_data","Data;Q^{2};Counts",400,0.,2.);
  TH1D *q2_0p1n0pi_data = new TH1D("q2_0p1n0pi_data","Data 0p1n0#pi^{#pm};Q^{2};Counts",400,0.,2.);
  TH1D *q2_0p1n0pi_MC = new TH1D("q2_0p1n0pi_MC","GENIE 0p1n0#pi^{#pm};Q^{2};Counts",400,0.,2.);//For Larry plot
  TH1D *q2_1pXn0pi_data = new TH1D("q2_1pXn0pi_data","Data 1pXn0#pi^{#pm};Q^{2};Counts",400,0.,2.);
  TH1D *q2_1pXn0pi_MC = new TH1D("q2_1pXn0pi_MC","GENIE 1pXn0#pi^{#pm};Q^{2};Counts",400,0.,2.);//For Larry plot
  TH1D *q2_1pXn1pim_data = new TH1D("q2_1pXn1pim_data","Data 1pXn1#pi^{-};Q^{2};Counts",400,0.,2.);
  TH1D *q2_0p1n1pip_data = new TH1D("q2_0p1n1pip_data","Data 0p1n1#pi^{+};Q^{2};Counts",400,0.,2.);
  TH1D *q2_2pXn0pi_data = new TH1D("q2_2pXn0pi_data","Data 2pXn0#pi^{#pm};Q^{2};Counts",400,0.,2.);
  TH1D *q2_1p1n0pi_data = new TH1D("q2_1p1n0pi_data","Data 1p1n0#pi^{#pm};Q^{2};Counts",400,0.,2.);
  TH1D *q2_1p1n1pi_data = new TH1D("q2_1p1n1pi_data","Data 1p1n1#pi^{#pm};Q^{2};Counts",400,0.,2.);
  TH1D *q2_1pXn2pi_data = new TH1D("q2_1pXn2pi_data","Data 1pXn2#pi^{#pm};Q^{2};Counts",400,0.,2.);
  TH1D *q2_0p1n2pi_data = new TH1D("q2_0p1n2pi_data","Data 0p1n2#pi^{#pm};Q^{2};Counts",400,0.,2.);
  
  //Calorimetric Energy distributions
  //TH1D *totrecoe_byChannel_data = new TH1D("totrecoe_byChannel_data","Data;Calorimetric Energy (GeV);Counts",100,0.,2.2);//This is ill defined
  TH1D *totrecoe_0p1n0pi_data = new TH1D("totrecoe_0p1n0pi_data","Data 0p1n0#pi^{#pm};Calorimetric Energy (GeV);Counts",100,0.,2.2);
  TH1D *totrecoe_0p1n0pi_MC = new TH1D("totrecoe_0p1n0pi_MC","GENIE 0p1n0#pi^{#pm};Calorimetric Energy (GeV);Counts",100,0.,2.2);//For Larry plot
  TH1D *totrecoe_1pXn0pi_data = new TH1D("totrecoe_1pXn0pi_data","Data 1pXn0#pi^{#pm};Calorimetric Energy (GeV);Counts",100,0.,2.2);
  TH1D *totrecoe_1pXn0pi_MC = new TH1D("totrecoe_1pXn0pi_MC","GENIE 1pXn0#pi^{#pm};Calorimetric Energy (GeV);Counts",100,0.,2.2);//For Larry plot
  TH1D *totrecoe_1pXn1pim_data = new TH1D("totrecoe_1pXn1pim_data","Data 1pXn1#pi^{-};Calorimetric Energy (GeV);Counts",100,0.,2.2);
  TH1D *totrecoe_0p1n1pip_data = new TH1D("totrecoe_0p1n1pip_data","Data 0p1n1#pi^{+};Calorimetric Energy (GeV);Counts",100,0.,2.2);
  TH1D *totrecoe_2pXn0pi_data = new TH1D("totrecoe_2pXn0pi_data","Data 2pXn0#pi^{#pm};Calorimetric Energy (GeV);Counts",100,0.,2.2);
  TH1D *totrecoe_1p1n0pi_data = new TH1D("totrecoe_1p1n0pi_data","Data 1p1n0#pi^{#pm};Calorimetric Energy (GeV);Counts",100,0.,2.2);
  TH1D *totrecoe_1p1n1pi_data = new TH1D("totrecoe_1p1n1pi_data","Data 1p1n1#pi^{#pm};Calorimetric Energy (GeV);Counts",100,0.,2.2);
  TH1D *totrecoe_1pXn2pi_data = new TH1D("totrecoe_1pXn2pi_data","Data 1pXn2#pi^{#pm};Calorimetric Energy (GeV);Counts",100,0.,2.2);
  TH1D *totrecoe_0p1n2pi_data = new TH1D("totrecoe_0p1n2pi_data","Data 0p1n2#pi^{#pm};Calorimetric Energy (GeV);Counts",100,0.,2.2);

  //W distributions
  TH1D *W_byChannel_data = new TH1D("W_byChannel_data","Data;W (GeV/c^{2};Counts",400,0.5,2.2);
  TH1D *W_0p1n0pi_data = new TH1D("W_0p1n0pi_data","Data 0p1n0#pi^{#pm};W (GeV/c^{2]);Counts",400,0.5,2.2);
  TH1D *W_0p1n0pi_MC = new TH1D("W_0p1n0pi_MC","GENIE 0p1n0#pi^{#pm};W (GeV/c^{2]);Counts",400,0.5,2.2);//For Larry plot
  TH1D *W_1pXn0pi_data = new TH1D("W_1pXn0pi_data","Data 1pXn0#pi^{#pm};W (GeV/c^{2]);Counts",400,0.5,2.2);
  TH1D *W_1pXn0pi_MC = new TH1D("W_1pXn0pi_MC","GENIE 1pXn0#pi^{#pm};W (GeV/c^{2]);Counts",400,0.5,2.2);//For Larry plot
  TH1D *W_1pXn1pim_data = new TH1D("W_1pXn1pim_data","Data 1pXn1#pi^{-};W (GeV/c^{2]);Counts",400,0.5,2.2);
  TH1D *W_0p1n1pip_data = new TH1D("W_0p1n1pip_data","Data 0p1n1#pi^{+};W (GeV/c^{2]);Counts",400,0.5,2.2);
  TH1D *W_2pXn0pi_data = new TH1D("W_2pXn0pi_data","Data 2pXn0#pi^{#pm};W (GeV/c^{2]);Counts",400,0.5,2.2);
  TH1D *W_1p1n0pi_data = new TH1D("W_1p1n0pi_data","Data 1p1n0#pi^{#pm};W (GeV/c^{2]);Counts",400,0.5,2.2);
  TH1D *W_1p1n1pi_data = new TH1D("W_1p1n1pi_data","Data 1p1n1#pi^{#pm};W (GeV/c^{2]);Counts",400,0.5,2.2);
  TH1D *W_1pXn2pi_data = new TH1D("W_1pXn2pi_data","Data 1pXn2#pi^{#pm};W (GeV/c^{2]);Counts",400,0.5,2.2);
  TH1D *W_0p1n2pi_data = new TH1D("W_0p1n2pi_data","Data 0p1n2#pi^{#pm};W (GeV/c^{2]);Counts",400,0.5,2.2);

  //JLB: 1D histograms for data
  //Data
  TH1D *mult_p_data = new TH1D("mult_p_data_data","Data;Multiplicity of p;Counts",5,-0.5,4.5);
  TH1D *mult_n_data = new TH1D("mult_n_data_data","Data;Multiplicity of n;Counts",5,-0.5,4.5);
  TH1D *mult_n_nocuts_data = new TH1D("mult_n_nocuts_data_data","Data: No #chi^{2} Cuts;Multiplicity of n;Counts",5,-0.5,4.5);
  TH1D *mult_pip_data = new TH1D("mult_pip_data","Data;Multiplicity of #pi^{+};Counts",5,-0.5,4.5);
  TH1D *mult_pim_data = new TH1D("mult_pim_data","Data;Multiplicity of #pi^{-};Counts",5,-0.5,4.5);
  //TH1D *q2_mc_h_data = new TH1D("q2_mc_h_data","Data: Q^{2} Distribution;Q^{2};Counts",400,0,2);
  TH1D *q2_data = new TH1D("q2_data","Data;Q^{2};Counts",400,0.,1.2);
  TH1D *p_chi2_cd_data = new TH1D("p_chi2_cd_data","Data: Proton CD #chi^{2} Distribution;Proton #chi^{2};Counts",100,-10.,10.);
  TH1D *n_chi2_cd_data = new TH1D("n_chi2_cd_data","Data: Neutron CD #chi^{2} Distribution;Neutron #chi^{2};Counts",100,-10.,10.);
  TH1D *n_chi2_nocuts_data = new TH1D("n_chi2_nocuts_data","Data: Neutron #chi^{2} Distribution, No Cuts;Neutron #chi^{2};Counts",100,-10.,10.);
  TH1D *pip_chi2_cd_data = new TH1D("pip_chi2_cd_data","Data: #pi^{+} CD #chi^{2} Distribution;#pi^{+} #chi^{2};Counts",100,-10.,10.);
  TH1D *pim_chi2_cd_data = new TH1D("pim_chi2_cd_data","Data: #pi^{-} CD #chi^{2} Distribution;#pi^{-} #chi^{2};Counts",100,-10.,10.);
  TH1D *p_chi2_fd_data = new TH1D("p_chi2_fd_data","Data: Proton FD #chi^{2} Distribution;Proton #chi^{2};Counts",100,-10.,10.);
  TH1D *n_chi2_fd_data = new TH1D("n_chi2_fd_data","Data: Neutron FD #chi^{2} Distribution;Neutron #chi^{2};Counts",100,-10.,10.);
  TH1D *pip_chi2_fd_data = new TH1D("pip_chi2_fd_data","Data: #pi^{+} FD #chi^{2} Distribution;#pi^{+} #chi^{2};Counts",100,-10.,10.);
  TH1D *pim_chi2_fd_data = new TH1D("pim_chi2_fd_data","Data: #pi^{-} FD #chi^{2} Distribution;#pi^{-} #chi^{2};Counts",100,-10.,10.);
  TH1D *energy_transfer_0_data = new TH1D("energy_transfer_0_data","Data;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *energy_transfer_1_data = new TH1D("energy_transfer_1_data","Data;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *energy_transfer_2_data = new TH1D("energy_transfer_2_data","Data;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *energy_transfer_3_data = new TH1D("energy_transfer_3_data","Data;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *invariant_mass_data = new TH1D("invariant_mass_data","Data;Mass (GeV/c^{2});Counts",100,0,2);
  TH1D *el_vz_h_data = new TH1D("el_vz_h_data","Data;Z-vertex (cm);Counts",1000,-10,10);
  TH1D *W_data = new TH1D("W_data","Data;W (GeV/c^{2};Counts",400,0.5,2.2);
  TH1D *theta_data = new TH1D("theta_data","Data;#theta;Counts",1000,0.,60.);
  TH1D *omega_data = new TH1D("omega_data","Data;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *p_proton_data = new TH1D("p_proton_data","Data;Momentum (GeV/c);Counts",1000,0.,4.);
  TH1D *p_miss_data = new TH1D("p_miss_data","Data;Missing Momentum (GeV/c);Counts",1000,-0.5,6.);
  TH1D *p_miss_theta_data = new TH1D("p_miss_theta_data","Data;#theta;Counts",1000,0.,200.);
  TH1D *xb_data = new TH1D("xb_data","Data;x_{B};Counts",1000,0.,3.);
  TH1D *y_data = new TH1D("y_data","Data;y = -q + #sqrt{#omega^{2} + 2 #omega M};Counts",1000,-0.5,1.);

  //JLB: 1D histograms for GENIE
  //MC Simulation
  TH1D *mult_p_MC = new TH1D("mult_p_data_MC","GENIE;Multiplicity of p;Counts",5,-0.5,4.5);
  TH1D *mult_n_MC = new TH1D("mult_n_data_MC","GENIE;Multiplicity of n;Counts",5,-0.5,4.5);
  TH1D *mult_n_nocuts_MC = new TH1D("mult_n_nocuts_MC","GENIE: No #chi^{2} Cuts;Multiplicity of n;Counts",5,-0.5,4.5);
  TH1D *mult_pip_MC = new TH1D("mult_pip_MC","GENIE;Multiplicity of #pi^{+};Counts",5,-0.5,4.5);
  TH1D *mult_pim_MC = new TH1D("mult_pim_MC","GENIE;Multiplicity of #pi^{-};Counts",5,-0.5,4.5);
  //TH1D *q2_mc_h_MC = new TH1D("q2_mc_h_MC","GENIE: Q^{2} Distribution;Q^{2};Counts",400,0,2);
  TH1D *q2_MC = new TH1D("q2_MC","GENIE;Q^{2};Counts",400,0.,1.2);
  TH1D *q2_qe_MC = new TH1D("q2_qe_MC","GENIE QE;Q^{2};Counts",400,0.,1.2);
  TH1D *q2_mec_MC = new TH1D("q2_mec_MC","GENIE MEC;Q^{2};Counts",400,0.,1.2);
  TH1D *q2_res_MC = new TH1D("q2_res_MC","GENIE RES;Q^{2};Counts",400,0.,1.2);
  TH1D *q2_dis_MC = new TH1D("q2_dis_MC","GENIE DIS;Q^{2};Counts",400,0.,1.2);
  TH1D *p_chi2_cd_MC = new TH1D("p_chi2_cd_MC","GENIE: Proton #chi^{2} Distribution;Proton #chi^{2};Counts",100,-10.,10.);
  TH1D *n_chi2_cd_MC = new TH1D("n_chi2_cd_MC","GENIE: Neutron #chi^{2} Distribution;Neutron #chi^{2};Counts",100,-10.,10.);\
  TH1D *n_chi2_nocuts_MC = new TH1D("n_chi2_nocuts_MC","GENIE: Neutron #chi^{2} Distribution, No Cuts;Neutron #chi^{2};Counts",100,-10.,10.);
  TH1D *pip_chi2_cd_MC = new TH1D("pip_chi2_cd_MC","GENIE: #pi^{+} #chi^{2} Distribution;#pi^{+} #chi^{2};Counts",100,-10.,10.);
  TH1D *pim_chi2_cd_MC = new TH1D("pim_chi2_cd_MC","GENIE: #pi^{-} #chi^{2} Distribution;#pi^{-} #chi^{2};Counts",100,-10.,10.);
  TH1D *p_chi2_fd_MC = new TH1D("p_chi2_fd_MC","GENIE: Proton #chi^{2} Distribution;Proton #chi^{2};Counts",100,-10.,10.);
  TH1D *n_chi2_fd_MC = new TH1D("n_chi2_fd_MC","GENIE: Neutron #chi^{2} Distribution;Neutron #chi^{2};Counts",100,-10.,10.);
  TH1D *pip_chi2_fd_MC = new TH1D("pip_chi2_fd_MC","GENIE: #pi^{+} #chi^{2} Distribution;#pi^{+} #chi^{2};Counts",100,-10.,10.);
  TH1D *pim_chi2_fd_MC = new TH1D("pim_chi2_fd_MC","GENIE: #pi^{-} #chi^{2} Distribution;#pi^{-} #chi^{2};Counts",100,-10.,10.);
  TH1D *energy_transfer_0_MC = new TH1D("energy_transfer_0_MC","GENIE;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *energy_transfer_1_MC = new TH1D("energy_transfer_1_MC","GENIE;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *energy_transfer_2_MC = new TH1D("energy_transfer_2_MC","GENIE;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *energy_transfer_3_MC = new TH1D("energy_transfer_3_MC","GENIE;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *energy_transfer_qe_0_MC = new TH1D("energy_transfer_qe_0_MC","GENIE QE;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *energy_transfer_qe_1_MC = new TH1D("energy_transfer_qe_1_MC","GENIE QE;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *energy_transfer_qe_2_MC = new TH1D("energy_transfer_qe_2_MC","GENIE QE;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *energy_transfer_qe_3_MC = new TH1D("energy_transfer_qe_3_MC","GENIE QE;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *energy_transfer_mec_0_MC = new TH1D("energy_transfer_mec_0_MC","GENIE MEC;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *energy_transfer_mec_1_MC = new TH1D("energy_transfer_mec_1_MC","GENIE MEC;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *energy_transfer_mec_2_MC = new TH1D("energy_transfer_mec_2_MC","GENIE MEC;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *energy_transfer_mec_3_MC = new TH1D("energy_transfer_mec_3_MC","GENIE MEC;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *energy_transfer_res_0_MC = new TH1D("energy_transfer_res_0_MC","GENIE RES;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *energy_transfer_res_1_MC = new TH1D("energy_transfer_res_1_MC","GENIE RES;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *energy_transfer_res_2_MC = new TH1D("energy_transfer_res_2_MC","GENIE RES;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *energy_transfer_res_3_MC = new TH1D("energy_transfer_res_3_MC","GENIE RES;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *energy_transfer_dis_0_MC = new TH1D("energy_transfer_dis_0_MC","GENIE DIS;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *energy_transfer_dis_1_MC = new TH1D("energy_transfer_dis_1_MC","GENIE DIS;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *energy_transfer_dis_2_MC = new TH1D("energy_transfer_dis_2_MC","GENIE DIS;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *energy_transfer_dis_3_MC = new TH1D("energy_transfer_dis_3_MC","GENIE DIS;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *invariant_mass_MC = new TH1D("invariant_mass_MC","GENIE;Mass (GeV/c^{2});Counts",100,0,2);
  TH1D *el_vz_h_MC = new TH1D("el_vz_h_MC","GENIE;Counts",1000,-10,10);
  TH1D *W_MC = new TH1D("W_MC","GENIE;W (GeV/c^{2};Counts",400,0.5,2.2);
  TH1D *W_qe_MC = new TH1D("W_qe_MC","GENIE QE;W (GeV/c^{2};Counts",400,0.5,2.2);
  TH1D *W_mec_MC = new TH1D("W_mec_MC","GENIE MEC;W (GeV/c^{2};Counts",400,0.5,2.2);
  TH1D *W_res_MC = new TH1D("W_res_MC","GENIE RES;W (GeV/c^{2};Counts",400,0.5,2.2);
  TH1D *W_resonances_1p2_to_2p1_MC = new TH1D("W_resonances_1p2_to_2p1_MC","GENIE Resonance Events with W #epsilon (1.2,2.1) GeV/c^{2};W (GeV/c^{2};Counts",400,0.5,2.2);
  TH1D *W_dis_MC = new TH1D("W_dis_MC","GENIE DIS;W (GeV/c^{2}",400,0.5,2.2);
  TH1D *theta_MC = new TH1D("theta_MC","GENIE;#theta;Counts",1000,0.,60.);
  TH1D *theta_qe_MC = new TH1D("theta_qe_MC","GENIE;#theta;Counts",1000,0.,60.);
  TH1D *theta_mec_MC = new TH1D("theta_mec_MC","GENIE MEC;#theta;Counts",1000,0.,60.);
  TH1D *theta_res_MC = new TH1D("theta_res_MC","GENIE RES;#theta;Counts",1000,0.,60.);
  TH1D *theta_dis_MC = new TH1D("theta_dis_MC","GENIE DIS;#theta;Counts",1000,0.,60.);
  TH1D *omega_MC = new TH1D("omega_MC","GENIE;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *omega_qe_MC = new TH1D("omega_qe_MC","GENIE QE;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *omega_mec_MC = new TH1D("omega_mec_MC","GENIE MEC;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *omega_res_MC = new TH1D("omega_RES_MC","GENIE RES;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *omega_dis_MC = new TH1D("omega_dis_MC","GENIE DIS;#omega (GeV);Counts",1000,0.,1.8);
  TH1D *p_proton_MC = new TH1D("p_proton_MC","GENIE;Momentum (GeV/c);Counts",1000,0.,4.);
  TH1D *p_proton_qe_MC = new TH1D("p_proton_qe_MC","GENIE QE;Momentum (GeV/c);Counts",1000,0.,4.);
  TH1D *p_proton_mec_MC = new TH1D("p_proton_mec_MC","GENIE MEC;Momentum (GeV/c);Counts",1000,0.,4.);
  TH1D *p_proton_res_MC = new TH1D("p_proton_res_MC","GENIE RES;Momentum (GeV/c);Counts",1000,0.,4.);
  TH1D *p_proton_dis_MC = new TH1D("p_proton_dis_MC","GENIE DIS;Momentum (GeV/c);Counts",1000,0.,4.);
  TH1D *p_miss_MC = new TH1D("p_miss_MC","GENIE;Missing Momentum (GeV/c);Counts",1000,-0.5,6.);
  TH1D *p_miss_theta_MC = new TH1D("p_miss_theta_MC","GENIE;#theta;Counts",1000,0.,200.);
  TH1D *processid_MC = new TH1D("processid_MC", "GENIE: Process ID Number;Process ID: 1 = QE, 2 = MEC, 3 = RES, 4 = DIS;Counts",4,0.5,4.5);
  TH1D *resid_MC = new TH1D("resid_MC","GENIE Resonance ID Numbers for Resonance Events;Resonance ID;Counts",21,0.,20.);
  TH1D *resonances_1p2_to_2p1_MC = new TH1D("resonances_1p2_to_2p1_MC","GENIE ID Numbers for Resonance Events with W #epsilon (1.2,2.1) GeV/c^{2};Resonance ID;Counts",21,0.,20.);
  TH1D *xb_MC = new TH1D("xb_MC","GENIE;x_{B};Counts",1000,0.,3.);
  TH1D *xb_qe_MC = new TH1D("xb_qe_MC","GENIE QE;x_{B};Counts",1000,0.,3.);
  TH1D *xb_mec_MC = new TH1D("xb_mec_MC","GENIE MEC;x_{B};Counts",1000,0.,3.);
  TH1D *xb_res_MC = new TH1D("xb_res_MC","GENIE RES;x_{B};Counts",1000,0.,3.);
  TH1D *xb_dis_MC = new TH1D("xb_dis_MC","GENIE DIS;x_{B};Counts",1000,0.,3.);
  TH1D *y_MC = new TH1D("y_MC","GENIE;y = -q + #sqrt{#omega^{2} + 2 #omega M};Counts",1000,-0.5,1.);
  TH1D *y_qe_MC = new TH1D("y_qe_MC","GENIE QE;y = -q + #sqrt{#omega^{2} + 2 #omega M};Counts",1000,-0.5,1.);
  TH1D *y_mec_MC = new TH1D("y_mec_MC","GENIE MEC;y = -q + #sqrt{#omega^{2} + 2 #omega M};Counts",1000,-0.5,1.);
  TH1D *y_res_MC = new TH1D("y_res_MC","GENIE RES;y = -q + #sqrt{#omega^{2} + 2 #omega M};Counts",1000,-0.5,1.);
  TH1D *y_dis_MC = new TH1D("y_dis_MC","GENIE DIS;y = -q + #sqrt{#omega^{2} + 2 #omega M};Counts",1000,-0.5,1.);

  //2D histograms for specific data/MC samples
  //Data
  TH2D *q2_theta_data = new TH2D("q2_theta_data","Data: Electron Angle vs. Q^{2};Q^{2};#theta",400,0.,1.2,400,0.,20.);
  TH2D *ecal_data = new TH2D("ecal_data","Data: Electron Momentum Sampling Fraction;Electron Momentum (GeV/c);Sampling Fraction",400,0.,1.2,1000,0.,0.5);
  TH2D *theta_mom_data = new TH2D("theta_mom_data","Data: Electron Angle vs. Momentum;Momentum (GeV/c);#theta",400,0.,3.,400,0.,20.);
  TH2D *q2_W_data = new TH2D("q2_W_data","Data: Q^{2} vs. W;W (GeV/c^{2});Q^{2}",400,0.,2.,400,0.,1.2);
  TH2D *W_theta_data = new TH2D("W_theta_data","Data: W vs #theta;W (GeV/c^{2});#theta",400,0.,2.,400,0.,20.);
  TH2D *omega_theta_data = new TH2D("omega_theta_data","Data: Energy Transfer vs. #theta;#omega (GeV);#theta",1000,0.,1.9,400,0.,20.);
  TH2D *mult_p_pis_data = new TH2D("mult_p_pis_data","Data: Multiplicity of Protons vs. Pions #pi^{#pm};Proton Multiplicity;#pi^{#pm} Multiplicity",5,-0.5,4.5,5,-0.5,4.5);
  TH2D *mult_p_n_data = new TH2D("mult_p_n_data","Data: Multiplicity of Protons vs. Neutrons;Proton Multiplicity;Neutron Multiplicity",5,-0.5,4.5,5,-0.5,4.5);
  TH2D *p_miss_theta_vs_omega_data = new TH2D("p_miss_theta_vs_omega_data","Data: Missing Momentum Angle vs. Energy Transfer;#omega (GeV); #theta",1000,0.,1.8,400,0.,190.);
  
  //MC Simulation
  TH2D *q2_theta_MC = new TH2D("q2_theta_MC","GENIE: Electron Angle vs. Q^{2};Q^{2};#theta",400,0.,1.2,400,0.,20.);
  TH2D *ecal_MC = new TH2D("ecal_MC","GENIE: Electron Momentum Sampling Fraction;Electron Momentum (GeV/c);Sampling Fraction",400,0.,1.2,1000,0.,0.5);
  TH2D *theta_mom_MC = new TH2D("theta_mom_MC","GENIE: Electron Angle vs. Momentum;Momentum (GeV/c);#theta",400,0.,3.,400,0.,20.);
  TH2D *q2_W_MC = new TH2D("q2_W_MC","GENIE: Q^{2} vs. W;W (GeV/c^{2});Q^{2}",400,0.,2.,400,0.,1.2);
  TH2D *W_theta_MC = new TH2D("W_theta_MC","GENIE: W vs #theta;W (GeV/c^{2});#theta",400,0.,2.,400,0.,20.);
  TH2D *omega_theta_MC = new TH2D("omega_theta_MC","GENIE: Energy Transfer vs. #theta;#omega (GeV);#theta",1000,0.,1.9,400,0.,20.);
  TH2D *mult_p_pis_MC = new TH2D("mult_p_pis_MC","GENIE: Multiplicity of Protons vs. Pions #pi^{#pm};Proton Multiplicity;#pi^{#pm} Multiplicity",5,-0.5,4.5,5,-0.5,4.5);
  TH2D *mult_p_n_MC = new TH2D("mult_p_n_MC","GENIE: Multiplicity of Protons vs. Neutrons;Proton Multiplicity;Neutron Multiplicity",5,-0.5,4.5,5,-0.5,4.5);
  TH2D *p_miss_theta_vs_omega_MC = new TH2D("p_miss_theta_vs_omega_MC","GENIE: Missing Momentum Angle vs. Energy Transfer;#omega (GeV); #theta",1000,0.,1.8,400,0.,190.);
  TH2D *W_resid_1p2_to_2p1_MC = new TH2D("W_resid","GENIE Resonance ID Numbers vs. W #epsilon (1.2,2.1) GeV/c^{2};W (GeV/c^{2});Resonance ID",400,1.15,2.15,21,0.,20.);
  
  //omega vs q3 plots
  TH2D *q0_vs_q3_data = new TH2D("q0_vs_q3_data","Data: q_{0} vs. |q_{3}|;|q_{3}| (GeV/c);#omega (GeV)",1000,0.,1.8,1000,0.,1.3);
  TH2D *q0_vs_q3_MC = new TH2D("q0_vs_q3_MC","GENIE: q_{0} vs. |q_{3}|;|q_{3}| (GeV/c);#omega (GeV)",1000,0.,1.8,1000,0.,1.3);
  TH2D *q0_vs_q3_1pXn0pi_data = new TH2D("q0_vs_q3_1pXn0pi_data","Data, 1pXn0pi: q_{0} vs. |q_{3}|;|q_{3}| (GeV/c);#omega (GeV)",1000,0.,1.8,1000,0.,1.3);
  TH2D *q0_vs_q3_1pXn0pi_MC = new TH2D("q0_vs_q3_1pXn0pi_MC","GENIE, 1pXn0pi: q_{0} vs. |q_{3}|;|q_{3}| (GeV/c);#omega (GeV)",1000,0.,1.8,1000,0.,1.3);
  TH2D *q0_vs_q3_2pXn0pi_data = new TH2D("q0_vs_q3_2pXn0pi_data","Data, 2pXn0pi: q_{0} vs. |q_{3}|;|q_{3}| (GeV/c);#omega (GeV)",1000,0.,1.8,1000,0.,1.3);
  TH2D *q0_vs_q3_2pXn0pi_MC = new TH2D("q0_vs_q3_2pXn0pi_MC","GENIE, 2pXn0pi: q_{0} vs. |q_{3}|;|q_{3}| (GeV/c);#omega (GeV)",1000,0.,1.8,1000,0.,1.3);

  //Primakoff production plots
  TH2D *invmass_vs_theta_0pXn1pip1pim_data = new TH2D("invmass_vs_theta_0pXn1pip1pim_data","Data: 0pXn1#pi^{+}1#pi^{-};#theta; 1#pi^{+}1#pi^{-} Invariant Mass (GeV/c^{2})",400,0.,20.,1000,0.,2.);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //2D Histograms of particular topologies without cuts on various kinematic variables (what should the cuts be?)
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //0p1pip
  //Potential nucleon reconstructed invariant mass vs. kinematic variables
  TH2D *nuc_inv_mass_W_0p_1pip_data = new TH2D("nuc_inv_mass_W_0p_1pip_data","Data: Potential Nucleon Invariant Mass vs. Invariant Mass, W, for 0p1#pi^{+};W (GeV/c^{2});Nucleon Invariant Mass (GeV/c^{2})",400,0.,3.5,400,0.,6.);
  TH2D *nuc_inv_mass_W_0p_1pip_MC = new TH2D("nuc_inv_mass_W_0p_1pip_MC","GENIE: Potential Nucleon Invariant Mass vs. Invariant Mass, W, for 0p1#pi^{+};W (GeV/c^{2});Nucleon Invariant Mass (GeV/c^{2})",400,0.,3.5,400,0.,6.);
  TH2D *nuc_inv_mass_omega_0p_1pip_data = new TH2D("nuc_inv_mass_omega_0p_1pip_data","Data: Potential Nucleon Invariant Mass vs. Energy Transfer for 0p1#pi^{+};#omega (GeV);Nucleon Invariant Mass (GeV/c^{2})",400,0.,3.5,400,0.,6.);
  TH2D *nuc_inv_mass_omega_0p_1pip_MC = new TH2D("nuc_inv_mass_omega_0p_1pip_MC","GENIE: Potential Nucleon Invariant Mass vs. Energy Transfer for 0p1#pi^{+};#omega (GeV);Nucleon Invariant Mass (GeV/c^{2})",400,0.,3.5,400,0.,6.);
  TH2D *nuc_inv_mass_q2_0p_1pip_data = new TH2D("nuc_inv_mass_q2_0p_1pip_data","Data: Potential Nucleon Invariant Mass vs. Q^{2} for 0p1#pi^{+};Q^{2};Nucleon Invariant Mass (GeV/c^{2})",400,0.,6.,400,0.,6.);
  TH2D *nuc_inv_mass_q2_0p_1pip_MC = new TH2D("nuc_inv_mass_q2_0p_1pip_MC","GENIE: Potential Nucleon Invariant Mass vs. Q^{2} for 0p1#pi^{+};Q^{2};Nucleon Invariant Mass (GeV/c^{2})",400,0.,6.,400,0.,6.);
  TH2D *nuc_inv_mass_xb_0p_1pip_data = new TH2D("nuc_inv_mass_xb_0p_1pip_data","Data: Potential Nucleon Invariant Mass vs. x_{B} for 0p1#pi^{+};x_{B};Nucleon Invariant Mass (GeV/c^{2})",400,0.,2.5,400,0.,6.);
  TH2D *nuc_inv_mass_xb_0p_1pip_MC = new TH2D("nuc_inv_mass_xb_0p_1pip_MC","GENIE: Potential Nucleon Invariant Mass vs. x_{B} for 0p1#pi^{+};x_{B};Nucleon Invariant Mass (GeV/c^{2})",400,0.,2.5,400,0.,6.);
  TH2D *nuc_inv_mass_y_0p_1pip_data = new TH2D("nuc_inv_mass_y_0p_1pip_data","Data: Potential Nucleon Invariant Mass vs. y for 0p1#pi^{+};y (GeV);Nucleon Invariant Mass (GeV/c^{2})",400,-0.5,1.,400,0.,6.);
  TH2D *nuc_inv_mass_y_0p_1pip_MC = new TH2D("nuc_inv_mass_y_0p_1pip_MC","GENIE: Potential Nucleon Invariant Mass vs. y for 0p1#pi^{+};y (GeV);Nucleon Invariant Mass (GeV/c^{2})",400,-0.5,1.,400,0.,6.);
  //Total Kinetic Energy vs. kinematic variables
  TH2D *totrecoe_W_0p_1pip_data = new TH2D("totrecoe_W_0p_1pip_data","Data: Total Kinetic Energy vs. Invariant Mass, W, for 0p1#pi^{+};W (GeV/c^{2});Total Kinetic Energy (GeV)",400,0.,3.5,400,0.,6.);
  TH2D *totrecoe_W_0p_1pip_MC = new TH2D("totrecoe_W_0p_1pip_MC","GENIE: Total Kinetic Energy vs. Invariant Mass, W, for 0p1#pi^{+};W (GeV/c^{2});Total Kinetic Energy (GeV)",400,0.,3.5,400,0.,6.);
  TH2D *totrecoe_omega_0p_1pip_data = new TH2D("totrecoe_omega_0p_1pip_data","Data: Total Kinetic Energy vs. Energy Transfer for 0p1#pi^{+};#omega (GeV);Total Kinetic Energy (GeV)",400,0.,3.5,400,0.,6.);
  TH2D *totrecoe_omega_0p_1pip_MC = new TH2D("totrecoe_omega_0p_1pip_MC","GENIE: Total Kinetic Energy vs. Energy Transfer for 0p1#pi^{+};#omega (GeV);Total Kinetic Energy (GeV)",400,0.,3.5,400,0.,6.);
  TH2D *totrecoe_q2_0p_1pip_data = new TH2D("totrecoe_q2_0p_1pip_data","Data: Total Kinetic Energy vs. Q^{2} for 0p1#pi^{+};Q^{2};Total Kinetic Energy (GeV)",400,0.,6.,400,0.,6.);
  TH2D *totrecoe_q2_0p_1pip_MC = new TH2D("totrecoe_q2_0p_1pip_MC","GENIE: Total Kinetic Energy vs. Q^{2} for 0p1#pi^{+};Q^{2};Total Kinetic Energy (GeV)",400,0.,6.,400,0.,6.);
  TH2D *totrecoe_xb_0p_1pip_data = new TH2D("totrecoe_xb_0p_1pip_data","Data: Total Kinetic Energy vs. x_{B} for 0p1#pi^{+};x_{B};Total Kinetic Energy (GeV)",400,0.,2.5,400,0.,6.);
  TH2D *totrecoe_xb_0p_1pip_MC = new TH2D("totrecoe_xb_0p_1pip_MC","GENIE: Total Kinetic Energy vs. x_{B} for 0p1#pi^{+};x_{B};Total Kinetic Energy (GeV)",400,0.,2.5,400,0.,6.);
  TH2D *totrecoe_y_0p_1pip_data = new TH2D("totrecoe_y_0p_1pip_data","Data: Total Kinetic Energy vs. y for 0p1#pi^{+};y (GeV);Total Kinetic Energy (GeV)",400,-0.5,1.,400,0.,6.);
  TH2D *totrecoe_y_0p_1pip_MC = new TH2D("totrecoe_y_0p_1pip_MC","GENIE: Total Kinetic Energy vs. y for 0p1#pi^{+};y (GeV);Total Kinetic Energy (GeV)",400,-0.5,1.,400,0.,6.);
  //Invariant mass vs. W, x, and y
  TH2D *W_q2_0p_1pip_data = new TH2D("W_q2_0p_1pip_data","Data: Invariant Mass, W vs. Q^{2} for 0p1#pi^{+};Q^{2};W (GeV/c^{2})",400,0.,6.,400,0.,3.5);
  TH2D *W_q2_0p_1pip_MC = new TH2D("W_q2_0p_1pip_MC","GENIE: Invariant Mass, W vs. Q^{2} for 0p1#pi^{+};Q^{2};W (GeV/c^{2})",400,0.,6.,400,0.,3.5);
  TH2D *W_xb_0p_1pip_data = new TH2D("W_xb_0p_1pip_data","Data: Invariant Mass, W vs. x_{B} for 0p1#pi^{+};x_{B};W (GeV/c^{2})",400,0.,2.5,400,0.,3.5);
  TH2D *W_xb_0p_1pip_MC = new TH2D("W_xb_0p_1pip_MC","GENIE: Invariant Mass, W vs. x_{B} for 0p1#pi^{+};x_{B};W (GeV/c^{2})",400,0.,2.5,400,0.,3.5);
  TH2D *W_y_0p_1pip_data = new TH2D("W_y_0p_1pip_data","Data: Invariant Mass, W vs. y for 0p1#pi^{+};y (GeV);W (GeV/c^{2})",400,-0.5,1.,400,0.,3.5);
  TH2D *W_y_0p_1pip_MC = new TH2D("W_y_0p_1pip_MC","GENIE: Invariant Mass, W vs. y for 0p1#pi^{+};y (GeV);W (GeV/c^{2})",400,-0.5,1.,400,0.,3.5);
  //Q^2 vs. x and y
  TH2D *q2_xb_0p_1pip_data = new TH2D("q2_xb_0p_1pip_data","Data: Q^{2} vs. x_{B} for 0p1#pi^{+};x_{B};Q^{2}",400,0.,2.5,400,0.,6.);
  TH2D *q2_xb_0p_1pip_MC = new TH2D("q2_xb_0p_1pip_MC","GENIE: Q^{2} vs. x_{B} for 0p1#pi^{+};x_{B};Q^{2}",400,0.,2.5,400,0.,6.);
  TH2D *q2_y_0p_1pip_data = new TH2D("q2_y_0p_1pip_data","Data: Q^{2} vs. y for 0p1#pi^{+};y (GeV);Q^{2}",400,-0.5,1.,400,0.,6.);
  TH2D *q2_y_0p_1pip_MC = new TH2D("q2_y_0p_1pip_MC","GENIE: Q^{2} vs. y for 0p1#pi^{+};y (GeV);Q^{2}",400,-0.5,1.,400,0.,6.);
  //xb vs. y
  TH2D *y_xb_0p_1pip_data = new TH2D("y_xb_0p_1pip_data","Data: x_{B} vs. y for 0p1#pi^{+};x_{B};y(GeV)",400,0.,2.5,400,-0.5,1.);
  TH2D *y_xb_0p_1pip_MC = new TH2D("y_xb_0p_1pip_MC","GENIE: x_{B} vs. y for 0p1#pi^{+};x_{B};y (GeV)",400,0.,2.5,400,-0.5,1.);
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //0p1pim
  //Potential nucleon reconstructed invariant mass vs. kinematic variables
  TH2D *nuc_inv_mass_W_0p_1pim_data = new TH2D("nuc_inv_mass_W_0p_1pim_data","Data: Potential Nucleon Invariant Mass vs. Invariant Mass, W, for 0p1#pi^{-};W (GeV/c^{2});Nucleon Invariant Mass (GeV/c^{2})",400,0.,3.5,400,0.,6.);
  TH2D *nuc_inv_mass_W_0p_1pim_MC = new TH2D("nuc_inv_mass_W_0p_1pim_MC","GENIE: Potential Nucleon Invariant Mass vs. Invariant Mass, W, for 0p1#pi^{-};W (GeV/c^{2});Nucleon Invariant Mass (GeV/c^{2})",400,0.,3.5,400,0.,6.);
  TH2D *nuc_inv_mass_omega_0p_1pim_data = new TH2D("nuc_inv_mass_omega_0p_1pim_data","Data: Potential Nucleon Invariant Mass vs. Energy Transfer for 0p1#pi^{-};#omega (GeV);Nucleon Invariant Mass (GeV/c^{2})",400,0.,3.5,400,0.,6.);
  TH2D *nuc_inv_mass_omega_0p_1pim_MC = new TH2D("nuc_inv_mass_omega_0p_1pim_MC","GENIE: Potential Nucleon Invariant Mass vs. Energy Transfer for 0p1#pi^{-};#omega (GeV);Nucleon Invariant Mass (GeV/c^{2})",400,0.,3.5,400,0.,6.);
  TH2D *nuc_inv_mass_q2_0p_1pim_data = new TH2D("nuc_inv_mass_q2_0p_1pim_data","Data: Potential Nucleon Invariant Mass vs. Q^{2} for 0p1#pi^{-};Q^{2};Nucleon Invariant Mass (GeV/c^{2})",400,0.,6.,400,0.,6.);
  TH2D *nuc_inv_mass_q2_0p_1pim_MC = new TH2D("nuc_inv_mass_q2_0p_1pim_MC","GENIE: Potential Nucleon Invariant Mass vs. Q^{2} for 0p1#pi^{-};Q^{2};Nucleon Invariant Mass (GeV/c^{2})",400,0.,6.,400,0.,6.);
  TH2D *nuc_inv_mass_xb_0p_1pim_data = new TH2D("nuc_inv_mass_xb_0p_1pim_data","Data: Potential Nucleon Invariant Mass vs. x_{B} for 0p1#pi^{-};x_{B};Nucleon Invariant Mass (GeV/c^{2})",400,0.,2.5,400,0.,6.);
  TH2D *nuc_inv_mass_xb_0p_1pim_MC = new TH2D("nuc_inv_mass_xb_0p_1pim_MC","GENIE: Potential Nucleon Invariant Mass vs. x_{B} for 0p1#pi^{-};x_{B};Nucleon Invariant Mass (GeV/c^{2})",400,0.,2.5,400,0.,6.);
  TH2D *nuc_inv_mass_y_0p_1pim_data = new TH2D("nuc_inv_mass_y_0p_1pim_data","Data: Potential Nucleon Invariant Mass vs. y for 0p1#pi^{-};y (GeV);Nucleon Invariant Mass (GeV/c^{2})",400,-0.5,1.,400,0.,6.);
  TH2D *nuc_inv_mass_y_0p_1pim_MC = new TH2D("nuc_inv_mass_y_0p_1pim_MC","GENIE: Potential Nucleon Invariant Mass vs. y for 0p1#pi^{-};y (GeV);Nucleon Invariant Mass (GeV/c^{2})",400,-0.5,1.,400,0.,6.);
  //Total Kinetic Energy vs. kinematic variables
  TH2D *totrecoe_W_0p_1pim_data = new TH2D("totrecoe_W_0p_1pim_data","Data: Total Kinetic Energy vs. Invariant Mass, W, for 0p1#pi^{-};W (GeV/c^{2});Total Kinetic Energy (GeV)",400,0.,3.5,400,0.,6.);
  TH2D *totrecoe_W_0p_1pim_MC = new TH2D("totrecoe_W_0p_1pim_MC","GENIE: Total Kinetic Energy vs. Invariant Mass, W, for 0p1#pi^{-};W (GeV/c^{2});Total Kinetic Energy (GeV)",400,0.,3.5,400,0.,6.);
  TH2D *totrecoe_omega_0p_1pim_data = new TH2D("totrecoe_omega_0p_1pim_data","Data: Total Kinetic Energy vs. Energy Transfer  for 0p1#pi^{-};#omega (GeV);Total Kinetic Energy (GeV)",400,0.,3.5,400,0.,6.);
  TH2D *totrecoe_omega_0p_1pim_MC = new TH2D("totrecoe_omega_0p_1pim_MC","GENIE: Total Kinetic Energy vs. Energy Transfer  for 0p1#pi^{-};#omega (GeV);Total Kinetic Energy (GeV)",400,0.,3.5,400,0.,6.);
  TH2D *totrecoe_q2_0p_1pim_data = new TH2D("totrecoe_q2_0p_1pim_data","Data: Total Kinetic Energy vs. Q^{2} for 0p1#pi^{-};Q^{2};Total Kinetic Energy (GeV)",400,0.,6.,400,0.,6.);
  TH2D *totrecoe_q2_0p_1pim_MC = new TH2D("totrecoe_q2_0p_1pim_MC","GENIE: Total Kinetic Energy vs. Q^{2} for 0p1#pi^{-};Q^{2};Total Kinetic Energy (GeV)",400,0.,6.,400,0.,6.);
  TH2D *totrecoe_xb_0p_1pim_data = new TH2D("totrecoe_xb_0p_1pim_data","Data: Total Kinetic Energy vs. x_{B} for 0p1#pi^{-};x_{B};Total Kinetic Energy (GeV)",400,0.,2.5,400,0.,6.);
  TH2D *totrecoe_xb_0p_1pim_MC = new TH2D("totrecoe_xb_0p_1pim_MC","GENIE: Total Kinetic Energy vs. x_{B} for 0p1#pi^{-};x_{B};Total Kinetic Energy (GeV)",400,0.,2.5,400,0.,6.);
  TH2D *totrecoe_y_0p_1pim_data = new TH2D("totrecoe_y_0p_1pim_data","Data: Total Kinetic Energy vs. y for 0p1#pi^{-};y (GeV);Total Kinetic Energy (GeV)",400,-0.5,1.,400,0.,6.);
  TH2D *totrecoe_y_0p_1pim_MC = new TH2D("totrecoe_y_0p_1pim_MC","GENIE: Total Kinetic Energy vs. y for 0p1#pi^{-};y (GeV);Total Kinetic Energy (GeV)",400,-0.5,1.,400,0.,6.);
  //Invariant mass vs. W, x, and y
  TH2D *W_q2_0p_1pim_data = new TH2D("W_q2_0p_1pim_data","Data: Invariant Mass, W vs. Q^{2} for 0p1#pi^{-};Q^{2};W (GeV/c^{2})",400,0.,6.,400,0.,3.5);
  TH2D *W_q2_0p_1pim_MC = new TH2D("W_q2_0p_1pim_MC","GENIE: Invariant Mass, W vs. Q^{2} for 0p1#pi^{-};Q^{2};W (GeV/c^{2})",400,0.,6.,400,0.,3.5);
  TH2D *W_xb_0p_1pim_data = new TH2D("W_xb_0p_1pim_data","Data: Invariant Mass, W vs. x_{B} for 0p1#pi^{-};x_{B};W (GeV/c^{2})",400,0.,2.5,400,0.,3.5);
  TH2D *W_xb_0p_1pim_MC = new TH2D("W_xb_0p_1pim_MC","GENIE: Invariant Mass, W vs. x_{B} for 0p1#pi^{-};x_{B};W (GeV/c^{2})",400,0.,2.5,400,0.,3.5);
  TH2D *W_y_0p_1pim_data = new TH2D("W_y_0p_1pim_data","Data: Invariant Mass, W vs. y for 0p1#pi^{-};y (GeV);W (GeV/c^{2})",400,-0.5,1.,400,0.,3.5);
  TH2D *W_y_0p_1pim_MC = new TH2D("W_y_0p_1pim_MC","GENIE: Invariant Mass, W vs. y for 0p1#pi^{-};y (GeV);W (GeV/c^{2})",400,-0.5,1.,400,0.,3.5);
  //Q^2 vs. x and y
  TH2D *q2_xb_0p_1pim_data = new TH2D("q2_xb_0p_1pim_data","Data: Q^{2} vs. x_{B} for 0p1#pi^{-};x_{B};Q^{2}",400,0.,2.5,400,0.,6.);
  TH2D *q2_xb_0p_1pim_MC = new TH2D("q2_xb_0p_1pim_MC","GENIE: Q^{2} vs. x_{B} for 0p1#pi^{-};x_{B};Q^{2}",400,0.,2.5,400,0.,6.);
  TH2D *q2_y_0p_1pim_data = new TH2D("q2_y_0p_1pim_data","Data: Q^{2} vs. y for 0p1#pi^{-};y (GeV);Q^{2}",400,-0.5,1.,400,0.,6.);
  TH2D *q2_y_0p_1pim_MC = new TH2D("q2_y_0p_1pim_MC","GENIE: Q^{2} vs. y for 0p1#pi^{-};y (GeV);Q^{2}",400,-0.5,1.,400,0.,6.);
  //xb vs. y
  TH2D *y_xb_0p_1pim_data = new TH2D("y_xb_0p_1pim_data","Data: x_{B} vs. y for 0p1#pi^{-};x_{B};y(GeV)",400,0.,2.5,400,-0.5,1.);
  TH2D *y_xb_0p_1pim_MC = new TH2D("y_xb_0p_1pim_MC","GENIE: x_{B} vs. y for 0p1#pi^{-};x_{B};y (GeV)",400,0.,2.5,400,-0.5,1.);
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //1p0pi (no pip or pim)
  //Potential nucleon reconstructed invariant mass vs. kinematic variables
  TH2D *nuc_inv_mass_W_1pXn0pi_data = new TH2D("nuc_inv_mass_W_1pXn0pi_data","Data: Potential Nucleon Invariant Mass vs. Invariant Mass, W, for 1p0#pi^{#pm};W (GeV/c^{2});Nucleon Invariant Mass (GeV/c^{2})",400,0.,3.5,400,0.,6.);
  TH2D *nuc_inv_mass_W_1pXn0pi_MC = new TH2D("nuc_inv_mass_W_1pXn0pi_MC","GENIE: Potential Nucleon Invariant Mass vs. Invariant Mass, W, for 1p0#pi^{#pm};W (GeV/c^{2});Nucleon Invariant Mass (GeV/c^{2})",400,0.,3.5,400,0.,6.);
  TH2D *nuc_inv_mass_omega_1pXn0pi_data = new TH2D("nuc_inv_mass_omega_1pXn0pi_data","Data: Potential Nucleon Invariant Mass vs. Energy Transfer for 1p0#pi^{#pm};#omega (GeV);Nucleon Invariant Mass (GeV/c^{2})",400,0.,3.5,400,0.,6.);
  TH2D *nuc_inv_mass_omega_1pXn0pi_MC = new TH2D("nuc_inv_mass_omega_1pXn0pi_MC","GENIE: Potential Nucleon Invariant Mass vs. Energy Transfer for 1p0#pi^{#pm};#omega (GeV);Nucleon Invariant Mass (GeV/c^{2})",400,0.,3.5,400,0.,6.);
  TH2D *nuc_inv_mass_q2_1pXn0pi_data = new TH2D("nuc_inv_mass_q2_1pXn0pi_data","Data: Potential Nucleon Invariant Mass vs. Q^{2} for 1p0#pi^{#pm};Q^{2};Nucleon Invariant Mass (GeV/c^{2})",400,0.,6.,400,0.,6.);
  TH2D *nuc_inv_mass_q2_1pXn0pi_MC = new TH2D("nuc_inv_mass_q2_1pXn0pi_MC","GENIE: Potential Nucleon Invariant Mass vs. Q^{2} for 1p0#pi^{#pm};Q^{2};Nucleon Invariant Mass (GeV/c^{2})",400,0.,6.,400,0.,6.);
  TH2D *nuc_inv_mass_xb_1pXn0pi_data = new TH2D("nuc_inv_mass_xb_1pXn0pi_data","Data: Potential Nucleon Invariant Mass vs. x_{B} for 1p0#pi^{#pm};x_{B};Nucleon Invariant Mass (GeV/c^{2})",400,0.,2.5,400,0.,6.);
  TH2D *nuc_inv_mass_xb_1pXn0pi_MC = new TH2D("nuc_inv_mass_xb_1pXn0pi_MC","GENIE: Potential Nucleon Invariant Mass vs. x_{B} for 1p0#pi^{#pm};x_{B};Nucleon Invariant Mass (GeV/c^{2})",400,0.,2.5,400,0.,6.);
  TH2D *nuc_inv_mass_y_1pXn0pi_data = new TH2D("nuc_inv_mass_y_1pXn0pi_data","Data: Potential Nucleon Invariant Mass vs. y for 1p0#pi^{#pm};y (GeV);Nucleon Invariant Mass (GeV/c^{2})",400,-0.5,1.,400,0.,6.);
  TH2D *nuc_inv_mass_y_1pXn0pi_MC = new TH2D("nuc_inv_mass_y_1pXn0pi_MC","GENIE: Potential Nucleon Invariant Mass vs. y for 1p0#pi^{#pm};y (GeV);Nucleon Invariant Mass (GeV/c^{2})",400,-0.5,1.,400,0.,6.);
  //Total Kinetic Energy vs. kinematic variables
  TH2D *totrecoe_W_1pXn0pi_data = new TH2D("totrecoe_W_1pXn0pi_data","Data: Total Kinetic Energy vs. Invariant Mass, W, for 1p0#pi^{#pm};W (GeV/c^{2});Total Kinetic Energy (GeV)",400,0.,3.5,400,0.,6.);
  TH2D *totrecoe_W_1pXn0pi_MC = new TH2D("totrecoe_W_1pXn0pi_MC","GENIE: Total Kinetic Energy vs. Invariant Mass, W, for 1p0#pi^{#pm};W (GeV/c^{2});Total Kinetic Energy (GeV)",400,0.,3.5,400,0.,6.);
  TH2D *totrecoe_omega_1pXn0pi_data = new TH2D("totrecoe_omega_1pXn0pi_data","Data: Total Kinetic Energy vs. Energy Transfer for 1p0#pi^{#pm};#omega (GeV);Total Kinetic Energy (GeV)",400,0.,3.5,400,0.,6.);
  TH2D *totrecoe_omega_1pXn0pi_MC = new TH2D("totrecoe_omega_1pXn0pi_MC","GENIE: Total Kinetic Energy vs. Energy Transfer for 1p0#pi^{#pm};#omega (GeV);Total Kinetic Energy (GeV)",400,0.,3.5,400,0.,6.);
  TH2D *totrecoe_q2_1pXn0pi_data = new TH2D("totrecoe_q2_1pXn0pi_data","Data: Total Kinetic Energy vs. Q^{2} for 1p0#pi^{#pm};Q^{2};Total Kinetic Energy (GeV)",400,0.,6.,400,0.,6.);
  TH2D *totrecoe_q2_1pXn0pi_MC = new TH2D("totrecoe_q2_1pXn0pi_MC","GENIE: Total Kinetic Energy vs. Q^{2} for 1p0#pi^{#pm};Q^{2};Total Kinetic Energy (GeV)",400,0.,6.,400,0.,6.);
  TH2D *totrecoe_xb_1pXn0pi_data = new TH2D("totrecoe_xb_1pXn0pi_data","Data: Total Kinetic Energy vs. x_{B} for 1p0#pi^{#pm};x_{B};Total Kinetic Energy (GeV)",400,0.,2.5,400,0.,6.);
  TH2D *totrecoe_xb_1pXn0pi_MC = new TH2D("totrecoe_xb_1pXn0pi_MC","GENIE: Total Kinetic Energy vs. x_{B} for 1p0#pi^{#pm};x_{B};Total Kinetic Energy (GeV)",400,0.,2.5,400,0.,6.);
  TH2D *totrecoe_y_1pXn0pi_data = new TH2D("totrecoe_y_1pXn0pi_data","Data: Total Kinetic Energy vs. y for 1p0#pi^{#pm};y (GeV);Total Kinetic Energy (GeV)",400,-0.5,1.,400,0.,6.);
  TH2D *totrecoe_y_1pXn0pi_MC = new TH2D("totrecoe_y_1pXn0pi_MC","GENIE: Total Kinetic Energy vs. y for 1p0#pi^{#pm};y (GeV);Total Kinetic Energy (GeV)",400,-0.5,1.,400,0.,6.);
  //Invariant mass vs. W, x, and y
  TH2D *W_q2_1pXn0pi_data = new TH2D("W_q2_1pXn0pi_data","Data: Invariant Mass, W vs. Q^{2} for 1p0#pi^{#pm};Q^{2};W (GeV/c^{2})",400,0.,6.,400,0.,3.5);
  TH2D *W_q2_1pXn0pi_MC = new TH2D("W_q2_1pXn0pi_MC","GENIE: Invariant Mass, W vs. Q^{2} for 1p0#pi^{#pm};Q^{2};W (GeV/c^{2})",400,0.,6.,400,0.,3.5);
  TH2D *W_xb_1pXn0pi_data = new TH2D("W_xb_1pXn0pi_data","Data: Invariant Mass, W vs. x_{B} for 1p0#pi^{#pm};x_{B};W (GeV/c^{2})",400,0.,2.5,400,0.,3.5);
  TH2D *W_xb_1pXn0pi_MC = new TH2D("W_xb_1pXn0pi_MC","GENIE: Invariant Mass, W vs. x_{B} for 1p0#pi^{#pm};x_{B};W (GeV/c^{2})",400,0.,2.5,400,0.,3.5);
  TH2D *W_y_1pXn0pi_data = new TH2D("W_y_1pXn0pi_data","Data: Invariant Mass, W vs. y for 1p0#pi^{#pm};y (GeV);W (GeV/c^{2})",400,-0.5,1.,400,0.,3.5);
  TH2D *W_y_1pXn0pi_MC = new TH2D("W_y_1pXn0pi_MC","GENIE: Invariant Mass, W vs. y for 1p0#pi^{#pm};y (GeV);W (GeV/c^{2})",400,-0.5,1.,400,0.,3.5);
  //Q^2 vs. x and y
  TH2D *q2_xb_1pXn0pi_data = new TH2D("q2_xb_1pXn0pi_data","Data: Q^{2} vs. x_{B} for 1p0#pi^{#pm};x_{B};Q^{2}",400,0.,2.5,400,0.,6.);
  TH2D *q2_xb_1pXn0pi_MC = new TH2D("q2_xb_1pXn0pi_MC","GENIE: Q^{2} vs. x_{B} for 1p0#pi^{#pm};x_{B};Q^{2}",400,0.,2.5,400,0.,6.);
  TH2D *q2_y_1pXn0pi_data = new TH2D("q2_y_1pXn0pi_data","Data: Q^{2} vs. y for 1p0#pi^{#pm};y (GeV);Q^{2}",400,-0.5,1.,400,0.,6.);
  TH2D *q2_y_1pXn0pi_MC = new TH2D("q2_y_1pXn0pi_MC","GENIE: Q^{2} vs. y for 1p0#pi^{#pm};y (GeV);Q^{2}",400,-0.5,1.,400,0.,6.);
  //xb vs. y
  TH2D *y_xb_1pXn0pi_data = new TH2D("y_xb_1pXn0pi_data","Data: x_{B} vs. y for 1p0#pi^{#pm};x_{B};y(GeV)",400,0.,2.5,400,-0.5,1.);
  TH2D *y_xb_1pXn0pi_MC = new TH2D("y_xb_1pXn0pi_MC","GENIE: x_{B} vs. y for 1p0#pi^{#pm};x_{B};y (GeV)",400,0.,2.5,400,-0.5,1.);
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //1n0pi (no pip or pim)
  //Potential nucleon reconstructed invariant mass vs. kinematic variables
  TH2D *nuc_inv_mass_W_1n_0pi_data = new TH2D("nuc_inv_mass_W_1n_0pi_data","Data: Potential Nucleon Invariant Mass vs. Invariant Mass, W, for 1n0#pi^{#pm};W (GeV/c^{2});Nucleon Invariant Mass (GeV/c^{2})",400,0.,3.5,400,0.,6.);
  TH2D *nuc_inv_mass_W_1n_0pi_MC = new TH2D("nuc_inv_mass_W_1n_0pi_MC","GENIE: Potential Nucleon Invariant Mass vs. Invariant Mass, W, for 1n0#pi^{#pm};W (GeV/c^{2});Nucleon Invariant Mass (GeV/c^{2})",400,0.,3.5,400,0.,6.);
  TH2D *nuc_inv_mass_omega_1n_0pi_data = new TH2D("nuc_inv_mass_omega_1n_0pi_data","Data: Potential Nucleon Invariant Mass vs. Energy Transfer for 1n0#pi^{#pm};#omega (GeV);Nucleon Invariant Mass (GeV/c^{2})",400,0.,3.5,400,0.,6.);
  TH2D *nuc_inv_mass_omega_1n_0pi_MC = new TH2D("nuc_inv_mass_omega_1n_0pi_MC","GENIE: Potential Nucleon Invariant Mass vs. Energy Transfer for 1n0#pi^{#pm};#omega (GeV);Nucleon Invariant Mass (GeV/c^{2})",400,0.,3.5,400,0.,6.);
  TH2D *nuc_inv_mass_q2_1n_0pi_data = new TH2D("nuc_inv_mass_q2_1n_0pi_data","Data: Potential Nucleon Invariant Mass vs. Q^{2} for 1n0#pi^{#pm};Q^{2};Nucleon Invariant Mass (GeV/c^{2})",400,0.,6.,400,0.,6.);
  TH2D *nuc_inv_mass_q2_1n_0pi_MC = new TH2D("nuc_inv_mass_q2_1n_0pi_MC","GENIE: Potential Nucleon Invariant Mass vs. Q^{2} for 1n0#pi^{#pm};Q^{2};Nucleon Invariant Mass (GeV/c^{2})",400,0.,6.,400,0.,6.);
  TH2D *nuc_inv_mass_xb_1n_0pi_data = new TH2D("nuc_inv_mass_xb_1n_0pi_data","Data: Potential Nucleon Invariant Mass vs. x_{B} for 1n0#pi^{#pm};x_{B};Nucleon Invariant Mass (GeV/c^{2})",400,0.,2.5,400,0.,6.);
  TH2D *nuc_inv_mass_xb_1n_0pi_MC = new TH2D("nuc_inv_mass_xb_1n_0pi_MC","GENIE: Potential Nucleon Invariant Mass vs. x_{B} for 1n0#pi^{#pm};x_{B};Nucleon Invariant Mass (GeV/c^{2})",400,0.,2.5,400,0.,6.);
  TH2D *nuc_inv_mass_y_1n_0pi_data = new TH2D("nuc_inv_mass_y_1n_0pi_data","Data: Potential Nucleon Invariant Mass vs. y for 1n0#pi^{#pm};y (GeV);Nucleon Invariant Mass (GeV/c^{2})",400,-0.5,1.,400,0.,6.);
  TH2D *nuc_inv_mass_y_1n_0pi_MC = new TH2D("nuc_inv_mass_y_1n_0pi_MC","GENIE: Potential Nucleon Invariant Mass vs. y for 1n0#pi^{#pm};y (GeV);Nucleon Invariant Mass (GeV/c^{2})",400,-0.5,1.,400,0.,6.);
  //Total Kinetic Energy vs. kinematic variables
  TH2D *totrecoe_W_1n_0pi_data = new TH2D("totrecoe_W_1n_0pi_data","Data: Total Kinetic Energy vs. Invariant Mass, W, for 1n0#pi^{#pm};W (GeV/c^{2});Total Kinetic Energy (GeV)",400,0.,3.5,400,0.,6.);
  TH2D *totrecoe_W_1n_0pi_MC = new TH2D("totrecoe_W_1n_0pi_MC","GENIE: Total Kinetic Energy vs. Invariant Mass, W, for 1n0#pi^{#pm};W (GeV/c^{2});Total Kinetic Energy (GeV)",400,0.,3.5,400,0.,6.);
  TH2D *totrecoe_omega_1n_0pi_data = new TH2D("totrecoe_omega_1n_0pi_data","Data: Total Kinetic Energy vs. Energy Transfer for 1n0#pi^{#pm};#omega (GeV);Total Kinetic Energy (GeV)",400,0.,3.5,400,0.,6.);
  TH2D *totrecoe_omega_1n_0pi_MC = new TH2D("totrecoe_omega_1n_0pi_MC","GENIE: Total Kinetic Energy vs. Energy Transfer for 1n0#pi^{#pm};#omega (GeV);Total Kinetic Energy (GeV)",400,0.,3.5,400,0.,6.);
  TH2D *totrecoe_q2_1n_0pi_data = new TH2D("totrecoe_q2_1n_0pi_data","Data: Total Kinetic Energy vs. Q^{2} for 1n0#pi^{#pm};Q^{2};Total Kinetic Energy (GeV)",400,0.,6.,400,0.,6.);
  TH2D *totrecoe_q2_1n_0pi_MC = new TH2D("totrecoe_q2_1n_0pi_MC","GENIE: Total Kinetic Energy vs. Q^{2} for 1n0#pi^{#pm};Q^{2};Total Kinetic Energy (GeV)",400,0.,6.,400,0.,6.);
  TH2D *totrecoe_xb_1n_0pi_data = new TH2D("totrecoe_xb_1n_0pi_data","Data: Total Kinetic Energy vs. x_{B} for 1n0#pi^{#pm};x_{B};Total Kinetic Energy (GeV)",400,0.,2.5,400,0.,6.);
  TH2D *totrecoe_xb_1n_0pi_MC = new TH2D("totrecoe_xb_1n_0pi_MC","GENIE: Total Kinetic Energy vs. x_{B} for 1n0#pi^{#pm};x_{B};Total Kinetic Energy (GeV)",400,0.,2.5,400,0.,6.);
  TH2D *totrecoe_y_1n_0pi_data = new TH2D("totrecoe_y_1n_0pi_data","Data: Total Kinetic Energy vs. y for 1n0#pi^{#pm};y (GeV);Total Kinetic Energy (GeV)",400,-0.5,1.,400,0.,6.);
  TH2D *totrecoe_y_1n_0pi_MC = new TH2D("totrecoe_y_1n_0pi_MC","GENIE: Total Kinetic Energy vs. y for 1n0#pi^{#pm};y (GeV);Total Kinetic Energy (GeV)",400,-0.5,1.,400,0.,6.);
  //Invariant mass vs. W, x, and y
  TH2D *W_q2_1n_0pi_data = new TH2D("W_q2_1n_0pi_data","Data: Invariant Mass, W vs. Q^{2} for 1n0#pi^{#pm};Q^{2};W (GeV/c^{2})",400,0.,6.,400,0.,3.5);
  TH2D *W_q2_1n_0pi_MC = new TH2D("W_q2_1n_0pi_MC","GENIE: Invariant Mass, W vs. Q^{2} for 1n0#pi^{#pm};Q^{2};W (GeV/c^{2})",400,0.,6.,400,0.,3.5);
  TH2D *W_xb_1n_0pi_data = new TH2D("W_xb_1n_0pi_data","Data: Invariant Mass, W vs. x_{B} for 1n0#pi^{#pm};x_{B};W (GeV/c^{2})",400,0.,2.5,400,0.,3.5);
  TH2D *W_xb_1n_0pi_MC = new TH2D("W_xb_1n_0pi_MC","GENIE: Invariant Mass, W vs. x_{B} for 1n0#pi^{#pm};x_{B};W (GeV/c^{2})",400,0.,2.5,400,0.,3.5);
  TH2D *W_y_1n_0pi_data = new TH2D("W_y_1n_0pi_data","Data: Invariant Mass, W vs. y for 1n0#pi^{#pm};y (GeV);W (GeV/c^{2})",400,-0.5,1.,400,0.,3.5);
  TH2D *W_y_1n_0pi_MC = new TH2D("W_y_1n_0pi_MC","GENIE: Invariant Mass, W vs. y for 1n0#pi^{#pm};y (GeV);W (GeV/c^{2})",400,-0.5,1.,400,0.,3.5);
  //Q^2 vs. x and y
  TH2D *q2_xb_1n_0pi_data = new TH2D("q2_xb_1n_0pi_data","Data: Q^{2} vs. x_{B} for 1n0#pi^{#pm};x_{B};Q^{2}",400,0.,2.5,400,0.,6.);
  TH2D *q2_xb_1n_0pi_MC = new TH2D("q2_xb_1n_0pi_MC","GENIE: Q^{2} vs. x_{B} for 1n0#pi^{#pm};x_{B};Q^{2}",400,0.,2.5,400,0.,6.);
  TH2D *q2_y_1n_0pi_data = new TH2D("q2_y_1n_0pi_data","Data: Q^{2} vs. y for 1n0#pi^{#pm};y (GeV);Q^{2}",400,-0.5,1.,400,0.,6.);
  TH2D *q2_y_1n_0pi_MC = new TH2D("q2_y_1n_0pi_MC","GENIE: Q^{2} vs. y for 1n0#pi^{#pm};y (GeV);Q^{2}",400,-0.5,1.,400,0.,6.);
  //xb vs. y
  TH2D *y_xb_1n_0pi_data = new TH2D("y_xb_1n_0pi_data","Data: x_{B} vs. y for 1n0#pi^{#pm};x_{B};y(GeV)",400,0.,2.5,400,-0.5,1.);
  TH2D *y_xb_1n_0pi_MC = new TH2D("y_xb_1n_0pi_MC","GENIE: x_{B} vs. y for 1n0#pi^{#pm};x_{B};y (GeV)",400,0.,2.5,400,-0.5,1.);
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //1p1pim
  //Potential nucleon reconstructed invariant mass vs. kinematic variables
  TH2D *nuc_inv_mass_W_1pXn1pim_data = new TH2D("nuc_inv_mass_W_1pXn1pim_data","Data: Potential Nucleon Invariant Mass vs. Invariant Mass, W, for 1p1#pi^{-};W (GeV/c^{2});Nucleon Invariant Mass (GeV/c^{2})",400,0.,3.5,400,0.,6.);
  TH2D *nuc_inv_mass_W_1pXn1pim_MC = new TH2D("nuc_inv_mass_W_1pXn1pim_MC","GENIE: Potential Nucleon Invariant Mass vs. Invariant Mass, W, for 1p1#pi^{-};W (GeV/c^{2});Nucleon Invariant Mass (GeV/c^{2})",400,0.,3.5,400,0.,6.);
  TH2D *nuc_inv_mass_omega_1pXn1pim_data = new TH2D("nuc_inv_mass_omega_1pXn1pim_data","Data: Potential Nucleon Invariant Mass vs. Energy Transfer for 1p1#pi^{-};#omega (GeV);Nucleon Invariant Mass (GeV/c^{2})",400,0.,3.5,400,0.,6.);
  TH2D *nuc_inv_mass_omega_1pXn1pim_MC = new TH2D("nuc_inv_mass_omega_1pXn1pim_MC","GENIE: Potential Nucleon Invariant Mass vs. Energy Transfer for 1p1#pi^{-};#omega (GeV);Nucleon Invariant Mass (GeV/c^{2})",400,0.,3.5,400,0.,6.);
  TH2D *nuc_inv_mass_q2_1pXn1pim_data = new TH2D("nuc_inv_mass_q2_1pXn1pim_data","Data: Potential Nucleon Invariant Mass vs. Q^{2} for 1p1#pi^{-};Q^{2};Nucleon Invariant Mass (GeV/c^{2})",400,0.,6.,400,0.,6.);
  TH2D *nuc_inv_mass_q2_1pXn1pim_MC = new TH2D("nuc_inv_mass_q2_1pXn1pim_MC","GENIE: Potential Nucleon Invariant Mass vs. Q^{2} for 1p1#pi^{-};Q^{2};Nucleon Invariant Mass (GeV/c^{2})",400,0.,6.,400,0.,6.);
  TH2D *nuc_inv_mass_xb_1pXn1pim_data = new TH2D("nuc_inv_mass_xb_1pXn1pim_data","Data: Potential Nucleon Invariant Mass vs. x_{B} for 1p1#pi^{-};x_{B};Nucleon Invariant Mass (GeV/c^{2})",400,0.,2.5,400,0.,6.);
  TH2D *nuc_inv_mass_xb_1pXn1pim_MC = new TH2D("nuc_inv_mass_xb_1pXn1pim_MC","GENIE: Potential Nucleon Invariant Mass vs. x_{B} for 1p1#pi^{-};x_{B};Nucleon Invariant Mass (GeV/c^{2})",400,0.,2.5,400,0.,6.);
  TH2D *nuc_inv_mass_y_1pXn1pim_data = new TH2D("nuc_inv_mass_y_1pXn1pim_data","Data: Potential Nucleon Invariant Mass vs. y for 1p1#pi^{-};y (GeV);Nucleon Invariant Mass (GeV/c^{2})",400,-0.5,1.,400,0.,6.);
  TH2D *nuc_inv_mass_y_1pXn1pim_MC = new TH2D("nuc_inv_mass_y_1pXn1pim_MC","GENIE: Potential Nucleon Invariant Mass vs. y for 1p1#pi^{-};y (GeV);Nucleon Invariant Mass (GeV/c^{2})",400,-0.5,1.,400,0.,6.);
  //Total Kinetic Energy vs. kinematic variables
  TH2D *totrecoe_W_1pXn1pim_data = new TH2D("totrecoe_W_1pXn1pim_data","Data: Total Kinetic Energy vs. Invariant Mass, W, for 1p1#pi^{-};W (GeV/c^{2});Total Kinetic Energy (GeV)",400,0.,3.5,400,0.,6.);
  TH2D *totrecoe_W_1pXn1pim_MC = new TH2D("totrecoe_W_1pXn1pim_MC","GENIE: Total Kinetic Energy vs. Invariant Mass, W, for 1p1#pi^{-};W (GeV/c^{2});Total Kinetic Energy (GeV)",400,0.,3.5,400,0.,6.);
  TH2D *totrecoe_omega_1pXn1pim_data = new TH2D("totrecoe_omega_1pXn1pim_data","Data: Total Kinetic Energy vs. Energy Transfer for 1p1#pi^{-};#omega (GeV);Total Kinetic Energy (GeV)",400,0.,3.5,400,0.,6.);
  TH2D *totrecoe_omega_1pXn1pim_MC = new TH2D("totrecoe_omega_1pXn1pim_MC","GENIE: Total Kinetic Energy vs. Energy Transfer for 1p1#pi^{-};#omega (GeV);Total Kinetic Energy (GeV)",400,0.,3.5,400,0.,6.);
  TH2D *totrecoe_q2_1pXn1pim_data = new TH2D("totrecoe_q2_1pXn1pim_data","Data: Total Kinetic Energy vs. Q^{2} for 1p1#pi^{-};Q^{2};Total Kinetic Energy (GeV)",400,0.,6.,400,0.,6.);
  TH2D *totrecoe_q2_1pXn1pim_MC = new TH2D("totrecoe_q2_1pXn1pim_MC","GENIE: Total Kinetic Energy vs. Q^{2} for 1p1#pi^{-};Q^{2};Total Kinetic Energy (GeV)",400,0.,6.,400,0.,6.);
  TH2D *totrecoe_xb_1pXn1pim_data = new TH2D("totrecoe_xb_1pXn1pim_data","Data: Total Kinetic Energy vs. x_{B} for 1p1#pi^{-};x_{B};Total Kinetic Energy (GeV)",400,0.,2.5,400,0.,6.);
  TH2D *totrecoe_xb_1pXn1pim_MC = new TH2D("totrecoe_xb_1pXn1pim_MC","GENIE: Total Kinetic Energy vs. x_{B} for 1p1#pi^{-};x_{B};Total Kinetic Energy (GeV)",400,0.,2.5,400,0.,6.);
  TH2D *totrecoe_y_1pXn1pim_data = new TH2D("totrecoe_y_1pXn1pim_data","Data: Total Kinetic Energy vs. y for 1p1#pi^{-};y (GeV);Total Kinetic Energy (GeV)",400,-0.5,1.,400,0.,6.);
  TH2D *totrecoe_y_1pXn1pim_MC = new TH2D("totrecoe_y_1pXn1pim_MC","GENIE: Total Kinetic Energy vs. y for 1p1#pi^{-};y (GeV);Total Kinetic Energy (GeV)",400,-0.5,1.,400,0.,6.);
  //Invariant mass vs. W, x, and y
  TH2D *W_q2_1pXn1pim_data = new TH2D("W_q2_1pXn1pim_data","Data: Invariant Mass, W vs. Q^{2} for 1p1#pi^{-};Q^{2};W (GeV/c^{2})",400,0.,6.,400,0.,3.5);
  TH2D *W_q2_1pXn1pim_MC = new TH2D("W_q2_1pXn1pim_MC","GENIE: Invariant Mass, W vs. Q^{2} for 1p1#pi^{-};Q^{2};W (GeV/c^{2})",400,0.,6.,400,0.,3.5);
  TH2D *W_xb_1pXn1pim_data = new TH2D("W_xb_1pXn1pim_data","Data: Invariant Mass, W vs. x_{B} for 1p1#pi^{-};x_{B};W (GeV/c^{2})",400,0.,2.5,400,0.,3.5);
  TH2D *W_xb_1pXn1pim_MC = new TH2D("W_xb_1pXn1pim_MC","GENIE: Invariant Mass, W vs. x_{B} for 1p1#pi^{-};x_{B};W (GeV/c^{2})",400,0.,2.5,400,0.,3.5);
  TH2D *W_y_1pXn1pim_data = new TH2D("W_y_1pXn1pim_data","Data: Invariant Mass, W vs. y for 1p1#pi^{-};y (GeV);W (GeV/c^{2})",400,-0.5,1.,400,0.,3.5);
  TH2D *W_y_1pXn1pim_MC = new TH2D("W_y_1pXn1pim_MC","GENIE: Invariant Mass, W vs. y for 1p1#pi^{-};y (GeV);W (GeV/c^{2})",400,-0.5,1.,400,0.,3.5);
  //Q^2 vs. x and y
  TH2D *q2_xb_1pXn1pim_data = new TH2D("q2_xb_1pXn1pim_data","Data: Q^{2} vs. x_{B} for 1p1#pi^{-};x_{B};Q^{2}",400,0.,2.5,400,0.,6.);
  TH2D *q2_xb_1pXn1pim_MC = new TH2D("q2_xb_1pXn1pim_MC","GENIE: Q^{2} vs. x_{B} for 1p1#pi^{-};x_{B};Q^{2}",400,0.,2.5,400,0.,6.);
  TH2D *q2_y_1pXn1pim_data = new TH2D("q2_y_1pXn1pim_data","Data: Q^{2} vs. y for 1p1#pi^{-};y (GeV);Q^{2}",400,-0.5,1.,400,0.,6.);
  TH2D *q2_y_1pXn1pim_MC = new TH2D("q2_y_1pXn1pim_MC","GENIE: Q^{2} vs. y for 1p1#pi^{-};y (GeV);Q^{2}",400,-0.5,1.,400,0.,6.);
  //xb vs. y
  TH2D *y_xb_1pXn1pim_data = new TH2D("y_xb_1pXn1pim_data","Data: x_{B} vs. y for 1p1#pi^{-};x_{B};y(GeV)",400,0.,2.5,400,-0.5,1.);
  TH2D *y_xb_1pXn1pim_MC = new TH2D("y_xb_1pXn1pim_MC","GENIE: x_{B} vs. y for 1p1#pi^{-};x_{B};y (GeV)",400,0.,2.5,400,-0.5,1.);
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //1n1pip
  //Potential nucleon reconstructed invariant mass vs. kinematic variables
  TH2D *nuc_inv_mass_W_1n_1pip_data = new TH2D("nuc_inv_mass_W_1n_1pip_data","Data: Potential Nucleon Invariant Mass vs. Invariant Mass, W, for 1n1#pi^{+};W (GeV/c^{2});Nucleon Invariant Mass (GeV/c^{2})",400,0.,3.5,400,0.,6.);
  TH2D *nuc_inv_mass_W_1n_1pip_MC = new TH2D("nuc_inv_mass_W_1n_1pip_MC","GENIE: Potential Nucleon Invariant Mass vs. Invariant Mass, W, for 1n1#pi^{+};W (GeV/c^{2});Nucleon Invariant Mass (GeV/c^{2})",400,0.,3.5,400,0.,6.);
  TH2D *nuc_inv_mass_omega_1n_1pip_data = new TH2D("nuc_inv_mass_omega_1n_1pip_data","Data: Potential Nucleon Invariant Mass vs. Energy Transfer for 1n1#pi^{+};#omega (GeV);Nucleon Invariant Mass (GeV/c^{2})",400,0.,3.5,400,0.,6.);
  TH2D *nuc_inv_mass_omega_1n_1pip_MC = new TH2D("nuc_inv_mass_omega_1n_1pip_MC","GENIE: Potential Nucleon Invariant Mass vs. Energy Transfer for 1n1#pi^{+};#omega (GeV);Nucleon Invariant Mass (GeV/c^{2})",400,0.,3.5,400,0.,6.);
  TH2D *nuc_inv_mass_q2_1n_1pip_data = new TH2D("nuc_inv_mass_q2_1n_1pip_data","Data: Potential Nucleon Invariant Mass vs. Q^{2} for 1n1#pi^{+};Q^{2};Nucleon Invariant Mass (GeV/c^{2})",400,0.,6.,400,0.,6.);
  TH2D *nuc_inv_mass_q2_1n_1pip_MC = new TH2D("nuc_inv_mass_q2_1n_1pip_MC","GENIE: Potential Nucleon Invariant Mass vs. Q^{2} for 1n1#pi^{+};Q^{2};Nucleon Invariant Mass (GeV/c^{2})",400,0.,6.,400,0.,6.);
  TH2D *nuc_inv_mass_xb_1n_1pip_data = new TH2D("nuc_inv_mass_xb_1n_1pip_data","Data: Potential Nucleon Invariant Mass vs. x_{B} for 1n1#pi^{+};x_{B};Nucleon Invariant Mass (GeV/c^{2})",400,0.,2.5,400,0.,6.);
  TH2D *nuc_inv_mass_xb_1n_1pip_MC = new TH2D("nuc_inv_mass_xb_1n_1pip_MC","GENIE: Potential Nucleon Invariant Mass vs. x_{B} for 1n1#pi^{+};x_{B};Nucleon Invariant Mass (GeV/c^{2})",400,0.,2.5,400,0.,6.);
  TH2D *nuc_inv_mass_y_1n_1pip_data = new TH2D("nuc_inv_mass_y_1n_1pip_data","Data: Potential Nucleon Invariant Mass vs. y for 1n1#pi^{+};y (GeV);Nucleon Invariant Mass (GeV/c^{2})",400,-0.5,1.,400,0.,6.);
  TH2D *nuc_inv_mass_y_1n_1pip_MC = new TH2D("nuc_inv_mass_y_1n_1pip_MC","GENIE: Potential Nucleon Invariant Mass vs. y for 1n1#pi^{+};y (GeV);Nucleon Invariant Mass (GeV/c^{2})",400,-0.5,1.,400,0.,6.);
  //Total Kinetic Energy vs. kinematic variables
  TH2D *totrecoe_W_1n_1pip_data = new TH2D("totrecoe_W_1n_1pip_data","Data: Total Kinetic Energy vs. Invariant Mass, W, for 1n1#pi^{+};W (GeV/c^{2});Total Kinetic Energy (GeV)",400,0.,3.5,400,0.,6.);
  TH2D *totrecoe_W_1n_1pip_MC = new TH2D("totrecoe_W_1n_1pip_MC","GENIE: Total Kinetic Energy vs. Invariant Mass, W, for 1n1#pi^{+};W (GeV/c^{2});Total Kinetic Energy (GeV)",400,0.,3.5,400,0.,6.);
  TH2D *totrecoe_omega_1n_1pip_data = new TH2D("totrecoe_omega_1n_1pip_data","Data: Total Kinetic Energy vs. Energy Transfer for 1n1#pi^{+};#omega (GeV);Total Kinetic Energy (GeV)",400,0.,3.5,400,0.,6.);
  TH2D *totrecoe_omega_1n_1pip_MC = new TH2D("totrecoe_omega_1n_1pip_MC","GENIE: Total Kinetic Energy vs. Energy Transfer for 1n1#pi^{+};#omega (GeV);Total Kinetic Energy (GeV)",400,0.,3.5,400,0.,6.);
  TH2D *totrecoe_q2_1n_1pip_data = new TH2D("totrecoe_q2_1n_1pip_data","Data: Total Kinetic Energy vs. Q^{2} for 1n1#pi^{+};Q^{2};Total Kinetic Energy (GeV)",400,0.,6.,400,0.,6.);
  TH2D *totrecoe_q2_1n_1pip_MC = new TH2D("totrecoe_q2_1n_1pip_MC","GENIE: Total Kinetic Energy vs. Q^{2} for 1n1#pi^{+};Q^{2};Total Kinetic Energy (GeV)",400,0.,6.,400,0.,6.);
  TH2D *totrecoe_xb_1n_1pip_data = new TH2D("totrecoe_xb_1n_1pip_data","Data: Total Kinetic Energy vs. x_{B} for 1n1#pi^{+};x_{B};Total Kinetic Energy (GeV)",400,0.,2.5,400,0.,6.);
  TH2D *totrecoe_xb_1n_1pip_MC = new TH2D("totrecoe_xb_1n_1pip_MC","GENIE: Total Kinetic Energy vs. x_{B} for 1n1#pi^{+};x_{B};Total Kinetic Energy (GeV)",400,0.,2.5,400,0.,6.);
  TH2D *totrecoe_y_1n_1pip_data = new TH2D("totrecoe_y_1n_1pip_data","Data: Total Kinetic Energy vs. y for 1n1#pi^{+};y (GeV);Total Kinetic Energy (GeV)",400,-0.5,1.,400,0.,6.);
  TH2D *totrecoe_y_1n_1pip_MC = new TH2D("totrecoe_y_1n_1pip_MC","GENIE: Total Kinetic Energy vs. y for 1n1#pi^{+};y (GeV);Total Kinetic Energy (GeV)",400,-0.5,1.,400,0.,6.);
  //Invariant mass vs. W, x, and y
  TH2D *W_q2_1n_1pip_data = new TH2D("W_q2_1n_1pip_data","Data: Invariant Mass, W vs. Q^{2} for 1n1#pi^{+};Q^{2};W (GeV/c^{2})",400,0.,6.,400,0.,3.5);
  TH2D *W_q2_1n_1pip_MC = new TH2D("W_q2_1n_1pip_MC","GENIE: Invariant Mass, W vs. Q^{2} for 1n1#pi^{+};Q^{2};W (GeV/c^{2})",400,0.,6.,400,0.,3.5);
  TH2D *W_xb_1n_1pip_data = new TH2D("W_xb_1n_1pip_data","Data: Invariant Mass, W vs. x_{B} for 1n1#pi^{+};x_{B};W (GeV/c^{2})",400,0.,2.5,400,0.,3.5);
  TH2D *W_xb_1n_1pip_MC = new TH2D("W_xb_1n_1pip_MC","GENIE: Invariant Mass, W vs. x_{B} for 1n1#pi^{+};x_{B};W (GeV/c^{2})",400,0.,2.5,400,0.,3.5);
  TH2D *W_y_1n_1pip_data = new TH2D("W_y_1n_1pip_data","Data: Invariant Mass, W vs. y for 1n1#pi^{+};y (GeV);W (GeV/c^{2})",400,-0.5,1.,400,0.,3.5);
  TH2D *W_y_1n_1pip_MC = new TH2D("W_y_1n_1pip_MC","GENIE: Invariant Mass, W vs. y for 1n1#pi^{+};y (GeV);W (GeV/c^{2})",400,-0.5,1.,400,0.,3.5);
  //Q^2 vs. x and y
  TH2D *q2_xb_1n_1pip_data = new TH2D("q2_xb_1n_1pip_data","Data: Q^{2} vs. x_{B} for 1n1#pi^{+};x_{B};Q^{2}",400,0.,2.5,400,0.,6.);
  TH2D *q2_xb_1n_1pip_MC = new TH2D("q2_xb_1n_1pip_MC","GENIE: Q^{2} vs. x_{B} for 1n1#pi^{+};x_{B};Q^{2}",400,0.,2.5,400,0.,6.);
  TH2D *q2_y_1n_1pip_data = new TH2D("q2_y_1n_1pip_data","Data: Q^{2} vs. y for 1n1#pi^{+};y (GeV);Q^{2}",400,-0.5,1.,400,0.,6.);
  TH2D *q2_y_1n_1pip_MC = new TH2D("q2_y_1n_1pip_MC","GENIE: Q^{2} vs. y for 1n1#pi^{+};y (GeV);Q^{2}",400,-0.5,1.,400,0.,6.);
  //xb vs. y
  TH2D *y_xb_1n_1pip_data = new TH2D("y_xb_1n_1pip_data","Data: x_{B} vs. y for 1n1#pi^{+};x_{B};y(GeV)",400,0.,2.5,400,-0.5,1.);
  TH2D *y_xb_1n_1pip_MC = new TH2D("y_xb_1n_1pip_MC","GENIE: x_{B} vs. y for 1n1#pi^{+};x_{B};y (GeV)",400,0.,2.5,400,-0.5,1.);
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //1D Histograms with cuts on various kinematic variables
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //QE single proton: only one proton without charged pions, y:(-0.300,0.300), x:(0.8,1.2), W:(0.80,1.050), q2:(0.0,3.0)
  //QE W distributions
  TH1D *W_qe_1p_cuts_data = new TH1D("W_qe_1p_cuts_data","Data w/QE Cuts;W (GeV/c^{2};Counts",200,0,1.5);
  TH1D *W_qe_1p_cuts_MC = new TH1D("W_qe_1p_cuts_MC","GENIE w/QE Cuts;W (GeV/c^{2};Counts",200,0,1.5);
  TH1D *W_qe_1p_truth_MC = new TH1D("W_qe_1p_truth_MC","GENIE w/QE Truth;W (GeV/c^{2};Counts",200,0,1.5);
  THStack *W_qe_1p_cuts_stacked = new THStack("W_qe_1p_cuts_stacked","QE Single Proton Events W Distribution Data to GENIE Comparison;W (GeV/c^{2};Counts");
  //QE Bjorken-x distributions
  TH1D *xb_qe_1p_cuts_data = new TH1D("xb_qe_1p_cuts_data","Data w/QE Cuts",200,0,1.5);
  TH1D *xb_qe_1p_cuts_MC = new TH1D("xb_qe_1p_cuts_MC","GENIE w/QE Cuts",200,0,1.5);
  TH1D *xb_qe_1p_truth_MC = new TH1D("xb_qe_1p_truth_MC","GENIE w/QE Truth",200,0,1.5);
  THStack *xb_qe_1p_cuts_stacked = new THStack("xb_qe_1p_cuts_stacked","QE Single Proton Events Bjorken-x Distribution Data to GENIE Comparison;x_{B};Counts");
  //QE single Proton momentum distribution
  TH1D *p_proton_qe_1p_cuts_data = new TH1D("p_proton_qe_1p_cuts_data","Data w/QE Cuts;Momentum (GeV/c);Counts",200,0.0,2.0);
  TH1D *p_proton_qe_1p_cuts_MC = new TH1D("p_proton_qe_1p_cuts_MC","GENIE w/QE Cuts;Momentum (GeV/c);Counts",200,0.0,2.0);
  TH1D *p_proton_qe_1p_truth_MC = new TH1D("p_proton_qe_1p_truth_MC","GENIE w/QE Truth;Momentum(GeV/c);Counts",200,0.0,2.0);
  THStack *p_proton_qe_1p_cuts_stacked = new THStack("p_proton_qe_1p_cuts_stacked","QE Single Proton Events Proton Momentum Distribution Data to GENIE Comparison;Momentum (GeV/c);Counts");
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //MEC single proton: only one proton without charged pions, y:(0.05,0.300), x:(0.8,1.2), W:(0.80,1.50), q2:(0.0,3.0)
  //MEC W distributions
  TH1D *W_mec_1p_cuts_data = new TH1D("W_mec_1p_cuts_data","Data w/MEC Cuts;W (GeV/c^{2};Counts",200,0,1.5);
  TH1D *W_mec_1p_cuts_MC = new TH1D("W_mec_1p_cuts_MC","GENIE w/MEC Cuts;W (GeV/c^{2};Counts",200,0,1.5);
  TH1D *W_mec_1p_truth_MC = new TH1D("W_mec_1p_truth_MC","GENIE w/MEC Truth;W (GeV/c^{2};Counts",200,0,1.5);
  THStack *W_mec_1p_cuts_stacked = new THStack("W_mec_1p_cuts_stacked","MEC Single Proton Events W Distribution Data to GENIE Comparison;W (GeV/c^{2};Counts");
  //MEC Bjorken-x distributions
  TH1D *xb_mec_1p_cuts_data = new TH1D("xb_mec_1p_cuts_data","Data w/MEC Cuts",200,0,1.5);
  TH1D *xb_mec_1p_cuts_MC = new TH1D("xb_mec_1p_cuts_MC","GENIE w/MEC Cuts",200,0,1.5);
  TH1D *xb_mec_1p_truth_MC = new TH1D("xb_mec_1p_truth_MC","GENIE w/MEC Truth",200,0,1.5);
  THStack *xb_mec_1p_cuts_stacked = new THStack("xb_mec_1p_cuts_stacked","MEC Single Proton Events Bjorken-x Distribution Data to GENIE Comparison;x_{B};Counts");
  //MEC single Proton momentum distribution
  TH1D *p_proton_mec_1p_cuts_data = new TH1D("p_proton_mec_1p_cuts_data","Data w/MEC Cuts;Momentum (GeV/c);Counts",200,0.0,2.0);
  TH1D *p_proton_mec_1p_cuts_MC = new TH1D("p_proton_mec_1p_cuts_MC","GENIE w/MEC Cuts;Momentum (GeV/c);Counts",200,0.0,2.0);
  TH1D *p_proton_mec_1p_truth_MC = new TH1D("p_proton_mec_1p_truth_MC","GENIE w/MEC Truth;Momentum(GeV/c);Counts",200,0.0,2.0);
  THStack *p_proton_mec_1p_cuts_stacked = new THStack("p_proton_mec_1p_cuts_stacked","MEC Single Proton Events Proton Momentum Distribution Data to GENIE Comparison;Momentum (GeV/c);Counts");
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //QE neutron: only one neutron without charged pions, y:(-300,300), x:(0.8,1.2), W:(0.800,1.050), q2:(0.0,3.0)

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //MEC two proton (proton-proton): only two protons without charged pions, y:(0.05,0.300), x:(0.8,1.2), W:(0.80,1.50), q2:(0.0,3.0)
  //MEC W distributions
  TH1D *W_mec_2p_cuts_data = new TH1D("W_mec_2p_cuts_data","Data w/MEC Cuts;W (GeV/c^{2};Counts",200,0,1.5);
  TH1D *W_mec_2p_cuts_MC = new TH1D("W_mec_2p_cuts_MC","GENIE w/MEC Cuts;W (GeV/c^{2};Counts",200,0,1.5);
  TH1D *W_mec_2p_truth_MC = new TH1D("W_mec_2p_truth_MC","GENIE w/MEC Truth;W (GeV/c^{2};Counts",200,0,1.5);
  THStack *W_mec_2p_cuts_stacked = new THStack("W_mec_2p_cuts_stacked","MEC Two Proton Events W Distribution Data to GENIE Comparison;W (GeV/c^{2};Counts");
  //MEC Bjorken-x distributions
  TH1D *xb_mec_2p_cuts_data = new TH1D("xb_mec_2p_cuts_data","Data w/MEC Cuts",200,0,1.5);
  TH1D *xb_mec_2p_cuts_MC = new TH1D("xb_mec_2p_cuts_MC","GENIE w/MEC Cuts",200,0,1.5);
  TH1D *xb_mec_2p_truth_MC = new TH1D("xb_mec_2p_truth_MC","GENIE w/MEC Truth",200,0,1.5);
  THStack *xb_mec_2p_cuts_stacked = new THStack("xb_mec_2p_cuts_stacked","MEC Two Proton Events Bjorken-x Distribution Data to GENIE Comparison;x_{B};Counts");
  //MEC single Proton momentum distribution
  TH1D *p_proton_mec_2p_cuts_data = new TH1D("p_proton_mec_2p_cuts_data","Data w/MEC Cuts;Momentum (GeV/c);Counts",200,0.0,2.0);
  TH1D *p_proton_mec_2p_cuts_MC = new TH1D("p_proton_mec_2p_cuts_MC","GENIE w/MEC Cuts;Momentum (GeV/c);Counts",200,0.0,2.0);
  TH1D *p_proton_mec_2p_truth_MC = new TH1D("p_proton_mec_2p_truth_MC","GENIE w/MEC Truth;Momentum(GeV/c);Counts",200,0.0,2.0);
  THStack *p_proton_mec_2p_cuts_stacked = new THStack("p_proton_mec_2p_cuts_stacked","MEC Two Proton Events Proton Momentum Distribution Data to GENIE Comparison;Momentum (GeV/c);Counts");
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //MEC proton-neutron

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //RES single charged pion
  //RES single pion W distributions
  TH1D *W_res_0p1pi_cuts_data = new TH1D("W_res_0p1pi_cuts_data","Data w/RES Cuts;W (GeV/c^{2};Counts",200,1.0,2.5);
  TH1D *W_res_0p1pi_cuts_MC = new TH1D("W_res_0p1pi_cuts_MC","GENIE w/RES Cuts;W (GeV/c^{2};Counts",200,1.0,2.5);
  TH1D *W_res_0p1pi_truth_MC = new TH1D("W_res_0p1pi_truth_MC","GENIE w/RES Truth;W (GeV/c^{2};Counts",200,1.0,2.5);
  THStack *W_res_0p1pi_cuts_stacked = new THStack("W_res_0p1pi_cuts_stacked","RES Single #pi^{#pm} (0p1#pi^{#pm}) Events W Distribution Data to GENIE Comparison;W (GeV/c^{2};Counts");
  //RES single pion Bjorken-x distributions
  TH1D *xb_res_0p1pi_cuts_data = new TH1D("xb_res_0p1pi_cuts_data","Data w/RES Cuts",200,0,1.5);
  TH1D *xb_res_0p1pi_cuts_MC = new TH1D("xb_res_0p1pi_cuts_MC","GENIE w/RES Cuts",200,0,1.5);
  TH1D *xb_res_0p1pi_truth_MC = new TH1D("xb_res_0p1pi_truth_MC","GENIE w/RES Truth",200,0,1.5);
  THStack *xb_res_0p1pi_cuts_stacked = new THStack("xb_res_0p1pi_cuts_stacked","RES Single #pi^{#pm} (0p1#pi^{#pm}) Events Bjorken-x Distribution Data to GENIE Comparison;x_{B};Counts");
  //RES single pion momentum distribution
  TH1D *p_pi_res_0p1pi_cuts_data = new TH1D("p_pi_res_0p1pi_cuts_data","Data w/RES Cuts;Momentum (GeV/c);Counts",200,0.0,2.0);
  TH1D *p_pi_res_0p1pi_cuts_MC = new TH1D("p_pi_res_0p1pi_cuts_MC","GENIE w/RES Cuts;Momentum (GeV/c);Counts",200,0.0,2.0);
  TH1D *p_pi_res_0p1pi_truth_MC = new TH1D("p_pi_res_0p1pi_truth_MC","GENIE w/RES Truth;Momentum(GeV/c);Counts",200,0.0,2.0);
  THStack *p_pi_res_0p1pi_cuts_stacked = new THStack("p_pi_res_0p1pi_cuts_stacked","RES Single #pi^{#pm} (0p1#pi^{#pm}) Events Pion Momentum Distribution Data to GENIE Comparison;Pion Momentum (GeV/c);Counts");
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //RES single charged pion with a single proton
  //W distributions
  TH1D *W_res_1p1pi_cuts_data = new TH1D("W_res_1p1pi_cuts_data","Data w/RES Cuts;W (GeV/c^{2};Counts",200,1.0,2.5);
  TH1D *W_res_1p1pi_cuts_MC = new TH1D("W_res_1p1pi_cuts_MC","GENIE w/RES Cuts;W (GeV/c^{2};Counts",200,1.0,2.5);
  TH1D *W_res_1p1pi_truth_MC = new TH1D("W_res_1p1pi_truth_MC","GENIE w/RES Truth;W (GeV/c^{2};Counts",200,1.0,2.5);
  THStack *W_res_1p1pi_cuts_stacked = new THStack("W_res_1p1pi_cuts_stacked","RES Single #pi^{#pm} with Single Proton (1p1#pi^{#pm}) Events W Distribution Data to GENIE Comparison;W (GeV/c^{2};Counts");
  //Bjorken-x distributions
  TH1D *xb_res_1p1pi_cuts_data = new TH1D("xb_res_1p1pi_cuts_data","Data w/RES Cuts",200,0,1.5);
  TH1D *xb_res_1p1pi_cuts_MC = new TH1D("xb_res_1p1pi_cuts_MC","GENIE w/RES Cuts",200,0,1.5);
  TH1D *xb_res_1p1pi_truth_MC = new TH1D("xb_res_1p1pi_truth_MC","GENIE w/RES Truth",200,0,1.5);
  THStack *xb_res_1p1pi_cuts_stacked = new THStack("xb_res_1p1pi_cuts_stacked","RES Single #pi^{#pm} with Single Proton (1p1#pi^{#pm}) Events Bjorken-x Distribution Data to GENIE Comparison;x_{B};Counts");
  //Pion momentum distribution
  TH1D *p_pi_res_1p1pi_cuts_data = new TH1D("p_pi_res_1p1pi_cuts_data","Data w/RES Cuts;Momentum (GeV/c);Counts",200,0.0,2.0);
  TH1D *p_pi_res_1p1pi_cuts_MC = new TH1D("p_pi_res_1p1pi_cuts_MC","GENIE w/RES Cuts;Momentum (GeV/c);Counts",200,0.0,2.0);
  TH1D *p_pi_res_1p1pi_truth_MC = new TH1D("p_pi_res_1p1pi_truth_MC","GENIE w/RES Truth;Momentum(GeV/c);Counts",200,0.0,2.0);
  THStack *p_pi_res_1p1pi_cuts_stacked = new THStack("p_pi_res_1p1pi_cuts_stacked","RES Single #pi^{#pm} with Single Proton (1p1#pi^{#pm}) Events Pion Momentum Distribution Data to GENIE Comparison;Pion Momentum (GeV/c);Counts");
  //Proton momentum distribution
  TH1D *p_proton_res_1p1pi_cuts_data = new TH1D("p_proton_res_1p1pi_cuts_data","Data w/RES Cuts;Momentum (GeV/c);Counts",200,0.0,2.0);
  TH1D *p_proton_res_1p1pi_cuts_MC = new TH1D("p_proton_res_1p1pi_cuts_MC","GENIE w/RES Cuts;Momentum (GeV/c);Counts",200,0.0,2.0);
  TH1D *p_proton_res_1p1pi_truth_MC = new TH1D("p_proton_res_1p1pi_truth_MC","GENIE w/RES Truth;Momentum(GeV/c);Counts",200,0.0,2.0);
  THStack *p_proton_res_1p1pi_cuts_stacked = new THStack("p_proton_res_1p1pi_cuts_stacked","RES Single #pi^{#pm} with Single Proton (1p1#pi^{#pm}) Events Proton Momentum Distribution Data to GENIE Comparison;Proton Momentum (GeV/c);Counts");
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //RES pion-neutron

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //RES pion-pion-neutron

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Counters
  int counter_data = 0, counter_MC = 0; //Good electrons
  int counter_qe_data = 0, counter_qe_MC = 0;//QE events
  int counter_mec_data = 0, counter_mec_MC = 0;//MEC events
  int counter_res_data = 0, counter_res_MC = 0, counter_1p2_to_2p1res_MC = 0;//Resonance events
  int counter_dis_data = 0, counter_dis_MC = 0;//DIS events
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Specific topology counters
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Without cuts (raw counts on topology only)
  int counter_0p1pip_raw_data = 0, counter_0p1pip_raw_MC = 0;
  int counter_0p1pim_raw_data = 0, counter_0p1pim_raw_MC = 0;
  int counter_1p0pi_raw_data = 0, counter_1p0pi_raw_MC = 0;
  int counter_1n0pi_raw_data = 0, counter_1n0pi_raw_MC = 0;
  int counter_1n1pip_raw_data = 0, counter_1n1pip_raw_MC = 0;
  int counter_1p1pim_raw_data = 0, counter_1p1pim_raw_MC = 0;
  int counter_2p0pi_raw_data = 0, counter_2p0pi_raw_MC = 0;
  int counter_1p1n_raw_data = 0, counter_1p1n_raw_MC = 0;
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //With cuts
  //Single proton QE (truth is available for processes within GENIE MC)
  int counter_qe_1p0pi_cuts_data = 0, counter_qe_1p0pi_cuts_MC = 0, counter_qe_1p0pi_truth_MC = 0;
  //Single proton MEC
  int counter_mec_1p0pi_cuts_data = 0, counter_mec_1p0pi_cuts_MC = 0, counter_mec_1p0pi_truth_MC = 0;
  //Two proton (proton-proton) MEC
  int counter_mec_2p0pi_cuts_data = 0, counter_mec_2p0pi_cuts_MC = 0, counter_mec_2p0pi_truth_MC = 0;
  //RES single charged pion (0p1p)
  int counter_res_0p1pi_cuts_data = 0, counter_res_0p1pi_cuts_MC = 0, counter_res_0p1pi_truth_MC = 0;
  //RES single charged pion with single proton (1p1p)
  int counter_res_1p1pi_cuts_data = 0, counter_res_1p1pi_cuts_MC = 0, counter_res_1p1pi_truth_MC = 0;

  //Cut strings
  TString qe_1p_cuts  = "protons.size()==1 && y>-0.300 && y<0.300 && x_b>0.8 && x_b<1.2 && W>0.8 && W<1.050 && q2>0. && q2<3.0";
  TString mec_1p_cuts = "protons.size()==1 && y>0.05 && y<0.300 && x_b>0.8 && x_b<1.2 && W>0.8 && W<1.5 && q2>0. && q2<3.0";
  TString mec_2p_cuts = "protons.size()==2 && y>0.05 && y<0.300 && x_b>0.8 && x_b<1.2 && W>0.8 && W<1.5 && q2>0. && q2<3.0";
  TString res_0p1pi_cuts = "protons.size()==0 && ((pips.size()==1 && pims.size()==0) || (pips.size()==0 && pims.size()==1)) && y>0.2 && y<0.6 && x_b>0.4 && x_b<0.8 && W>1.050 && W<2.0 && q2>0.8 && q2<6.0";
  TString res_1p1pi_cuts = "protons.size()==1 && ((pips.size()==1 && pims.size()==0) || (pips.size()==0 && pims.size()==1)) && y>0.2 && y<0.6 && x_b>0.4 && x_b<0.8 && W>1.050 && W<2.0 && q2>0.8 && q2<6.0";

  double p_chi2_fd_cut_lo = -4.;
  double p_chi2_fd_cut_hi = 4.;
  double p_chi2_cd_cut_lo = -4.;
  double p_chi2_cd_cut_hi = 10.;

  double n_chi2_fd_cut_lo = -2.;//These are junk, neutrons dont have chi2 distributions
  double n_chi2_fd_cut_hi = 2.;
  double n_chi2_cd_cut_lo = -3.;
  double n_chi2_cd_cut_hi = 3.;

  double pip_chi2_fd_cut_lo = -2.;
  double pip_chi2_fd_cut_hi = 4.5;
  double pip_chi2_cd_cut_lo = -2.;
  double pip_chi2_cd_cut_hi = 6.;

  double pim_chi2_fd_cut_lo = -2.;
  double pim_chi2_fd_cut_hi = 4.;
  double pim_chi2_cd_cut_lo = -1.;
  double pim_chi2_cd_cut_hi = 5.;

  //std::cout << "I'm above the data" << std::endl;
  //////////////////////////////////////////////////////////////////////////////
  //         DATA
  //////////////////////////////////////////////////////////////////////////////
  while(chain_data.Next())
    {
      beam_energy_data = beamE;
      beam_data.SetE(beam_energy_data);
      beam_data.SetPz(beam_energy_data);
      //std::cout << "I'm in the while loop" << std::endl;
      if( beamE > 1e-4)
	      {
          if( counter_data == 0)
            cout<<"Beam energy set manually: "<<beamE<<endl;
          beam_data.SetE(beamE);
          beam_data.SetPz(beamE);
	      }

      //TLorentzVector el_mc;
      c12_data->mcparts()->setEntry(0);

      //el_mc.SetXYZM(c12_data->mcparts()->getPx(), c12_data->mcparts()->getPy(),c12_data->mcparts()->getPz(),db->GetParticle(11)->Mass());
      
      // get particles by type
      auto electrons=c12_data->getByID(11);
      auto protons=c12_data->getByID(2212);
      auto neutrons=c12_data->getByID(2112);
      auto pips=c12_data->getByID(211);
      auto pims=c12_data->getByID(-211); 
      
      //std::cout << electrons.size() << std::endl;

      if(electrons.size()==1)
	      {	  
          double energy =  (electrons[0]->cal(PCAL)->getEnergy() +  electrons[0]->cal(ECIN)->getEnergy() +  electrons[0]->cal(ECOUT)->getEnergy())/electrons[0]->getP();
          
          bool ecal = ( (electrons[0]->cal(ECIN)->getLv() >= 14 && electrons[0]->cal(ECIN)->getLw() >= 14) && (energy > 0.18 && energy < 0.28) && electrons[0]->getP() > 0.5 && electrons[0]->getP() < 10 );

          SetLorentzVector(el,electrons[0]);

          TLorentzVector q = beam_data - el; //photon  4-vector  
          //Cut out all high radiative correction events
          if (q.E()>omega_threshold_data){continue;}    

          double q2        = -q.M2();
          double x_b       = q2/(2 * mass_p * q.E() ); //Bjorken-x 
          double y         = -q.P() + sqrt( q.E()*q.E() + 2*q.E()*mass_p); //y (Bjorken-y), the scaling variable; from Jourdan, eq. 20
          double vz_e      = electrons[0]->par()->getVz();
          double W         = sqrt(mass_p*mass_p - q2 + 2*q.E()*mass_p); //Invariant mass, assuming electron off of a standing proton
          //std::cout << W << std::endl;
          //std::cout << vz_e << std::endl;
          el_vz_h_data->Fill(vz_e);
          
          //Electron quality cuts
          bool targ_1 = abs( vz_e - (-1.48)) < 1 * 1.25;
          bool targ_2 = abs( vz_e - (-6.3))  < 1 * 1.25;

          //ecal_data->Fill(electrons[0]->getP(),energy/electrons[0]->getP());
          //////////Modified
              //if( !(ecal && (targ_1 || targ_2) ) )//empty rgk run

          //if( !(ecal && vz_e < vz_max && vz_e > vz_min))
          //////////Modified
          //continue;

          //el_vz_h_data->Fill(vz_e);
          ecal_data->Fill(electrons[0]->getP(),energy);

          //remove FT (seems to be still in the MC; FT is not used in RGM)
          //if( electrons[0]->getRegion() == 1000)
            //continue;

          //Electron kinematics monitoring
          //(e,e') cross sections

          theta_data->Fill(el.Theta()*TMath::RadToDeg());
          theta_mom_data->Fill(el.P(),el.Theta()*TMath::RadToDeg());
          q2_data->Fill(q2);
          q2_theta_data->Fill(q2,el.Theta()*TMath::RadToDeg());
          W_data->Fill(W);
          xb_data->Fill(x_b);
          y_data->Fill(y);
          q2_W_data->Fill(W,q2);
          omega_data->Fill(q.E());
          omega_theta_data->Fill(q.E(),el.Theta()*TMath::RadToDeg());
          W_theta_data->Fill(W,el.Theta()*TMath::RadToDeg());
          //Plots for snowmass paper
          q0_vs_q3_data->Fill(q.Pz(),q.E());
          if (protons.size()==1 && pips.size()==0 && pims.size()==0){q0_vs_q3_1pXn0pi_data->Fill(q.Pz(),q.E());}
          if (protons.size()==2 && pips.size()==0 && pims.size()==0){q0_vs_q3_2pXn0pi_data->Fill(q.Pz(),q.E());}

          //Primakoff idea
          if (protons.size()==0 && pips.size()==1 && pims.size()==1)
            {
              SetLorentzVector(pip,pips[0]);
              SetLorentzVector(pim,pims[0]);
              TLorentzVector pions = pip + pim;
              invmass_vs_theta_0pXn1pip1pim_data->Fill(pions.Theta()*TMath::RadToDeg(),pions.M());
            }
          

          if(el.Theta()*TMath::RadToDeg() < 10 && el.Theta()*TMath::RadToDeg() > 5)
            energy_transfer_0_data->Fill(q.E());
          else if(el.Theta()*TMath::RadToDeg() < 15  && el.Theta()*TMath::RadToDeg() >= 10)
            energy_transfer_1_data->Fill(q.E());
          else if(el.Theta()*TMath::RadToDeg() < 25  && el.Theta()*TMath::RadToDeg() >= 15)
            energy_transfer_2_data->Fill(q.E());
          else if(el.Theta()*TMath::RadToDeg() < 35  && el.Theta()*TMath::RadToDeg() >= 25)
            energy_transfer_3_data->Fill(q.E());

          //By sectors
          W_bySector_data->Fill(W);
          q2_bySector_data->Fill(q2);
          //Sector 1
          if (electrons[0]->trk(DC)->getSector() == 1)
            {
              W_S1_data->Fill(W);
              q2_S1_data->Fill(q2);
              W_theta_S1_data->Fill(W,el.Theta()*TMath::RadToDeg());
            }
          //Etc...
          if (electrons[0]->trk(DC)->getSector() == 2)
            {
              W_S2_data->Fill(W);
              q2_S2_data->Fill(q2);
              W_theta_S2_data->Fill(W,el.Theta()*TMath::RadToDeg());
            }
          if (electrons[0]->trk(DC)->getSector() == 3)
            {
              W_S3_data->Fill(W);
              q2_S3_data->Fill(q2);
              W_theta_S3_data->Fill(W,el.Theta()*TMath::RadToDeg());
            }
          if (electrons[0]->trk(DC)->getSector() == 4)
            {
              W_S4_data->Fill(W);
              q2_S4_data->Fill(q2);
              W_theta_S4_data->Fill(W,el.Theta()*TMath::RadToDeg());
            }
          if (electrons[0]->trk(DC)->getSector() == 5)
            {
              W_S5_data->Fill(W);
              q2_S5_data->Fill(q2);
              W_theta_S5_data->Fill(W,el.Theta()*TMath::RadToDeg());
            }
          if (electrons[0]->trk(DC)->getSector() == 6)
            {
              W_S6_data->Fill(W);
              q2_S6_data->Fill(q2);
              W_theta_S6_data->Fill(W,el.Theta()*TMath::RadToDeg());
            }

          //MULTIPLICITIES 
          //Monitor the proton, neutron, and pion multiplicities 
          int num_p = 0;
          int num_n = 0;
          int num_n_nocuts = 0;
          int num_pip = 0;
          int num_pim = 0;
          
          for(int iPr = 0; iPr < protons.size(); iPr++)
            {
              double p_chi2_cd_v = 999999;
              double p_chi2_fd_v = 999999;
              
              if(protons[iPr]->getRegion()==CD)
                {
                  p_chi2_cd_data->Fill(protons[iPr]->par()->getChi2Pid());
                  p_chi2_cd_v = protons[iPr]->par()->getChi2Pid();
                }
            
              if(protons[iPr]->getRegion()==FD)
                {
                  p_chi2_fd_data->Fill(protons[iPr]->par()->getChi2Pid());
                  p_chi2_fd_v = protons[iPr]->par()->getChi2Pid();
                }
              
              if( ( ((p_chi2_fd_v < p_chi2_fd_cut_hi) && (p_chi2_fd_v > p_chi2_fd_cut_lo)) || ((p_chi2_cd_v < p_chi2_cd_cut_hi) && (p_chi2_cd_v > p_chi2_cd_cut_lo)) ) )
                {
                  num_p++;
                }
            }

          for(int iNr = 0; iNr < neutrons.size(); iNr++)
            {
              //num_n_nocuts++;
              n_chi2_nocuts_data->Fill(neutrons[iNr]->par()->getChi2Pid());

              double n_chi2_cd_v = 999999;
              double n_chi2_fd_v = 999999;
              
              if(neutrons[iNr]->getRegion()==CD)
                {
                  n_chi2_cd_data->Fill(neutrons[iNr]->par()->getChi2Pid());
                  n_chi2_cd_v = neutrons[iNr]->par()->getChi2Pid();
                }
            
              if(neutrons[iNr]->getRegion()==FD)
                {
                  n_chi2_fd_data->Fill(neutrons[iNr]->par()->getChi2Pid());
                  n_chi2_fd_v = neutrons[iNr]->par()->getChi2Pid();
                }
              
              if( ( ((n_chi2_fd_v < n_chi2_fd_cut_hi) && (n_chi2_fd_v > n_chi2_fd_cut_lo)) || ((n_chi2_cd_v < n_chi2_cd_cut_hi) && (n_chi2_cd_v > n_chi2_cd_cut_lo)) ) )
                {
                  num_n++;
                }
            }
	  
          for(int iPim = 0; iPim < pims.size(); iPim++)
            {
              double pim_chi2_cd_v = 999999;
              double pim_chi2_fd_v = 999999;
              
              if(pims[iPim]->getRegion()==CD)
                {
                  pim_chi2_cd_data->Fill(pims[iPim]->par()->getChi2Pid());
                  pim_chi2_cd_v = pims[iPim]->par()->getChi2Pid();
                }
            
              if(pims[iPim]->getRegion()==FD)
                {
                  pim_chi2_fd_data->Fill(pims[iPim]->par()->getChi2Pid());
                  pim_chi2_fd_v = pims[iPim]->par()->getChi2Pid();
                }
              
              if( ( ((pim_chi2_fd_v < pim_chi2_fd_cut_hi) && (pim_chi2_fd_v > pim_chi2_fd_cut_lo)) || ((pim_chi2_cd_v < pim_chi2_cd_cut_hi) && (pim_chi2_cd_v > pim_chi2_cd_cut_lo)) ) )
                {
                  num_pim++;
                }
            }

          for(int iPip = 0; iPip < pips.size(); iPip++)
            {
              double pip_chi2_cd_v = 999999;
              double pip_chi2_fd_v = 999999;
              
              if(pips[iPip]->getRegion()==CD)
                {
                  pip_chi2_cd_data->Fill(pips[iPip]->par()->getChi2Pid());
                  pip_chi2_cd_v = pips[iPip]->par()->getChi2Pid();
                }           
            
              if(pips[iPip]->getRegion()==FD)
                {
                  pip_chi2_fd_data->Fill(pips[iPip]->par()->getChi2Pid());
                  pip_chi2_fd_v = pips[iPip]->par()->getChi2Pid();
                }
            
              if( ( ((pip_chi2_fd_v < pip_chi2_fd_cut_hi) && (pip_chi2_fd_v > pip_chi2_fd_cut_lo)) || ((pip_chi2_cd_v < pip_chi2_cd_cut_hi) && (pip_chi2_cd_v > pip_chi2_cd_cut_lo)) ) )
                {
                  num_pip++;
                }
            }
        
          mult_p_data->Fill(num_p);
          mult_n_data->Fill(num_n);
          mult_n_nocuts_data->Fill(neutrons.size());
          mult_pim_data->Fill(num_pim);
          mult_pip_data->Fill(num_pip);
          mult_p_pis_data->Fill(num_p,num_pim+num_pip);
          mult_p_n_data->Fill(num_p,neutrons.size());
        
          //INVARIANT MASS
          //Calculating the Invariant mass of the p + pi- system 
          for(int iPr = 0; iPr < protons.size(); iPr++)
            {
              for(int iPim = 0; iPim < pims.size(); iPim++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  SetLorentzVector(pim,pims[iPim]);
                  
                  double vz_e_p = electrons[0]->par()->getVz() - protons[iPr]->par()->getVz();
                  TLorentzVector inv_mass = pr + pim; //missing 4-vector
                  
                  invariant_mass_data->Fill(inv_mass.M());
                }
            } 

          //Proton variables
          for(int iPr = 0; iPr < protons.size(); iPr++)
            {
              SetLorentzVector(pr,protons[iPr]);
              double p_proton = pr.P();
              //Get p_miss vector stuff
              TLorentzVector p_miss = pr - q;
              double p_miss_event = p_miss.P();
              double p_miss_theta = p_miss.Theta()*TMath::RadToDeg();
              //Fill histograms
              p_proton_data->Fill(p_proton);
              p_miss_data->Fill(p_miss_event);
              p_miss_theta_data->Fill(p_miss_theta);
              p_miss_theta_vs_omega_data->Fill(q.E(),p_miss_theta);
            }


          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          //1D Data histograms of variables as a function of reconstructed channel
          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          q2_byChannel_data->Fill(q2);
          W_byChannel_data->Fill(W);
          //0p1n0pi
          if ( 
               protons.size()==0
            && neutrons.size()==1
            /*&& (   //Must check that each neutron has a good chi2 value range in either the forward or central detectors
                   ((neutrons[0]->getRegion()==FD) && (neutrons[0]->par()->getChi2Pid() < n_chi2_fd_cut_hi) && (neutrons[0]->par()->getChi2Pid() > n_chi2_fd_cut_lo))
                || ((neutrons[0]->getRegion()==CD) && (neutrons[0]->par()->getChi2Pid() < n_chi2_cd_cut_hi) && (neutrons[0]->par()->getChi2Pid() > n_chi2_cd_cut_lo))
               )*/
            && pips.size()==0
            && pims.size()==0 
             )
            {
              double proton_totrecoe = 0., neutron_totrecoe = 0., pim_totrecoe = 0., pip_totrecoe=0.;
              for(int iNr = 0; iNr < neutrons.size(); iNr++)
                {
                   SetLorentzVector(nr,neutrons[iNr]);
                  neutron_totrecoe = neutron_totrecoe + nr.E() - mass_n;
                }
              double totrecoe = proton_totrecoe + neutron_totrecoe + pim_totrecoe + pip_totrecoe + el.E();
              q2_0p1n0pi_data->Fill(q2);
              totrecoe_0p1n0pi_data->Fill(totrecoe);
              W_0p1n0pi_data->Fill(W);
            }
          //1pXn0pi--no neutron/neutral particle cuts
          if ( 
               protons.size()==1
            && (   //Must check that each proton has a good chi2 value range in either the forward or central detectors
                   ((protons[0]->getRegion()==FD) && (protons[0]->par()->getChi2Pid() < p_chi2_fd_cut_hi) && (protons[0]->par()->getChi2Pid() > p_chi2_fd_cut_lo))
                || ((protons[0]->getRegion()==CD) && (protons[0]->par()->getChi2Pid() < p_chi2_cd_cut_hi) && (protons[0]->par()->getChi2Pid() > p_chi2_cd_cut_lo))
               )
            && pips.size()==0
            && pims.size()==0 
             )
            {
              double proton_totrecoe = 0., neutron_totrecoe = 0., pim_totrecoe = 0., pip_totrecoe=0.;
              for(int iPr = 0; iPr < protons.size(); iPr++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  proton_totrecoe = proton_totrecoe + pr.E() - mass_p;
                }
              double totrecoe = proton_totrecoe + neutron_totrecoe + pim_totrecoe + pip_totrecoe + el.E();
              q2_1pXn0pi_data->Fill(q2);
              totrecoe_1pXn0pi_data->Fill(totrecoe);
              W_1pXn0pi_data->Fill(W);
              //Find invariant mass of the system (help to see if it comes from a true event?)
              TLorentzVector inv_mass_1pXn0pi = pr + pip + pim - q;
              nuc_inv_mass_W_1pXn0pi_data->Fill(W,inv_mass_1pXn0pi.M());
              nuc_inv_mass_omega_1pXn0pi_data->Fill(q.E(),inv_mass_1pXn0pi.M());
              nuc_inv_mass_q2_1pXn0pi_data->Fill(q2,inv_mass_1pXn0pi.M());
              nuc_inv_mass_xb_1pXn0pi_data->Fill(x_b,inv_mass_1pXn0pi.M());
              nuc_inv_mass_y_1pXn0pi_data->Fill(y,inv_mass_1pXn0pi.M());
              //Fill 2D histograms
              totrecoe_W_1pXn0pi_data->Fill(W,totrecoe);
              totrecoe_omega_1pXn0pi_data->Fill(q.E(),totrecoe);
              totrecoe_q2_1pXn0pi_data->Fill(q2,totrecoe);
              totrecoe_xb_1pXn0pi_data->Fill(x_b,totrecoe);
              totrecoe_y_1pXn0pi_data->Fill(y,totrecoe);
              W_q2_1pXn0pi_data->Fill(q2,W);
              W_xb_1pXn0pi_data->Fill(x_b,W);
              W_y_1pXn0pi_data->Fill(y,W);
              q2_xb_1pXn0pi_data->Fill(x_b,q2);
              q2_y_1pXn0pi_data->Fill(y,q2);
              y_xb_1pXn0pi_data->Fill(x_b,y);
            }
          //1pXn1pim--no neutron/neutral particle cuts
          if ( 
               protons.size()==1
            && (   //Must check that each proton has a good chi2 value range in either the forward or central detectors
                   ((protons[0]->getRegion()==FD) && (protons[0]->par()->getChi2Pid() < p_chi2_fd_cut_hi) && (protons[0]->par()->getChi2Pid() > p_chi2_fd_cut_lo))
                || ((protons[0]->getRegion()==CD) && (protons[0]->par()->getChi2Pid() < p_chi2_cd_cut_hi) && (protons[0]->par()->getChi2Pid() > p_chi2_cd_cut_lo))
               )
            && pips.size()==0
            && pims.size()==1 
            && (   //Must check that each pim has a good chi2 value range in either the forward or central detectors
                   ((pims[0]->getRegion()==FD) && (pims[0]->par()->getChi2Pid() < pim_chi2_fd_cut_hi) && (pims[0]->par()->getChi2Pid() > pim_chi2_fd_cut_lo))
                || ((pims[0]->getRegion()==CD) && (pims[0]->par()->getChi2Pid() < pim_chi2_cd_cut_hi) && (pims[0]->par()->getChi2Pid() > pim_chi2_cd_cut_lo))
               )
             )
            {              
              double proton_totrecoe = 0., neutron_totrecoe = 0., pim_totrecoe = 0., pip_totrecoe=0.;
              for(int iPr = 0; iPr < protons.size(); iPr++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  proton_totrecoe = proton_totrecoe + pr.E() - mass_p;
                }
              for(int iPim = 0; iPim < pims.size(); iPim++)
                {
                  SetLorentzVector(pim,pims[iPim]);
                  pim_totrecoe = pim_totrecoe + pim.E();
                }
              double totrecoe = proton_totrecoe + neutron_totrecoe + pim_totrecoe + pip_totrecoe + el.E();
              q2_1pXn1pim_data->Fill(q2);
              totrecoe_1pXn1pim_data->Fill(totrecoe);
              W_1pXn1pim_data->Fill(W);
              //Find invariant mass of the system (help to see if it comes from a true event?)
              TLorentzVector inv_mass_1pXn1pim = pr + pip + pim - q;
              nuc_inv_mass_W_1pXn1pim_data->Fill(W,inv_mass_1pXn1pim.M());
              nuc_inv_mass_omega_1pXn1pim_data->Fill(q.E(),inv_mass_1pXn1pim.M());
              nuc_inv_mass_q2_1pXn1pim_data->Fill(q2,inv_mass_1pXn1pim.M());
              nuc_inv_mass_xb_1pXn1pim_data->Fill(x_b,inv_mass_1pXn1pim.M());
              nuc_inv_mass_y_1pXn1pim_data->Fill(y,inv_mass_1pXn1pim.M());
              //Fill 2D histograms
              totrecoe_W_1pXn1pim_data->Fill(W,totrecoe);
              totrecoe_omega_1pXn1pim_data->Fill(q.E(),totrecoe);
              totrecoe_q2_1pXn1pim_data->Fill(q2,totrecoe);
              totrecoe_xb_1pXn1pim_data->Fill(x_b,totrecoe);
              totrecoe_y_1pXn1pim_data->Fill(y,totrecoe);
              W_q2_1pXn1pim_data->Fill(q2,W);
              W_xb_1pXn1pim_data->Fill(x_b,W);
              W_y_1pXn1pim_data->Fill(y,W);
              q2_xb_1pXn1pim_data->Fill(x_b,q2);
              q2_y_1pXn1pim_data->Fill(y,q2);
              y_xb_1pXn1pim_data->Fill(x_b,y);
            }
          //0p1n1pip
          if ( 
               protons.size()==0
            && neutrons.size()==1
            /*&& (   //Must check that each neutron has a good chi2 value range in either the forward or central detectors
                   ((neutrons[0]->getRegion()==FD) && (neutrons[0]->par()->getChi2Pid() < n_chi2_fd_cut_hi) && (neutrons[0]->par()->getChi2Pid() > n_chi2_fd_cut_lo))
                || ((neutrons[0]->getRegion()==CD) && (neutrons[0]->par()->getChi2Pid() < n_chi2_cd_cut_hi) && (neutrons[0]->par()->getChi2Pid() > n_chi2_cd_cut_lo))
               )*/
            && pips.size()==1
            && (   //Must check that each pim has a good chi2 value range in either the forward or central detectors
                   ((pips[0]->getRegion()==FD) && (pips[0]->par()->getChi2Pid() < pip_chi2_fd_cut_hi) && (pips[0]->par()->getChi2Pid() > pip_chi2_fd_cut_lo))
                || ((pips[0]->getRegion()==CD) && (pips[0]->par()->getChi2Pid() < pip_chi2_cd_cut_hi) && (pips[0]->par()->getChi2Pid() > pip_chi2_cd_cut_lo))
               )
            && pims.size()==0 
             )
            {
              double proton_totrecoe = 0., neutron_totrecoe = 0., pim_totrecoe = 0., pip_totrecoe=0.;
              for(int iNr = 0; iNr < neutrons.size(); iNr++)
                {
                  SetLorentzVector(nr,neutrons[iNr]);
                  neutron_totrecoe = neutron_totrecoe + nr.E() - mass_n;
                }
              for(int iPip = 0; iPip < pips.size(); iPip++)
                {
                  SetLorentzVector(pip,pips[iPip]);
                  pip_totrecoe = pip_totrecoe + pip.E();
                }
              double totrecoe = proton_totrecoe + neutron_totrecoe + pim_totrecoe + pip_totrecoe + el.E();
              q2_0p1n1pip_data->Fill(q2);
              totrecoe_0p1n1pip_data->Fill(totrecoe);
              W_0p1n1pip_data->Fill(W);
            }
          //2pXn0pi--no neutron/neutral particle cuts
          if ( 
               protons.size()==2
            && (   //Must check that the first proton has a good chi2 value range in either the forward or central detectors
                   ((protons[0]->getRegion()==FD) && (protons[0]->par()->getChi2Pid() < p_chi2_fd_cut_hi) && (protons[0]->par()->getChi2Pid() > p_chi2_fd_cut_lo))
                || ((protons[0]->getRegion()==CD) && (protons[0]->par()->getChi2Pid() < p_chi2_cd_cut_hi) && (protons[0]->par()->getChi2Pid() > p_chi2_cd_cut_lo))
               )
            && (   //Must also check that the second proton has a good chi2 value range in either the forward or central detectors
                   ((protons[1]->getRegion()==FD) && (protons[1]->par()->getChi2Pid() < p_chi2_fd_cut_hi) && (protons[1]->par()->getChi2Pid() > p_chi2_fd_cut_lo))
                || ((protons[1]->getRegion()==CD) && (protons[1]->par()->getChi2Pid() < p_chi2_cd_cut_hi) && (protons[1]->par()->getChi2Pid() > p_chi2_cd_cut_lo))
               )
            && pips.size()==0
            && pims.size()==0 
             )
            {
              double proton_totrecoe = 0., neutron_totrecoe = 0., pim_totrecoe = 0., pip_totrecoe=0.;
              for(int iPr = 0; iPr < protons.size(); iPr++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  proton_totrecoe = proton_totrecoe + pr.E() - mass_p;
                }
              double totrecoe = proton_totrecoe + neutron_totrecoe + pim_totrecoe + pip_totrecoe + el.E();
              q2_2pXn0pi_data->Fill(q2);
              totrecoe_2pXn0pi_data->Fill(totrecoe);
              W_2pXn0pi_data->Fill(W);
            }
          //1p1n0pi
          if ( 
               protons.size()==1
            && (   //Must check that each proton has a good chi2 value range in either the forward or central detectors
                   ((protons[0]->getRegion()==FD) && (protons[0]->par()->getChi2Pid() < p_chi2_fd_cut_hi) && (protons[0]->par()->getChi2Pid() > p_chi2_fd_cut_lo))
                || ((protons[0]->getRegion()==CD) && (protons[0]->par()->getChi2Pid() < p_chi2_cd_cut_hi) && (protons[0]->par()->getChi2Pid() > p_chi2_cd_cut_lo))
               )
            && neutrons.size()==1
            /*&& (   //Must check that each neutron has a good chi2 value range in either the forward or central detectors
                   ((neutrons[0]->getRegion()==FD) && (neutrons[0]->par()->getChi2Pid() < n_chi2_fd_cut_hi) && (neutrons[0]->par()->getChi2Pid() > n_chi2_fd_cut_lo))
                || ((neutrons[0]->getRegion()==CD) && (neutrons[0]->par()->getChi2Pid() < n_chi2_cd_cut_hi) && (neutrons[0]->par()->getChi2Pid() > n_chi2_cd_cut_lo))
               )*/
            && pips.size()==0
            && pims.size()==0 
             )
            {
              double proton_totrecoe = 0., neutron_totrecoe = 0., pim_totrecoe = 0., pip_totrecoe=0.;
              for(int iPr = 0; iPr < protons.size(); iPr++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  proton_totrecoe = proton_totrecoe + pr.E() - mass_p;
                }
              for(int iNr = 0; iNr < neutrons.size(); iNr++)
                {
                  SetLorentzVector(nr,neutrons[iNr]);
                  neutron_totrecoe = neutron_totrecoe + nr.E() - mass_n;
                }
              double totrecoe = proton_totrecoe + neutron_totrecoe + pim_totrecoe + pip_totrecoe + el.E();
              q2_1p1n0pi_data->Fill(q2);
              totrecoe_1p1n0pi_data->Fill(totrecoe);
              W_1p1n0pi_data->Fill(W);
            }
          //1p1n1pi
          //if ( protons.size()==1 && neutrons.size()==1 && ((pips.size()==1 && pims.size()==0) || (pips.size()==0 && pims.size()==1)) )
          if ( 
               protons.size()==1
            && (   //Must check that each proton has a good chi2 value range in either the forward or central detectors
                   ((protons[0]->getRegion()==FD) && (protons[0]->par()->getChi2Pid() < p_chi2_fd_cut_hi) && (protons[0]->par()->getChi2Pid() > p_chi2_fd_cut_lo))
                || ((protons[0]->getRegion()==CD) && (protons[0]->par()->getChi2Pid() < p_chi2_cd_cut_hi) && (protons[0]->par()->getChi2Pid() > p_chi2_cd_cut_lo))
               )
            && neutrons.size()==1
            /*&& (   //Must check that each neutron has a good chi2 value range in either the forward or central detectors
                   ((neutrons[0]->getRegion()==FD) && (neutrons[0]->par()->getChi2Pid() < n_chi2_fd_cut_hi) && (neutrons[0]->par()->getChi2Pid() > n_chi2_fd_cut_lo))
                || ((neutrons[0]->getRegion()==CD) && (neutrons[0]->par()->getChi2Pid() < n_chi2_cd_cut_hi) && (neutrons[0]->par()->getChi2Pid() > n_chi2_cd_cut_lo))
               )*/
            && (
                   ( pips.size()==1 && pims.size()==0
                   && (   //Must check that each pim has a good chi2 value range in either the forward or central detectors
                         ((pips[0]->getRegion()==FD) && (pips[0]->par()->getChi2Pid() < pip_chi2_fd_cut_hi) && (pips[0]->par()->getChi2Pid() > pip_chi2_fd_cut_lo))
                      || ((pips[0]->getRegion()==CD) && (pips[0]->par()->getChi2Pid() < pip_chi2_cd_cut_hi) && (pips[0]->par()->getChi2Pid() > pip_chi2_cd_cut_lo))
                      )
                   )
                || ( pips.size()==0 && pims.size()==1
                   && (   //Must check that each pim has a good chi2 value range in either the forward or central detectors
                         ((pims[0]->getRegion()==FD) && (pims[0]->par()->getChi2Pid() < pip_chi2_fd_cut_hi) && (pims[0]->par()->getChi2Pid() > pip_chi2_fd_cut_lo))
                      || ((pims[0]->getRegion()==CD) && (pims[0]->par()->getChi2Pid() < pip_chi2_cd_cut_hi) && (pims[0]->par()->getChi2Pid() > pip_chi2_cd_cut_lo))
                      )
                   )
               )
             )
            {
              double proton_totrecoe = 0., neutron_totrecoe = 0., pim_totrecoe = 0., pip_totrecoe=0.;
              for(int iPr = 0; iPr < protons.size(); iPr++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  proton_totrecoe = proton_totrecoe + pr.E() - mass_p;
                }
              for(int iNr = 0; iNr < neutrons.size(); iNr++)
                {
                  SetLorentzVector(nr,neutrons[iNr]);
                  neutron_totrecoe = neutron_totrecoe + nr.E() - mass_n;
                }
              for(int iPip = 0; iPip < pips.size(); iPip++)
                {
                  SetLorentzVector(pip,pips[iPip]);
                  pip_totrecoe = pip_totrecoe + pip.E();
                }
              for(int iPim = 0; iPim < pims.size(); iPim++)
                {
                  SetLorentzVector(pim,pims[iPim]);
                  pim_totrecoe = pim_totrecoe + pim.E();
                }
              double totrecoe = proton_totrecoe + neutron_totrecoe + pim_totrecoe + pip_totrecoe + el.E();
              q2_1p1n1pi_data->Fill(q2);
              totrecoe_1p1n1pi_data->Fill(totrecoe);
              W_1p1n1pi_data->Fill(W);
            }
          //1pXn2pi--no neutron/neutral particle cuts
          //if ( protons.size()==1 && ((pips.size()==2 && pims.size()==0) || (pips.size()==1 && pims.size()==1) || (pips.size()==0 && pims.size()==2)) )
          if (
               protons.size()==1
            && (   //Must check that each proton has a good chi2 value range in either the forward or central detectors
                   ((protons[0]->getRegion()==FD) && (protons[0]->par()->getChi2Pid() < p_chi2_fd_cut_hi) && (protons[0]->par()->getChi2Pid() > p_chi2_fd_cut_lo))
                || ((protons[0]->getRegion()==CD) && (protons[0]->par()->getChi2Pid() < p_chi2_cd_cut_hi) && (protons[0]->par()->getChi2Pid() > p_chi2_cd_cut_lo))
               )
            && (
                  ( pips.size()==2 && pims.size()==0
                  && (   //Must check that each pim has a good chi2 value range in either the forward or central detectors
                       (
                           ((pips[0]->getRegion()==FD) && (pips[0]->par()->getChi2Pid() < pip_chi2_fd_cut_hi) && (pips[0]->par()->getChi2Pid() > pip_chi2_fd_cut_lo))
                        || ((pips[0]->getRegion()==CD) && (pips[0]->par()->getChi2Pid() < pip_chi2_cd_cut_hi) && (pips[0]->par()->getChi2Pid() > pip_chi2_cd_cut_lo))
                       )
                       &&
                       (
                           ((pips[1]->getRegion()==FD) && (pips[1]->par()->getChi2Pid() < pip_chi2_fd_cut_hi) && (pips[1]->par()->getChi2Pid() > pip_chi2_fd_cut_lo))
                        || ((pips[1]->getRegion()==CD) && (pips[1]->par()->getChi2Pid() < pip_chi2_cd_cut_hi) && (pips[1]->par()->getChi2Pid() > pip_chi2_cd_cut_lo))
                       )
                     )
                  )
               || ( pips.size()==1 && pims.size()==1
                  && (   //Must check that each pim has a good chi2 value range in either the forward or central detectors
                       (
                           ((pips[0]->getRegion()==FD) && (pips[0]->par()->getChi2Pid() < pip_chi2_fd_cut_hi) && (pips[0]->par()->getChi2Pid() > pip_chi2_fd_cut_lo))
                        || ((pips[0]->getRegion()==CD) && (pips[0]->par()->getChi2Pid() < pip_chi2_cd_cut_hi) && (pips[0]->par()->getChi2Pid() > pip_chi2_cd_cut_lo))
                       )
                       &&
                       (
                           ((pims[0]->getRegion()==FD) && (pims[0]->par()->getChi2Pid() < pip_chi2_fd_cut_hi) && (pims[0]->par()->getChi2Pid() > pip_chi2_fd_cut_lo))
                        || ((pims[0]->getRegion()==CD) && (pims[0]->par()->getChi2Pid() < pip_chi2_cd_cut_hi) && (pims[0]->par()->getChi2Pid() > pip_chi2_cd_cut_lo))
                       )
                     )
                  )
               || ( pips.size()==0 && pims.size()==2
                  && (   //Must check that each pim has a good chi2 value range in either the forward or central detectors
                       (
                           ((pims[0]->getRegion()==FD) && (pims[0]->par()->getChi2Pid() < pip_chi2_fd_cut_hi) && (pims[0]->par()->getChi2Pid() > pip_chi2_fd_cut_lo))
                        || ((pims[0]->getRegion()==CD) && (pims[0]->par()->getChi2Pid() < pip_chi2_cd_cut_hi) && (pims[0]->par()->getChi2Pid() > pip_chi2_cd_cut_lo))
                       )
                       &&
                       (
                           ((pims[1]->getRegion()==FD) && (pims[1]->par()->getChi2Pid() < pip_chi2_fd_cut_hi) && (pims[1]->par()->getChi2Pid() > pip_chi2_fd_cut_lo))
                        || ((pims[1]->getRegion()==CD) && (pims[1]->par()->getChi2Pid() < pip_chi2_cd_cut_hi) && (pims[1]->par()->getChi2Pid() > pip_chi2_cd_cut_lo))
                       )
                     )
                  )
               )
             )
            {
              double proton_totrecoe = 0., neutron_totrecoe = 0., pim_totrecoe = 0., pip_totrecoe=0.;
              for(int iPr = 0; iPr < protons.size(); iPr++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  proton_totrecoe = proton_totrecoe + pr.E() - mass_p;
                }
              for(int iNr = 0; iNr < neutrons.size(); iNr++)
                {
                  SetLorentzVector(nr,neutrons[iNr]);
                  neutron_totrecoe = neutron_totrecoe + nr.E() - mass_n;
                }
              for(int iPip = 0; iPip < pips.size(); iPip++)
                {
                  SetLorentzVector(pip,pips[iPip]);
                  pip_totrecoe = pip_totrecoe + pip.E();
                }
              for(int iPim = 0; iPim < pims.size(); iPim++)
                {
                  SetLorentzVector(pim,pims[iPim]);
                  pim_totrecoe = pim_totrecoe + pim.E();
                }
              double totrecoe = proton_totrecoe + neutron_totrecoe + pim_totrecoe + pip_totrecoe + el.E();
              q2_1pXn2pi_data->Fill(q2);
              totrecoe_1pXn2pi_data->Fill(totrecoe);
              W_1pXn2pi_data->Fill(W);
            }
          //0p1n2pi
          if (
               protons.size()==0
            && neutrons.size()==1
            /*&& (   //Must check that each neutron has a good chi2 value range in either the forward or central detectors
                   ((neutrons[0]->getRegion()==FD) && (neutrons[0]->par()->getChi2Pid() < n_chi2_fd_cut_hi) && (neutrons[0]->par()->getChi2Pid() > n_chi2_fd_cut_lo))
                || ((neutrons[0]->getRegion()==CD) && (neutrons[0]->par()->getChi2Pid() < n_chi2_cd_cut_hi) && (neutrons[0]->par()->getChi2Pid() > n_chi2_cd_cut_lo))
               )*/
            && (
                  ( pips.size()==2 && pims.size()==0
                  && (   //Must check that each pim has a good chi2 value range in either the forward or central detectors
                       (
                           ((pips[0]->getRegion()==FD) && (pips[0]->par()->getChi2Pid() < pip_chi2_fd_cut_hi) && (pips[0]->par()->getChi2Pid() > pip_chi2_fd_cut_lo))
                        || ((pips[0]->getRegion()==CD) && (pips[0]->par()->getChi2Pid() < pip_chi2_cd_cut_hi) && (pips[0]->par()->getChi2Pid() > pip_chi2_cd_cut_lo))
                       )
                       &&
                       (
                           ((pips[1]->getRegion()==FD) && (pips[1]->par()->getChi2Pid() < pip_chi2_fd_cut_hi) && (pips[1]->par()->getChi2Pid() > pip_chi2_fd_cut_lo))
                        || ((pips[1]->getRegion()==CD) && (pips[1]->par()->getChi2Pid() < pip_chi2_cd_cut_hi) && (pips[1]->par()->getChi2Pid() > pip_chi2_cd_cut_lo))
                       )
                     )
                  )
               || ( pips.size()==1 && pims.size()==1
                  && (   //Must check that each pim has a good chi2 value range in either the forward or central detectors
                       (
                           ((pips[0]->getRegion()==FD) && (pips[0]->par()->getChi2Pid() < pip_chi2_fd_cut_hi) && (pips[0]->par()->getChi2Pid() > pip_chi2_fd_cut_lo))
                        || ((pips[0]->getRegion()==CD) && (pips[0]->par()->getChi2Pid() < pip_chi2_cd_cut_hi) && (pips[0]->par()->getChi2Pid() > pip_chi2_cd_cut_lo))
                       )
                       &&
                       (
                           ((pims[0]->getRegion()==FD) && (pims[0]->par()->getChi2Pid() < pip_chi2_fd_cut_hi) && (pims[0]->par()->getChi2Pid() > pip_chi2_fd_cut_lo))
                        || ((pims[0]->getRegion()==CD) && (pims[0]->par()->getChi2Pid() < pip_chi2_cd_cut_hi) && (pims[0]->par()->getChi2Pid() > pip_chi2_cd_cut_lo))
                       )
                     )
                  )
               || ( pips.size()==0 && pims.size()==2
                  && (   //Must check that each pim has a good chi2 value range in either the forward or central detectors
                       (
                           ((pims[0]->getRegion()==FD) && (pims[0]->par()->getChi2Pid() < pip_chi2_fd_cut_hi) && (pims[0]->par()->getChi2Pid() > pip_chi2_fd_cut_lo))
                        || ((pims[0]->getRegion()==CD) && (pims[0]->par()->getChi2Pid() < pip_chi2_cd_cut_hi) && (pims[0]->par()->getChi2Pid() > pip_chi2_cd_cut_lo))
                       )
                       &&
                       (
                           ((pims[1]->getRegion()==FD) && (pims[1]->par()->getChi2Pid() < pip_chi2_fd_cut_hi) && (pims[1]->par()->getChi2Pid() > pip_chi2_fd_cut_lo))
                        || ((pims[1]->getRegion()==CD) && (pims[1]->par()->getChi2Pid() < pip_chi2_cd_cut_hi) && (pims[1]->par()->getChi2Pid() > pip_chi2_cd_cut_lo))
                       )
                     )
                  )
               )
             )
            {
              double proton_totrecoe = 0., neutron_totrecoe = 0., pim_totrecoe = 0., pip_totrecoe=0.;
              for(int iNr = 0; iNr < neutrons.size(); iNr++)
                {
                  SetLorentzVector(nr,neutrons[iNr]);
                  neutron_totrecoe = neutron_totrecoe + nr.E() - mass_n;
                }
              for(int iPip = 0; iPip < pips.size(); iPip++)
                {
                  SetLorentzVector(pip,pips[iPip]);
                  pip_totrecoe = pip_totrecoe + pip.E();
                }
              for(int iPim = 0; iPim < pims.size(); iPim++)
                {
                  SetLorentzVector(pim,pims[iPim]);
                  pim_totrecoe = pim_totrecoe + pim.E();
                }
              double totrecoe = proton_totrecoe + neutron_totrecoe + pim_totrecoe + pip_totrecoe + el.E();
              q2_0p1n2pi_data->Fill(q2);
              totrecoe_0p1n2pi_data->Fill(totrecoe);
              W_0p1n2pi_data->Fill(W);
            }

          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          //2D Histograms of particular topologies without cuts (what should the cuts be)
          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          //0p1pip
          if (protons.size()==0 && pips.size()==1 && pims.size()==0)
            {
              //Find invariant mass of the system (help to see if it comes from a true event?)
              TLorentzVector inv_mass_0p_1pip = pr + pip + pim - q;
              nuc_inv_mass_W_0p_1pip_data->Fill(W,inv_mass_0p_1pip.M());
              nuc_inv_mass_omega_0p_1pip_data->Fill(q.E(),inv_mass_0p_1pip.M());
              nuc_inv_mass_q2_0p_1pip_data->Fill(q2,inv_mass_0p_1pip.M());
              nuc_inv_mass_xb_0p_1pip_data->Fill(x_b,inv_mass_0p_1pip.M());
              nuc_inv_mass_y_0p_1pip_data->Fill(y,inv_mass_0p_1pip.M());
              //Run over particle energies
              double proton_totrecoe = 0., pim_totrecoe = 0., pip_totrecoe=0.;
              for(int iPr = 0; iPr < protons.size(); iPr++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  proton_totrecoe = proton_totrecoe + pr.E() - mass_p;
                }
              for(int iPim = 0; iPim < pims.size(); iPim++)
                {
                  SetLorentzVector(pim,pims[iPim]);
                  pim_totrecoe = pim_totrecoe + pim.E();
                }
              for(int iPip = 0; iPip < pips.size(); iPip++)
                {
                  SetLorentzVector(pip,pips[iPip]);
                  pip_totrecoe = pip_totrecoe + pip.E() - mass_pip;
                }
              //Find the Total Kinetic Energy of the particles in the final state
              double totrecoe = proton_totrecoe + pim_totrecoe + pip_totrecoe + el.E();
              //Fill 2D histograms
              totrecoe_W_0p_1pip_data->Fill(W,totrecoe);
              totrecoe_omega_0p_1pip_data->Fill(q.E(),totrecoe);
              totrecoe_q2_0p_1pip_data->Fill(q2,totrecoe);
              totrecoe_xb_0p_1pip_data->Fill(x_b,totrecoe);
              totrecoe_y_0p_1pip_data->Fill(y,totrecoe);
              W_q2_0p_1pip_data->Fill(q2,W);
              W_xb_0p_1pip_data->Fill(x_b,W);
              W_y_0p_1pip_data->Fill(y,W);
              q2_xb_0p_1pip_data->Fill(x_b,q2);
              q2_y_0p_1pip_data->Fill(y,q2);
              y_xb_0p_1pip_data->Fill(x_b,y);
            }
          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          //0p1pim
          if (protons.size()==0 && pips.size()==0 && pims.size()==1)
            {
              //Find invariant mass of the system (help to see if it comes from a true event?)
              TLorentzVector inv_mass_0p_1pim = pr + pip + pim - q;
              nuc_inv_mass_W_0p_1pim_data->Fill(W,inv_mass_0p_1pim.M());
              nuc_inv_mass_omega_0p_1pim_data->Fill(q.E(),inv_mass_0p_1pim.M());
              nuc_inv_mass_q2_0p_1pim_data->Fill(q2,inv_mass_0p_1pim.M());
              nuc_inv_mass_xb_0p_1pim_data->Fill(x_b,inv_mass_0p_1pim.M());
              nuc_inv_mass_y_0p_1pim_data->Fill(y,inv_mass_0p_1pim.M());
              //Run over particle energies
              double proton_totrecoe = 0., pim_totrecoe = 0., pip_totrecoe=0.;
              for(int iPr = 0; iPr < protons.size(); iPr++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  proton_totrecoe = proton_totrecoe + pr.E() - mass_p;
                }
              for(int iPim = 0; iPim < pims.size(); iPim++)
                {
                  SetLorentzVector(pim,pims[iPim]);
                  pim_totrecoe = pim_totrecoe + pim.E();
                }
              for(int iPip = 0; iPip < pips.size(); iPip++)
                {
                  SetLorentzVector(pip,pips[iPip]);
                  pip_totrecoe = pip_totrecoe + pip.E() - mass_pip;
                }
              double totrecoe = proton_totrecoe + pim_totrecoe + pip_totrecoe + el.E();
              //Fill 2D histograms
              totrecoe_W_0p_1pim_data->Fill(W,totrecoe);
              totrecoe_omega_0p_1pim_data->Fill(q.E(),totrecoe);
              totrecoe_q2_0p_1pim_data->Fill(q2,totrecoe);
              totrecoe_xb_0p_1pim_data->Fill(x_b,totrecoe);
              totrecoe_y_0p_1pim_data->Fill(y,totrecoe);
              W_q2_0p_1pim_data->Fill(q2,W);
              W_xb_0p_1pim_data->Fill(x_b,W);
              W_y_0p_1pim_data->Fill(y,W);
              q2_xb_0p_1pim_data->Fill(x_b,q2);
              q2_y_0p_1pim_data->Fill(y,q2);
              y_xb_0p_1pim_data->Fill(x_b,y);
            }
          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          //1p0pi
          /*if (protons.size()==1 && pips.size()==0 && pims.size()==0)
            {
              //Find invariant mass of the system (help to see if it comes from a true event?)
              TLorentzVector inv_mass_1pXn0pi = pr + pip + pim - q;
              nuc_inv_mass_W_1pXn0pi_data->Fill(W,inv_mass_1pXn0pi.M());
              nuc_inv_mass_omega_1pXn0pi_data->Fill(q.E(),inv_mass_1pXn0pi.M());
              nuc_inv_mass_q2_1pXn0pi_data->Fill(q2,inv_mass_1pXn0pi.M());
              nuc_inv_mass_xb_1pXn0pi_data->Fill(x_b,inv_mass_1pXn0pi.M());
              nuc_inv_mass_y_1pXn0pi_data->Fill(y,inv_mass_1pXn0pi.M());
              //Run over particle energies
              double proton_totrecoe = 0., pim_totrecoe = 0., pip_totrecoe=0.;
              for(int iPr = 0; iPr < protons.size(); iPr++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  proton_totrecoe = proton_totrecoe + pr.E() - mass_p;
                }
              for(int iPim = 0; iPim < pims.size(); iPim++)
                {
                  SetLorentzVector(pim,pims[iPim]);
                  pim_totrecoe = pim_totrecoe + pim.E();
                }
              for(int iPip = 0; iPip < pips.size(); iPip++)
                {
                  SetLorentzVector(pip,pips[iPip]);
                  pip_totrecoe = pip_totrecoe + pip.E() - mass_pip;
                }
              double totrecoe = proton_totrecoe + pim_totrecoe + pip_totrecoe + el.E();
              //Fill 2D histograms
              totrecoe_W_1pXn0pi_data->Fill(W,totrecoe);
              totrecoe_omega_1pXn0pi_data->Fill(q.E(),totrecoe);
              totrecoe_q2_1pXn0pi_data->Fill(q2,totrecoe);
              totrecoe_xb_1pXn0pi_data->Fill(x_b,totrecoe);
              totrecoe_y_1pXn0pi_data->Fill(y,totrecoe);
              W_q2_1pXn0pi_data->Fill(q2,W);
              W_xb_1pXn0pi_data->Fill(x_b,W);
              W_y_1pXn0pi_data->Fill(y,W);
              q2_xb_1pXn0pi_data->Fill(x_b,q2);
              q2_y_1pXn0pi_data->Fill(y,q2);
              y_xb_1pXn0pi_data->Fill(x_b,y);
            }*/
          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          //1p1pim
          /*if (protons.size()==1 && pips.size()==0 && pims.size()==1)
            {
              //Find invariant mass of the system (help to see if it comes from a true event?)
              TLorentzVector inv_mass_1pXn1pim = pr + pip + pim - q;
              nuc_inv_mass_W_1pXn1pim_data->Fill(W,inv_mass_1pXn1pim.M());
              nuc_inv_mass_omega_1pXn1pim_data->Fill(q.E(),inv_mass_1pXn1pim.M());
              nuc_inv_mass_q2_1pXn1pim_data->Fill(q2,inv_mass_1pXn1pim.M());
              nuc_inv_mass_xb_1pXn1pim_data->Fill(x_b,inv_mass_1pXn1pim.M());
              nuc_inv_mass_y_1pXn1pim_data->Fill(y,inv_mass_1pXn1pim.M());
              //Run over particle energies
              double proton_totrecoe = 0., pim_totrecoe = 0., pip_totrecoe=0.;
              for(int iPr = 0; iPr < protons.size(); iPr++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  proton_totrecoe = proton_totrecoe + pr.E() - mass_p;
                }
              for(int iPim = 0; iPim < pims.size(); iPim++)
                {
                  SetLorentzVector(pim,pims[iPim]);
                  pim_totrecoe = pim_totrecoe + pim.E();
                }
              for(int iPip = 0; iPip < pips.size(); iPip++)
                {
                  SetLorentzVector(pip,pips[iPip]);
                  pip_totrecoe = pip_totrecoe + pip.E() - mass_pip;
                }
              double totrecoe = proton_totrecoe + pim_totrecoe + pip_totrecoe + el.E();
              //Fill 2D histograms
              totrecoe_W_1pXn1pim_data->Fill(W,totrecoe);
              totrecoe_omega_1pXn1pim_data->Fill(q.E(),totrecoe);
              totrecoe_q2_1pXn1pim_data->Fill(q2,totrecoe);
              totrecoe_xb_1pXn1pim_data->Fill(x_b,totrecoe);
              totrecoe_y_1pXn1pim_data->Fill(y,totrecoe);
              W_q2_1pXn1pim_data->Fill(q2,W);
              W_xb_1pXn1pim_data->Fill(x_b,W);
              W_y_1pXn1pim_data->Fill(y,W);
              q2_xb_1pXn1pim_data->Fill(x_b,q2);
              q2_y_1pXn1pim_data->Fill(y,q2);
              y_xb_1pXn1pim_data->Fill(x_b,y);
            }*/
          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          //1D Histograms of particular topologies with cuts
          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          //QE single proton, no pions: only one proton without any charged pions, y:(-0.300,0.300), x_b:(0.8,1.2), W:(0.800,1.050), q2:(0.0,3.0)
          if (protons.size()==1 && pips.size()==0 && pims.size()==0 && y>-0.300 && y<0.300 && x_b>0.8 && x_b<1.2 && W>0.8 && W<1.050 && q2>0. && q2<3.0)
            {
              W_qe_1p_cuts_data->Fill(W);
              xb_qe_1p_cuts_data->Fill(x_b);
              counter_qe_1p0pi_cuts_data++;
              for(int iPr = 0; iPr < protons.size(); iPr++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  double p_proton = pr.P();
                  p_proton_qe_1p_cuts_data->Fill(p_proton);
                }
            }          

          //QE neutron: only one neutron, y:(-300,300), x:(0.8,1.2), W:(0.800,1.050), q2:(0.0,3.0)


          //MEC single proton, no pions
          if (protons.size()==1 && pips.size()==0 && pims.size()==0 && y>0.05 && y<0.300 && x_b>0.8 && x_b<1.2 && W>0.8 && W<1.5 && q2>0. && q2<3.0)
            {
              W_mec_1p_cuts_data->Fill(W);
              xb_mec_1p_cuts_data->Fill(x_b);
              counter_mec_1p0pi_cuts_data++;
              for(int iPr = 0; iPr < protons.size(); iPr++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  double p_proton = pr.P();
                  p_proton_mec_1p_cuts_data->Fill(p_proton);
                }
            }

          //MEC two proton (proton-proton), no pions
          if (protons.size()==2 && pips.size()==0 && pims.size()==0 && y>0.05 && y<0.300 && x_b>0.8 && x_b<1.2 && W>0.8 && W<1.5 && q2>0. && q2<3.0)
            {
              W_mec_2p_cuts_data->Fill(W);
              xb_mec_2p_cuts_data->Fill(x_b);
              counter_mec_2p0pi_cuts_data++;
              for(int iPr = 0; iPr < protons.size(); iPr++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  double p_proton = pr.P();
                  p_proton_mec_2p_cuts_data->Fill(p_proton);
                }
            }

          //MEC proton-neutron


          //RES single charged pion
          if (protons.size()==0 && ((pips.size()==1 && pims.size()==0) || (pips.size()==0 && pims.size()==1)) && y>0.2 && y<0.6 && x_b>0.4 && x_b<0.8 && W>1.050 && W<2.0 && q2>0.8 && q2<6.0)
            {
              W_res_0p1pi_cuts_data->Fill(W);
              xb_res_0p1pi_cuts_data->Fill(x_b);
              counter_res_0p1pi_cuts_data++;
              for(int iPr = 0; iPr < pims.size()+pips.size(); iPr++)
                {
                  if (pips.size()==1 && pims.size()==0)
                    {
                      SetLorentzVector(pip,pips[iPr]);
                      double p_pip = pip.P();
                      p_pi_res_0p1pi_cuts_data->Fill(p_pip);
                    }
                  if (pips.size()==0 && pims.size()==1)
                    {
                      SetLorentzVector(pim,pims[iPr]);
                      double p_pim = pim.P();
                      p_pi_res_0p1pi_cuts_data->Fill(p_pim);
                    }
                }
            }

          //RES single charged pion with single proton (1p1pi)
          if (protons.size()==1 && ((pips.size()==1 && pims.size()==0) || (pips.size()==0 && pims.size()==1)) && y>0.2 && y<0.6 && x_b>0.4 && x_b<0.8 && W>1.050 && W<2.0 && q2>0.8 && q2<6.0)
            {
              W_res_1p1pi_cuts_data->Fill(W);
              xb_res_1p1pi_cuts_data->Fill(x_b);
              counter_res_1p1pi_cuts_data++;
              for(int iPr = 0; iPr < protons.size(); iPr++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  double p_proton = pr.P();
                  p_proton_res_1p1pi_cuts_data->Fill(p_proton);
                }
              for(int iPr = 0; iPr < pims.size()+pips.size(); iPr++)
                {
                  if (pips.size()==1 && pims.size()==0)
                    {
                      SetLorentzVector(pip,pips[iPr]);
                      double p_pip = pip.P();
                      p_pi_res_1p1pi_cuts_data->Fill(p_pip);
                    }
                  if (pips.size()==0 && pims.size()==1)
                    {
                      SetLorentzVector(pim,pims[iPr]);
                      double p_pim = pim.P();
                      p_pi_res_1p1pi_cuts_data->Fill(p_pim);
                    }
                }
            }

          //RES pion-neutron


          //RES pion-pion-neutron

          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

          counter_data++;
	      }
    }

  mult_p_data->SetLineColor(kBlack); mult_p_stacked->Add(mult_p_data); //mult_p->Clear();
  //mult_n_data->SetLineColor(kBlack); mult_n_data->SetLineStyle(1); mult_n_stacked->Add(mult_n_data);
  mult_n_nocuts_data->SetLineColor(kBlack); mult_n_nocuts_data->SetLineStyle(9); mult_n_stacked->Add(mult_n_nocuts_data);
  mult_pip_data->SetLineColor(kBlack); mult_pip_stacked->Add(mult_pip_data); //mult_pip->Clear();
  mult_pim_data->SetLineColor(kBlack); mult_pim_stacked->Add(mult_pim_data); //mult_pim->Clear();
  q2_data->SetLineColor(kBlack); q2_stacked->Add(q2_data); //q2_data->Clear();
  p_chi2_cd_data->SetLineColor(kBlack); p_chi2_cd_data->SetLineStyle(1); p_chi2_stacked->Add(p_chi2_cd_data); //p_chi2_cd_data->Clear();
  n_chi2_cd_data->SetLineColor(kBlack); n_chi2_cd_data->SetLineStyle(1); n_chi2_stacked->Add(n_chi2_cd_data);
  n_chi2_nocuts_data->SetLineColor(kGreen); n_chi2_nocuts_data->SetLineStyle(1); n_chi2_stacked->Add(n_chi2_nocuts_data);
  pip_chi2_cd_data->SetLineColor(kBlack); pip_chi2_cd_data->SetLineStyle(1); pip_chi2_stacked->Add(pip_chi2_cd_data); //pip_chi2_cd_data->Clear();
  pim_chi2_cd_data->SetLineColor(kBlack); pim_chi2_cd_data->SetLineStyle(1); pim_chi2_stacked->Add(pim_chi2_cd_data); //pim_chi2_cd_data->Clear();
  p_chi2_fd_data->SetLineColor(kBlack); p_chi2_fd_data->SetLineStyle(9); p_chi2_stacked->Add(p_chi2_fd_data); //p_chi2_fd_data->Clear();
  n_chi2_fd_data->SetLineColor(kBlack); n_chi2_fd_data->SetLineStyle(9); n_chi2_stacked->Add(n_chi2_fd_data);
  pip_chi2_fd_data->SetLineColor(kBlack); pip_chi2_fd_data->SetLineStyle(9); pip_chi2_stacked->Add(pip_chi2_fd_data); //pip_chi2_fd_data->Clear();
  pim_chi2_fd_data->SetLineColor(kBlack); pim_chi2_fd_data->SetLineStyle(9); pim_chi2_stacked->Add(pim_chi2_fd_data); //pim_chi2_fd_data->Clear();
  energy_transfer_0_data->SetLineColor(kBlack); energy_transfer_0_stacked->Add(energy_transfer_0_data); //energy_transfer_0_data->Clear();
  energy_transfer_1_data->SetLineColor(kBlack); energy_transfer_1_stacked->Add(energy_transfer_1_data); //energy_transfer_1_data->Clear();
  energy_transfer_2_data->SetLineColor(kBlack); energy_transfer_2_stacked->Add(energy_transfer_2_data); //energy_transfer_2_data->Clear();
  energy_transfer_3_data->SetLineColor(kBlack); energy_transfer_3_stacked->Add(energy_transfer_3_data); //energy_transfer_3_data->Clear();
  invariant_mass_data->SetLineColor(kBlack); invariant_mass_stacked->Add(invariant_mass_data); //invariant_mass_data->Clear();
  el_vz_h_data->SetLineColor(kBlack); el_vz_h_stacked->Add(el_vz_h_data); //el_vz_h_data->Clear();
  W_data->SetLineColor(kBlack); W_stacked->Add(W_data); //W_data->Clear();
  theta_data->SetLineColor(kBlack); theta_stacked->Add(theta_data);
  omega_data->SetLineColor(kBlack); omega_stacked->Add(omega_data);
  p_proton_data->SetLineColor(kBlack); p_proton_stacked->Add(p_proton_data);
  p_miss_data->SetLineColor(kBlack); p_miss_stacked->Add(p_miss_data);
  p_miss_theta_data->SetLineColor(kBlack); p_miss_theta_stacked->Add(p_miss_theta_data);
  xb_data->SetLineColor(kBlack); xb_stacked->Add(xb_data);
  y_data->SetLineColor(kBlack); y_stacked->Add(y_data);

  //1D histograms for data for particular channels
  q2_byChannel_data->SetLineColor(1); q2_byChannel_stacked->Add(q2_byChannel_data);
  q2_0p1n0pi_data->SetLineColor(2); q2_0p1n0pi_data->SetMarkerColor(2); q2_0p1n0pi_data->Sumw2(); q2_byChannel_stacked->Add(q2_0p1n0pi_data);
  q2_1pXn0pi_data->SetLineColor(3); q2_1pXn0pi_data->SetMarkerColor(3); q2_1pXn0pi_data->Sumw2(); q2_byChannel_stacked->Add(q2_1pXn0pi_data);
  q2_1pXn1pim_data->SetLineColor(4); q2_byChannel_stacked->Add(q2_1pXn1pim_data);
  q2_0p1n1pip_data->SetLineColor(5); q2_byChannel_stacked->Add(q2_0p1n1pip_data);
  q2_2pXn0pi_data->SetLineColor(6); q2_byChannel_stacked->Add(q2_2pXn0pi_data);
  q2_1p1n0pi_data->SetLineColor(7); q2_byChannel_stacked->Add(q2_1p1n0pi_data);
  q2_1p1n1pi_data->SetLineColor(8); q2_byChannel_stacked->Add(q2_1p1n1pi_data);
  q2_1pXn2pi_data->SetLineColor(9); q2_byChannel_stacked->Add(q2_1pXn2pi_data);
  q2_0p1n2pi_data->SetLineColor(46); q2_byChannel_stacked->Add(q2_0p1n2pi_data);
  totrecoe_0p1n0pi_data->SetLineColor(2); totrecoe_0p1n0pi_data->SetMarkerColor(2); totrecoe_0p1n0pi_data->Sumw2(); totrecoe_byChannel_stacked->Add(totrecoe_0p1n0pi_data);
  totrecoe_1pXn0pi_data->SetLineColor(3); totrecoe_1pXn0pi_data->SetMarkerColor(3); totrecoe_1pXn0pi_data->Sumw2(); totrecoe_byChannel_stacked->Add(totrecoe_1pXn0pi_data);
  totrecoe_1pXn1pim_data->SetLineColor(4); totrecoe_byChannel_stacked->Add(totrecoe_1pXn1pim_data);
  totrecoe_0p1n1pip_data->SetLineColor(5); totrecoe_byChannel_stacked->Add(totrecoe_0p1n1pip_data);
  totrecoe_2pXn0pi_data->SetLineColor(6); totrecoe_byChannel_stacked->Add(totrecoe_2pXn0pi_data);
  totrecoe_1p1n0pi_data->SetLineColor(7); totrecoe_byChannel_stacked->Add(totrecoe_1p1n0pi_data);
  totrecoe_1p1n1pi_data->SetLineColor(8); totrecoe_byChannel_stacked->Add(totrecoe_1p1n1pi_data);
  totrecoe_1pXn2pi_data->SetLineColor(9); totrecoe_byChannel_stacked->Add(totrecoe_1pXn2pi_data);
  totrecoe_0p1n2pi_data->SetLineColor(46); totrecoe_byChannel_stacked->Add(totrecoe_0p1n2pi_data);
  W_byChannel_data->SetLineColor(1); W_byChannel_stacked->Add(W_byChannel_data);
  W_0p1n0pi_data->SetLineColor(2); W_0p1n0pi_data->SetMarkerColor(2); W_0p1n0pi_data->Sumw2(); W_byChannel_stacked->Add(W_0p1n0pi_data);
  W_1pXn0pi_data->SetLineColor(3); W_1pXn0pi_data->SetMarkerColor(3);  W_1pXn0pi_data->Sumw2(); W_byChannel_stacked->Add(W_1pXn0pi_data);
  W_1pXn1pim_data->SetLineColor(4); W_byChannel_stacked->Add(W_1pXn1pim_data);
  W_0p1n1pip_data->SetLineColor(5); W_byChannel_stacked->Add(W_0p1n1pip_data);
  W_2pXn0pi_data->SetLineColor(6); W_byChannel_stacked->Add(W_2pXn0pi_data);
  W_1p1n0pi_data->SetLineColor(7); W_byChannel_stacked->Add(W_1p1n0pi_data);
  W_1p1n1pi_data->SetLineColor(8); W_byChannel_stacked->Add(W_1p1n1pi_data);
  W_1pXn2pi_data->SetLineColor(9); W_byChannel_stacked->Add(W_1pXn2pi_data);
  W_0p1n2pi_data->SetLineColor(46); W_byChannel_stacked->Add(W_0p1n2pi_data);

  //1D Plots by sector
  //W
  W_bySector_data->SetLineColor(1); W_bySector_data_stacked->Add(W_bySector_data);
  W_S1_data->SetLineColor(2); W_S1_data->Sumw2(); W_S1_data->Scale(6.); W_bySector_data_stacked->Add(W_S1_data);
  W_S2_data->SetLineColor(3); W_S2_data->Sumw2(); W_S2_data->Scale(6.); W_bySector_data_stacked->Add(W_S2_data);
  W_S3_data->SetLineColor(4); W_S3_data->Sumw2(); W_S3_data->Scale(6.); W_bySector_data_stacked->Add(W_S3_data);
  W_S4_data->SetLineColor(5); W_S4_data->Sumw2(); W_S4_data->Scale(6.); W_bySector_data_stacked->Add(W_S4_data);
  W_S5_data->SetLineColor(6); W_S5_data->Sumw2(); W_S5_data->Scale(6.); W_bySector_data_stacked->Add(W_S5_data);
  W_S6_data->SetLineColor(7); W_S6_data->Sumw2(); W_S6_data->Scale(6.); W_bySector_data_stacked->Add(W_S6_data);
  //Q^2
  q2_bySector_data->SetLineColor(1); q2_bySector_data_stacked->Add(q2_bySector_data);
  q2_S1_data->SetLineColor(2); q2_S1_data->Sumw2(); q2_S1_data->Scale(6.); q2_bySector_data_stacked->Add(q2_S1_data);
  q2_S2_data->SetLineColor(3); q2_S2_data->Sumw2(); q2_S2_data->Scale(6.); q2_bySector_data_stacked->Add(q2_S2_data);
  q2_S3_data->SetLineColor(4); q2_S3_data->Sumw2(); q2_S3_data->Scale(6.); q2_bySector_data_stacked->Add(q2_S3_data);
  q2_S4_data->SetLineColor(5); q2_S4_data->Sumw2(); q2_S4_data->Scale(6.); q2_bySector_data_stacked->Add(q2_S4_data);
  q2_S5_data->SetLineColor(6); q2_S5_data->Sumw2(); q2_S5_data->Scale(6.); q2_bySector_data_stacked->Add(q2_S5_data);
  q2_S6_data->SetLineColor(7); q2_S6_data->Sumw2(); q2_S6_data->Scale(6.); q2_bySector_data_stacked->Add(q2_S6_data);

  //Print out integrals of particular channels
  outfile_byChannel << "Total counts for run 15672 of Ar-40 @ 2.1 GeV : " << q2_byChannel_data->Integral() << " per integrated luminosity: 1.35 fb^-1" << std::endl;
  outfile_byChannel << "Total expected counts for Argon-40 : " << (q2_byChannel_data->Integral())*(1.35/0.04) << " per integrated luminosity: 1.35 fb^-1" << std::endl;
  outfile_byChannel << "Total expected counts for Argon-40 per 1 fb^-1 : " << (q2_byChannel_data->Integral())/0.04 << std::endl;
  outfile_byChannel << "Channel    |    Ar40 #15429   | Exp. Ar40 Total |   Ar40 per 1 fb^-1"   << std::endl;
  outfile_byChannel << "----------------------------------------------------------------------" << std::endl;
  outfile_byChannel << "0p1n0pi    |    " << q2_0p1n0pi_data->Integral()  << "   |   " << (q2_0p1n0pi_data->Integral())*(1.35/0.04)  << "   |   " << (q2_0p1n0pi_data->Integral())/0.04  << std::endl;
  outfile_byChannel << "1pXn0pi    |    " << q2_1pXn0pi_data->Integral()  << "   |   " << (q2_1pXn0pi_data->Integral())*(1.35/0.04)  << "   |   " << (q2_1pXn0pi_data->Integral())/0.04  << std::endl;
  outfile_byChannel << "1pXn1pim   |    " << q2_1pXn1pim_data->Integral() << "   |   " << (q2_1pXn1pim_data->Integral())*(1.35/0.04) << "   |   " << (q2_1pXn1pim_data->Integral())/0.04 << std::endl;
  outfile_byChannel << "0p1n1pip   |    " << q2_0p1n1pip_data->Integral() << "   |   " << (q2_0p1n1pip_data->Integral())*(1.35/0.04) << "   |   " << (q2_0p1n1pip_data->Integral())/0.04 << std::endl;
  outfile_byChannel << "2pXn0pi    |    " << q2_2pXn0pi_data->Integral()  << "   |   " << (q2_2pXn0pi_data->Integral())*(1.35/0.04)  << "   |   " << (q2_2pXn0pi_data->Integral())/0.04  << std::endl;
  outfile_byChannel << "1p1n0pi    |    " << q2_1p1n0pi_data->Integral()  << "   |   " << (q2_1p1n0pi_data->Integral())*(1.35/0.04)  << "   |   " << (q2_1p1n0pi_data->Integral())/0.04  << std::endl;
  outfile_byChannel << "1p1n1pi    |    " << q2_1p1n1pi_data->Integral()  << "   |   " << (q2_1p1n1pi_data->Integral())*(1.35/0.04)  << "   |   " << (q2_1p1n1pi_data->Integral())/0.04  << std::endl;
  outfile_byChannel << "1pXn2pi    |    " << q2_1pXn2pi_data->Integral()  << "   |   " << (q2_1pXn2pi_data->Integral())*(1.35/0.04)  << "   |   " << (q2_1pXn2pi_data->Integral())/0.04  << std::endl;
  outfile_byChannel << "0p1n2pi    |    " << q2_0p1n2pi_data->Integral()  << "   |   " << (q2_0p1n2pi_data->Integral())*(1.35/0.04)  << "   |   " << (q2_0p1n2pi_data->Integral())/0.04  << std::endl;

  //Histograms with kinematics cuts for particular topologies
  //Single proton QE
  W_qe_1p_cuts_data->SetLineColor(kBlack); W_qe_1p_cuts_stacked->Add(W_qe_1p_cuts_data);
  xb_qe_1p_cuts_data->SetLineColor(kBlack); xb_qe_1p_cuts_stacked->Add(xb_qe_1p_cuts_data);
  p_proton_qe_1p_cuts_data->SetLineColor(kBlack); p_proton_qe_1p_cuts_stacked->Add(p_proton_qe_1p_cuts_data);
  //Single proton MEC
  W_mec_1p_cuts_data->SetLineColor(kBlack); W_mec_1p_cuts_stacked->Add(W_mec_1p_cuts_data);
  xb_mec_1p_cuts_data->SetLineColor(kBlack); xb_mec_1p_cuts_stacked->Add(xb_mec_1p_cuts_data);
  p_proton_mec_1p_cuts_data->SetLineColor(kBlack); p_proton_mec_1p_cuts_stacked->Add(p_proton_mec_1p_cuts_data);
  //Two proton (proton-proton) MEC
  W_mec_2p_cuts_data->SetLineColor(kBlack); W_mec_2p_cuts_stacked->Add(W_mec_2p_cuts_data);
  xb_mec_2p_cuts_data->SetLineColor(kBlack); xb_mec_2p_cuts_stacked->Add(xb_mec_2p_cuts_data);
  p_proton_mec_2p_cuts_data->SetLineColor(kBlack); p_proton_mec_2p_cuts_stacked->Add(p_proton_mec_2p_cuts_data);
  //RES single charged pion (0p1pi)
  W_res_0p1pi_cuts_data->SetLineColor(kBlack); W_res_0p1pi_cuts_stacked->Add(W_res_0p1pi_cuts_data);
  xb_res_0p1pi_cuts_data->SetLineColor(kBlack); xb_res_0p1pi_cuts_stacked->Add(xb_res_0p1pi_cuts_data);
  p_pi_res_0p1pi_cuts_data->SetLineColor(kBlack); p_pi_res_0p1pi_cuts_stacked->Add(p_pi_res_0p1pi_cuts_data);
  //RES single charged pion with single proton (1p1pi)
  W_res_1p1pi_cuts_data->SetLineColor(kBlack); W_res_1p1pi_cuts_stacked->Add(W_res_1p1pi_cuts_data);
  xb_res_1p1pi_cuts_data->SetLineColor(kBlack); xb_res_1p1pi_cuts_stacked->Add(xb_res_1p1pi_cuts_data);
  p_pi_res_1p1pi_cuts_data->SetLineColor(kBlack); p_pi_res_1p1pi_cuts_stacked->Add(p_pi_res_1p1pi_cuts_data);
  p_proton_res_1p1pi_cuts_data->SetLineColor(kBlack); p_proton_res_1p1pi_cuts_stacked->Add(p_proton_res_1p1pi_cuts_data);

  //////////////////////////////////////////////////////////////////////////////
  //         Monte Carlo (GENIE)
  //////////////////////////////////////////////////////////////////////////////
  while(chain_MC.Next())
    {
      beam_energy_MC = beamE;
      beam_MC.SetE(beam_energy_MC);
      beam_MC.SetPz(beam_energy_MC);

      if( beamE > 1e-4)
	      {
	        if( counter_MC == 0)
	        cout<<"Beam energy set manually: "<<beamE<<endl;
	        beam_MC.SetE(beamE);
	        beam_MC.SetPz(beamE);
	      }

      //TLorentzVector el_mc;
      c12_MC->mcparts()->setEntry(0);

      //el_mc.SetXYZM(c12_MC->mcparts()->getPx(), c12_MC->mcparts()->getPy(),c12_MC->mcparts()->getPz(),db->GetParticle(11)->Mass());

      // get particles by type
      auto electrons=c12_MC->getByID(11);
      auto protons=c12_MC->getByID(2212);
      auto neutrons=c12_MC->getByID(2112);
      auto pips=c12_MC->getByID(211);
      auto pims=c12_MC->getByID(-211);

      //Make sure that we have a good electron and that we are not seeing the weird delta function around zero in MC
      if(electrons.size()==1 && (electrons[0]->par()->getVz() < -0.00001 || electrons[0]->par()->getVz() > 0.00001) )
	      {
          double energy =  (electrons[0]->cal(PCAL)->getEnergy() +  electrons[0]->cal(ECIN)->getEnergy() +  electrons[0]->cal(ECOUT)->getEnergy())/electrons[0]->getP();
          
          bool ecal = ( (electrons[0]->cal(ECIN)->getLv() >= 14 && electrons[0]->cal(ECIN)->getLw() >= 14) && (energy > 0.18 && energy < 0.28) && electrons[0]->getP() > 0.5 && electrons[0]->getP() < 10 );

          SetLorentzVector(el,electrons[0]);

          TLorentzVector q = beam_MC - el; //photon 4-vector
          //Cut out all high radiative correction events
          if (q.E()>omega_threshold_MC){continue;}

          //Calculate some important kinematic variables
          double q2        = -q.M2();
          double x_b       = q2/(2 * mass_p * q.E() ); //Bjorken-x 
          double y         = -q.P() + sqrt( q.E()*q.E() + 2*q.E()*mass_p); //y (Bjorken-y), the scaling variable; from Jourdan, eq. 20
          double vz_e      = electrons[0]->par()->getVz();
          double W         = sqrt(mass_p*mass_p - q2 + 2*q.E()*mass_p); //Invariant mass, assuming electron off of a standing proton

          double processid = c12_MC->mcevent()->getWeight();//code = 1,2,3,4 = type = qel,mec,res,dis
          double resid = c12_MC->mcevent()->getPtarget();//targP = target polarization, which is not relevant = resid

          //std::cout << vz_e << std::endl;
          el_vz_h_MC->Fill(vz_e);
          
          //Electron quality cuts
          bool targ_1 = abs( vz_e - (-1.48)) < 1 * 1.25;
          bool targ_2 = abs( vz_e - (-6.3))  < 1 * 1.25;

          //	  ecal_MC->Fill(electrons[0]->getP(),energy/electrons[0]->getP());
              /////////Modified
          //if( !(ecal && (targ_1 || targ_2) ) )//empty rgk run

          //if( !(ecal && vz_e < vz_max && vz_e > vz_min))
          ////////Modified
          //continue;

          //el_vz_h->Fill(vz_e);
          ecal_MC->Fill(electrons[0]->getP(),energy);

          //remove FT (seems to be still in the MC; FT is not used in RGM)
          //if( electrons[0]->getRegion() == 1000)
            //continue;

          //Electron kinematics monitoring
          //(e,e') cross sections

          //Counters for various simulation processes
          if (processid==1){counter_qe_MC++;}
          if (processid==2){counter_mec_MC++;}
          if (processid==3){counter_res_MC++;}
          if (processid==4){counter_dis_MC++;}

          //Fill various histograms
          theta_MC->Fill(el.Theta()*TMath::RadToDeg());
          if (processid==1){theta_qe_MC->Fill(el.Theta()*TMath::RadToDeg());}
          if (processid==2){theta_mec_MC->Fill(el.Theta()*TMath::RadToDeg());}
          if (processid==3){theta_res_MC->Fill(el.Theta()*TMath::RadToDeg());}
          if (processid==4){theta_dis_MC->Fill(el.Theta()*TMath::RadToDeg());}
          theta_mom_MC->Fill(el.P(),el.Theta()*TMath::RadToDeg());
          q2_MC->Fill(q2);
          if (processid==1){q2_qe_MC->Fill(q2);}
          if (processid==2){q2_mec_MC->Fill(q2);}
          if (processid==3){q2_res_MC->Fill(q2);}
          if (processid==4){q2_dis_MC->Fill(q2);}
          q2_theta_MC->Fill(q2,el.Theta()*TMath::RadToDeg());
          W_MC->Fill(W);
          if (processid==1){W_qe_MC->Fill(W);}
          if (processid==2){W_mec_MC->Fill(W);}
          if (processid==3){W_res_MC->Fill(W);}
          if (processid==4){W_dis_MC->Fill(W);}    
          xb_MC->Fill(x_b);
          if (processid==1){xb_qe_MC->Fill(x_b);}
          if (processid==2){xb_mec_MC->Fill(x_b);}
          if (processid==3){xb_res_MC->Fill(x_b);}
          if (processid==4){xb_dis_MC->Fill(x_b);}
          y_MC->Fill(y);
          if (processid==1){y_qe_MC->Fill(y);}
          if (processid==2){y_mec_MC->Fill(y);}
          if (processid==3){y_res_MC->Fill(y);}
          if (processid==4){y_dis_MC->Fill(y);}
          q2_W_MC->Fill(W,q2);
          omega_MC->Fill(q.E());
          if (processid==1){omega_qe_MC->Fill(q.E());}
          if (processid==2){omega_mec_MC->Fill(q.E());}
          if (processid==3){omega_res_MC->Fill(q.E());}
          if (processid==4){omega_dis_MC->Fill(q.E());}
          omega_theta_MC->Fill(q.E(),el.Theta()*TMath::RadToDeg());
          W_theta_MC->Fill(W,el.Theta()*TMath::RadToDeg());
          q0_vs_q3_MC->Fill(q.Pz(),q.E());
          if (protons.size()==1 && pips.size()==0 && pims.size()==0){q0_vs_q3_1pXn0pi_MC->Fill(q.Pz(),q.E());}
          if (protons.size()==2 && pips.size()==0 && pims.size()==0){q0_vs_q3_2pXn0pi_MC->Fill(q.Pz(),q.E());}

          processid_MC->Fill(processid);
          if (processid==3 && resid>=0){resid_MC->Fill(resid);}

          //By sectors
          W_bySector_MC->Fill(W);
          q2_bySector_MC->Fill(q2);
          //Sector 1
          if (electrons[0]->trk(DC)->getSector() == 1)
            {
              W_S1_MC->Fill(W);
              q2_S1_MC->Fill(q2);
              W_theta_S1_MC->Fill(W,el.Theta()*TMath::RadToDeg());
            }
          //Etc...
          if (electrons[0]->trk(DC)->getSector() == 2)
            {
              W_S2_MC->Fill(W);
              q2_S2_MC->Fill(q2);
              W_theta_S2_MC->Fill(W,el.Theta()*TMath::RadToDeg());
            }
          if (electrons[0]->trk(DC)->getSector() == 3)
            {
              W_S3_MC->Fill(W);
              q2_S3_MC->Fill(q2);
              W_theta_S3_MC->Fill(W,el.Theta()*TMath::RadToDeg());
            }
          if (electrons[0]->trk(DC)->getSector() == 4)
            {
              W_S4_MC->Fill(W);
              q2_S4_MC->Fill(q2);
              W_theta_S4_MC->Fill(W,el.Theta()*TMath::RadToDeg());
            }
          if (electrons[0]->trk(DC)->getSector() == 5)
            {
              W_S5_MC->Fill(W);
              q2_S5_MC->Fill(q2);
              W_theta_S5_MC->Fill(W,el.Theta()*TMath::RadToDeg());
            }
          if (electrons[0]->trk(DC)->getSector() == 6)
            {
              W_S6_MC->Fill(W);
              q2_S6_MC->Fill(q2);
              W_theta_S6_MC->Fill(W,el.Theta()*TMath::RadToDeg());
            }

          //Find resonances in the MC (make sure to cut out all the -99s) with W between 1.2 and 2.1 GeV/c^2
          if (W > 1.2 && W < 2.1 && processid==3 && resid>=0)
            {
              W_resonances_1p2_to_2p1_MC->Fill(W);
              resonances_1p2_to_2p1_MC->Fill(resid);
              W_resid_1p2_to_2p1_MC->Fill(W,resid);
              counter_1p2_to_2p1res_MC++;
            }

          if(el.Theta()*TMath::RadToDeg() < 10 && el.Theta()*TMath::RadToDeg() > 5)
            {
              energy_transfer_0_MC->Fill(q.E());
              if (processid==1){energy_transfer_qe_0_MC->Fill(q.E());}
              if (processid==2){energy_transfer_mec_0_MC->Fill(q.E());}
              if (processid==3){energy_transfer_res_0_MC->Fill(q.E());}
              if (processid==4){energy_transfer_dis_0_MC->Fill(q.E());}              
            }
          else if(el.Theta()*TMath::RadToDeg() < 15  && el.Theta()*TMath::RadToDeg() >= 10)
            {
              energy_transfer_1_MC->Fill(q.E());
              if (processid==1){energy_transfer_qe_1_MC->Fill(q.E());}
              if (processid==2){energy_transfer_mec_1_MC->Fill(q.E());}
              if (processid==3){energy_transfer_res_1_MC->Fill(q.E());}
              if (processid==4){energy_transfer_dis_1_MC->Fill(q.E());}
            }
          else if(el.Theta()*TMath::RadToDeg() < 25  && el.Theta()*TMath::RadToDeg() >= 15)
            {
              energy_transfer_2_MC->Fill(q.E());
              if (processid==1){energy_transfer_qe_2_MC->Fill(q.E());}
              if (processid==2){energy_transfer_mec_2_MC->Fill(q.E());}
              if (processid==3){energy_transfer_res_2_MC->Fill(q.E());}
              if (processid==4){energy_transfer_dis_2_MC->Fill(q.E());}
            }
          else if(el.Theta()*TMath::RadToDeg() < 35  && el.Theta()*TMath::RadToDeg() >= 25)
            {
              energy_transfer_3_MC->Fill(q.E());
              if (processid==1){energy_transfer_qe_3_MC->Fill(q.E());}
              if (processid==2){energy_transfer_mec_3_MC->Fill(q.E());}
              if (processid==3){energy_transfer_res_3_MC->Fill(q.E());}
              if (processid==4){energy_transfer_dis_3_MC->Fill(q.E());}
            }
          
          //By sectors
          W_bySector_MC->Fill(W);
          q2_bySector_MC->Fill(q2);
          //Sector 1
          if (electrons[0]->trk(DC)->getSector() == 1)
            {
              W_S1_MC->Fill(W);
              q2_S1_MC->Fill(q2);
            }
          //Etc...
          if (electrons[0]->trk(DC)->getSector() == 2)
            {
              W_S2_MC->Fill(W);
              q2_S2_MC->Fill(q2);
            }
          if (electrons[0]->trk(DC)->getSector() == 3)
            {
              W_S3_MC->Fill(W);
              q2_S3_MC->Fill(q2);
            }
          if (electrons[0]->trk(DC)->getSector() == 4)
            {
              W_S4_MC->Fill(W);
              q2_S4_MC->Fill(q2);
            }
          if (electrons[0]->trk(DC)->getSector() == 5)
            {
              W_S5_MC->Fill(W);
              q2_S5_MC->Fill(q2);
            }
          if (electrons[0]->trk(DC)->getSector() == 6)
            {
              W_S6_MC->Fill(W);
              q2_S6_MC->Fill(q2);
            }
          
          //MULTIPLICITIES 
          //Monitor the proton, neutron, and pion multiplicities 
          int num_p = 0;
          int num_n = 0;
          int num_n_nocuts = 0;
          int num_pip = 0;
          int num_pim = 0;
          
          for(int iPr = 0; iPr < protons.size(); iPr++)
            {
              double p_chi2_cd_v = 999999;
              double p_chi2_fd_v = 999999;
              
              if(protons[iPr]->getRegion()==CD)
                {
                  p_chi2_cd_MC->Fill(protons[iPr]->par()->getChi2Pid());
                  p_chi2_cd_v = protons[iPr]->par()->getChi2Pid();
                }
            
              if(protons[iPr]->getRegion()==FD)
                {
                  p_chi2_fd_MC->Fill(protons[iPr]->par()->getChi2Pid());
                  p_chi2_fd_v = protons[iPr]->par()->getChi2Pid();
                }
              
              if( ( ((p_chi2_fd_v < p_chi2_fd_cut_hi) && (p_chi2_fd_v > p_chi2_fd_cut_lo)) || ((p_chi2_cd_v < p_chi2_cd_cut_hi) && (p_chi2_cd_v > p_chi2_cd_cut_lo)) ) )
                {
                  num_p++;
                }
            }

          for(int iNr = 0; iNr < neutrons.size(); iNr++)
            {
              //num_n_nocuts++;
              n_chi2_nocuts_MC->Fill(neutrons[iNr]->par()->getChi2Pid());

              double n_chi2_cd_v = 999999;
              double n_chi2_fd_v = 999999;
              
              if(neutrons[iNr]->getRegion()==CD)
                {
                  n_chi2_cd_MC->Fill(neutrons[iNr]->par()->getChi2Pid());
                  n_chi2_cd_v = neutrons[iNr]->par()->getChi2Pid();
                }
            
              if(neutrons[iNr]->getRegion()==FD)
                {
                  n_chi2_fd_MC->Fill(neutrons[iNr]->par()->getChi2Pid());
                  n_chi2_fd_v = neutrons[iNr]->par()->getChi2Pid();
                }
              
              if( ( ((n_chi2_fd_v < n_chi2_fd_cut_hi) && (n_chi2_fd_v > n_chi2_fd_cut_lo)) || ((n_chi2_cd_v < n_chi2_cd_cut_hi) && (n_chi2_cd_v > n_chi2_cd_cut_lo)) ) )
                {
                  num_n++;
                }
            }
	  
          for(int iPim = 0; iPim < pims.size(); iPim++)
            {
              double pim_chi2_cd_v = 999999;
              double pim_chi2_fd_v = 999999;
              
              if(pims[iPim]->getRegion()==CD)
                {
                  pim_chi2_cd_MC->Fill(pims[iPim]->par()->getChi2Pid());
                  pim_chi2_cd_v = pims[iPim]->par()->getChi2Pid();
                }
            
              if(pims[iPim]->getRegion()==FD)
                {
                  pim_chi2_fd_MC->Fill(pims[iPim]->par()->getChi2Pid());
                  pim_chi2_fd_v = pims[iPim]->par()->getChi2Pid();
                }
              
              if( ( ((pim_chi2_fd_v < pim_chi2_fd_cut_hi) && (pim_chi2_fd_v > pim_chi2_fd_cut_lo)) || ((pim_chi2_cd_v < pim_chi2_cd_cut_hi) && (pim_chi2_cd_v > pim_chi2_cd_cut_lo)) ) )
                {
                  num_pim++;
                }
            }

          for(int iPip = 0; iPip < pips.size(); iPip++)
            {
              double pip_chi2_cd_v = 999999;
              double pip_chi2_fd_v = 999999;
              
              if(pips[iPip]->getRegion()==CD)
                {
                  pip_chi2_cd_MC->Fill(pips[iPip]->par()->getChi2Pid());
                  pip_chi2_cd_v = pips[iPip]->par()->getChi2Pid();
                }           
            
              if(pips[iPip]->getRegion()==FD)
                {
                  pip_chi2_fd_MC->Fill(pips[iPip]->par()->getChi2Pid());
                  pip_chi2_fd_v = pips[iPip]->par()->getChi2Pid();
                }
            
              if( ( ((pip_chi2_fd_v < pip_chi2_fd_cut_hi) && (pip_chi2_fd_v > pip_chi2_fd_cut_lo)) || ((pip_chi2_cd_v < pip_chi2_cd_cut_hi) && (pip_chi2_cd_v > pip_chi2_cd_cut_lo)) ) )
                {
                  num_pip++;
                }
            }
        
          mult_p_MC->Fill(num_p);
          mult_n_MC->Fill(num_n);
          mult_n_nocuts_MC->Fill(neutrons.size());
          mult_pim_MC->Fill(num_pim);
          mult_pip_MC->Fill(num_pip);
          mult_p_pis_MC->Fill(num_p,num_pim+num_pip);
          mult_p_n_MC->Fill(num_p,neutrons.size());
          
	  
          //INVARIANT MASS
          //Calculating the Invariant mass of the p + pi- system 
          for(int iPr = 0; iPr < protons.size(); iPr++)
            {
              for(int iPim = 0; iPim < pims.size(); iPim++)
                {	
                  SetLorentzVector(pr,protons[iPr]);
                  SetLorentzVector(pim,pims[iPim]);
                  
                  double vz_e_p = electrons[0]->par()->getVz() - protons[iPr]->par()->getVz();
                  TLorentzVector inv_mass = pr + pim; //missing 4-vector
                  
                  invariant_mass_MC->Fill(inv_mass.M());
                }
	          } 

          //Proton variables
          for(int iPr = 0; iPr < protons.size(); iPr++)
            {
              SetLorentzVector(pr,protons[iPr]);
              double p_proton = pr.P();
              //Get p_miss vector stuff
              TLorentzVector p_miss = pr - q;//4-vector subtraction
              double p_miss_event = p_miss.P();
              double p_miss_theta = p_miss.Theta()*TMath::RadToDeg();
              //Fill histograms
              p_proton_MC->Fill(p_proton);
              p_miss_MC->Fill(p_miss_event);
              p_miss_theta_MC->Fill(p_miss_theta);
              p_miss_theta_vs_omega_MC->Fill(q.E(),p_miss_theta);
            }

          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          //2D Histograms of particular topologies without cuts (what should the cuts be)
          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          //0p1pip
          if (protons.size()==0 && pips.size()==1 && pims.size()==0)
            {
              //Find invariant mass of the system (help to see if it comes from a true event?)
              TLorentzVector inv_mass_0p_1pip = pr + pip + pim - q;
              nuc_inv_mass_W_0p_1pip_MC->Fill(W,inv_mass_0p_1pip.M());
              nuc_inv_mass_omega_0p_1pip_MC->Fill(q.E(),inv_mass_0p_1pip.M());
              nuc_inv_mass_q2_0p_1pip_MC->Fill(q2,inv_mass_0p_1pip.M());
              nuc_inv_mass_xb_0p_1pip_MC->Fill(x_b,inv_mass_0p_1pip.M());
              nuc_inv_mass_y_0p_1pip_MC->Fill(y,inv_mass_0p_1pip.M());
              //Run over particle energies
              double proton_totrecoe = 0., pim_totrecoe = 0., pip_totrecoe=0.;
              for(int iPr = 0; iPr < protons.size(); iPr++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  proton_totrecoe = proton_totrecoe + pr.E() - mass_p;
                }
              for(int iPim = 0; iPim < pims.size(); iPim++)
                {
                  SetLorentzVector(pim,pims[iPim]);
                  pim_totrecoe = pim_totrecoe + pim.E();
                }
              for(int iPip = 0; iPip < pips.size(); iPip++)
                {
                  SetLorentzVector(pip,pips[iPip]);
                  pip_totrecoe = pip_totrecoe + pip.E() - mass_pip;
                }
              double totrecoe = proton_totrecoe + pim_totrecoe + pip_totrecoe + el.E();
              //Fill 2D histograms
              totrecoe_W_0p_1pip_MC->Fill(W,totrecoe);
              totrecoe_omega_0p_1pip_MC->Fill(q.E(),totrecoe);
              totrecoe_q2_0p_1pip_MC->Fill(q2,totrecoe);
              totrecoe_xb_0p_1pip_MC->Fill(x_b,totrecoe);
              totrecoe_y_0p_1pip_MC->Fill(y,totrecoe);
              W_q2_0p_1pip_MC->Fill(q2,W);
              W_xb_0p_1pip_MC->Fill(x_b,W);
              W_y_0p_1pip_MC->Fill(y,W);
              q2_xb_0p_1pip_MC->Fill(x_b,q2);
              q2_y_0p_1pip_MC->Fill(y,q2);
              y_xb_0p_1pip_MC->Fill(x_b,y);
            }
          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          //0p1pim
          if (protons.size()==0 && pips.size()==0 && pims.size()==1)
            {
              //Find invariant mass of the system (help to see if it comes from a true event?)
              TLorentzVector inv_mass_0p_1pim = pr + pip + pim - q;
              nuc_inv_mass_W_0p_1pim_MC->Fill(W,inv_mass_0p_1pim.M());
              nuc_inv_mass_omega_0p_1pim_MC->Fill(q.E(),inv_mass_0p_1pim.M());
              nuc_inv_mass_q2_0p_1pim_MC->Fill(q2,inv_mass_0p_1pim.M());
              nuc_inv_mass_xb_0p_1pim_MC->Fill(x_b,inv_mass_0p_1pim.M());
              nuc_inv_mass_y_0p_1pim_MC->Fill(y,inv_mass_0p_1pim.M());
              //Run over particle energies
              double proton_totrecoe = 0., pim_totrecoe = 0., pip_totrecoe=0.;
              for(int iPr = 0; iPr < protons.size(); iPr++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  proton_totrecoe = proton_totrecoe + pr.E() - mass_p;
                }
              for(int iPim = 0; iPim < pims.size(); iPim++)
                {
                  SetLorentzVector(pim,pims[iPim]);
                  pim_totrecoe = pim_totrecoe + pim.E();
                }
              for(int iPip = 0; iPip < pips.size(); iPip++)
                {
                  SetLorentzVector(pip,pips[iPip]);
                  pip_totrecoe = pip_totrecoe + pip.E() - mass_pip;
                }
              double totrecoe = proton_totrecoe + pim_totrecoe + pip_totrecoe + el.E();
              //Fill 2D histograms
              totrecoe_W_0p_1pim_MC->Fill(W,totrecoe);
              totrecoe_omega_0p_1pim_MC->Fill(q.E(),totrecoe);
              totrecoe_q2_0p_1pim_MC->Fill(q2,totrecoe);
              totrecoe_xb_0p_1pim_MC->Fill(x_b,totrecoe);
              totrecoe_y_0p_1pim_MC->Fill(y,totrecoe);
              W_q2_0p_1pim_MC->Fill(q2,W);
              W_xb_0p_1pim_MC->Fill(x_b,W);
              W_y_0p_1pim_MC->Fill(y,W);
              q2_xb_0p_1pim_MC->Fill(x_b,q2);
              q2_y_0p_1pim_MC->Fill(y,q2);
              y_xb_0p_1pim_MC->Fill(x_b,y);
            }
          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          //1p0pi
          /*if (protons.size()==1 && pips.size()==0 && pims.size()==0)
            {
              //Find invariant mass of the system (help to see if it comes from a true event?)
              TLorentzVector inv_mass_1pXn0pi = pr + pip + pim - q;
              nuc_inv_mass_W_1pXn0pi_MC->Fill(W,inv_mass_1pXn0pi.M());
              nuc_inv_mass_omega_1pXn0pi_MC->Fill(q.E(),inv_mass_1pXn0pi.M());
              nuc_inv_mass_q2_1pXn0pi_MC->Fill(q2,inv_mass_1pXn0pi.M());
              nuc_inv_mass_xb_1pXn0pi_MC->Fill(x_b,inv_mass_1pXn0pi.M());
              nuc_inv_mass_y_1pXn0pi_MC->Fill(y,inv_mass_1pXn0pi.M());
              //Run over particle energies
              double proton_totrecoe = 0., pim_totrecoe = 0., pip_totrecoe=0.;
              for(int iPr = 0; iPr < protons.size(); iPr++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  proton_totrecoe = proton_totrecoe + pr.E() - mass_p;
                }
              for(int iPim = 0; iPim < pims.size(); iPim++)
                {
                  SetLorentzVector(pim,pims[iPim]);
                  pim_totrecoe = pim_totrecoe + pim.E();
                }
              for(int iPip = 0; iPip < pips.size(); iPip++)
                {
                  SetLorentzVector(pip,pips[iPip]);
                  pip_totrecoe = pip_totrecoe + pip.E() - mass_pip;
                }
              double totrecoe = proton_totrecoe + pim_totrecoe + pip_totrecoe + el.E();
              //Fill 2D histograms
              totrecoe_W_1pXn0pi_MC->Fill(W,totrecoe);
              totrecoe_omega_1pXn0pi_MC->Fill(q.E(),totrecoe);
              totrecoe_q2_1pXn0pi_MC->Fill(q2,totrecoe);
              totrecoe_xb_1pXn0pi_MC->Fill(x_b,totrecoe);
              totrecoe_y_1pXn0pi_MC->Fill(y,totrecoe);
              W_q2_1pXn0pi_MC->Fill(q2,W);
              W_xb_1pXn0pi_MC->Fill(x_b,W);
              W_y_1pXn0pi_MC->Fill(y,W);
              q2_xb_1pXn0pi_MC->Fill(x_b,q2);
              q2_y_1pXn0pi_MC->Fill(y,q2);
              y_xb_1pXn0pi_MC->Fill(x_b,y);
            }*/
          //1pXn0pi--no neutron/neutral particle cuts
          if ( 
               protons.size()==1
            && (   //Must check that each proton has a good chi2 value range in either the forward or central detectors
                   ((protons[0]->getRegion()==FD) && (protons[0]->par()->getChi2Pid() < p_chi2_fd_cut_hi) && (protons[0]->par()->getChi2Pid() > p_chi2_fd_cut_lo))
                || ((protons[0]->getRegion()==CD) && (protons[0]->par()->getChi2Pid() < p_chi2_cd_cut_hi) && (protons[0]->par()->getChi2Pid() > p_chi2_cd_cut_lo))
               )
            && pips.size()==0
            && pims.size()==0 
             )
            {
              double proton_totrecoe = 0., neutron_totrecoe = 0., pim_totrecoe = 0., pip_totrecoe=0.;
              for(int iPr = 0; iPr < protons.size(); iPr++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  proton_totrecoe = proton_totrecoe + pr.E() - mass_p;
                }
              double totrecoe = proton_totrecoe + neutron_totrecoe + pim_totrecoe + pip_totrecoe + el.E();
              //q2_1pXn0pi_MC->Fill(q2);
              //totrecoe_1pXn0pi_MC->Fill(totrecoe);
              //W_1pXn0pi_MC->Fill(W);
              //Find invariant mass of the system (help to see if it comes from a true event?)
              TLorentzVector inv_mass_1pXn0pi = pr + pip + pim - q;
              nuc_inv_mass_W_1pXn0pi_MC->Fill(W,inv_mass_1pXn0pi.M());
              nuc_inv_mass_omega_1pXn0pi_MC->Fill(q.E(),inv_mass_1pXn0pi.M());
              nuc_inv_mass_q2_1pXn0pi_MC->Fill(q2,inv_mass_1pXn0pi.M());
              nuc_inv_mass_xb_1pXn0pi_MC->Fill(x_b,inv_mass_1pXn0pi.M());
              nuc_inv_mass_y_1pXn0pi_MC->Fill(y,inv_mass_1pXn0pi.M());
              //Fill 2D histograms
              totrecoe_W_1pXn0pi_MC->Fill(W,totrecoe);
              totrecoe_omega_1pXn0pi_MC->Fill(q.E(),totrecoe);
              totrecoe_q2_1pXn0pi_MC->Fill(q2,totrecoe);
              totrecoe_xb_1pXn0pi_MC->Fill(x_b,totrecoe);
              totrecoe_y_1pXn0pi_MC->Fill(y,totrecoe);
              W_q2_1pXn0pi_MC->Fill(q2,W);
              W_xb_1pXn0pi_MC->Fill(x_b,W);
              W_y_1pXn0pi_MC->Fill(y,W);
              q2_xb_1pXn0pi_MC->Fill(x_b,q2);
              q2_y_1pXn0pi_MC->Fill(y,q2);
              y_xb_1pXn0pi_MC->Fill(x_b,y);
            }
          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          //1p1pim
          /*if (protons.size()==1 && pips.size()==0 && pims.size()==1)
            {
              //Find invariant mass of the system (help to see if it comes from a true event?)
              TLorentzVector inv_mass_1pXn1pim = pr + pip + pim - q;
              nuc_inv_mass_W_1pXn1pim_MC->Fill(W,inv_mass_1pXn1pim.M());
              nuc_inv_mass_omega_1pXn1pim_MC->Fill(q.E(),inv_mass_1pXn1pim.M());
              nuc_inv_mass_q2_1pXn1pim_MC->Fill(q2,inv_mass_1pXn1pim.M());
              nuc_inv_mass_xb_1pXn1pim_MC->Fill(x_b,inv_mass_1pXn1pim.M());
              nuc_inv_mass_y_1pXn1pim_MC->Fill(y,inv_mass_1pXn1pim.M());
              //Run over particle energies
              double proton_totrecoe = 0., pim_totrecoe = 0., pip_totrecoe=0.;
              for(int iPr = 0; iPr < protons.size(); iPr++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  proton_totrecoe = proton_totrecoe + pr.E() - mass_p;
                }
              for(int iPim = 0; iPim < pims.size(); iPim++)
                {
                  SetLorentzVector(pim,pims[iPim]);
                  pim_totrecoe = pim_totrecoe + pim.E();
                }
              for(int iPip = 0; iPip < pips.size(); iPip++)
                {
                  SetLorentzVector(pip,pips[iPip]);
                  pip_totrecoe = pip_totrecoe + pip.E() - mass_pip;
                }
              double totrecoe = proton_totrecoe + pim_totrecoe + pip_totrecoe + el.E();
              //Fill 2D histograms
              totrecoe_W_1pXn1pim_MC->Fill(W,totrecoe);
              totrecoe_omega_1pXn1pim_MC->Fill(q.E(),totrecoe);
              totrecoe_q2_1pXn1pim_MC->Fill(q2,totrecoe);
              totrecoe_xb_1pXn1pim_MC->Fill(x_b,totrecoe);
              totrecoe_y_1pXn1pim_MC->Fill(y,totrecoe);
              W_q2_1pXn1pim_MC->Fill(q2,W);
              W_xb_1pXn1pim_MC->Fill(x_b,W);
              W_y_1pXn1pim_MC->Fill(y,W);
              q2_xb_1pXn1pim_MC->Fill(x_b,q2);
              q2_y_1pXn1pim_MC->Fill(y,q2);
              y_xb_1pXn1pim_MC->Fill(x_b,y);
            }*/
            //1pXn1pim--no neutron/neutral particle cuts
          if ( 
               protons.size()==1
            && (   //Must check that each proton has a good chi2 value range in either the forward or central detectors
                   ((protons[0]->getRegion()==FD) && (protons[0]->par()->getChi2Pid() < p_chi2_fd_cut_hi) && (protons[0]->par()->getChi2Pid() > p_chi2_fd_cut_lo))
                || ((protons[0]->getRegion()==CD) && (protons[0]->par()->getChi2Pid() < p_chi2_cd_cut_hi) && (protons[0]->par()->getChi2Pid() > p_chi2_cd_cut_lo))
               )
            && pips.size()==0
            && pims.size()==1 
            && (   //Must check that each pim has a good chi2 value range in either the forward or central detectors
                   ((pims[0]->getRegion()==FD) && (pims[0]->par()->getChi2Pid() < pim_chi2_fd_cut_hi) && (pims[0]->par()->getChi2Pid() > pim_chi2_fd_cut_lo))
                || ((pims[0]->getRegion()==CD) && (pims[0]->par()->getChi2Pid() < pim_chi2_cd_cut_hi) && (pims[0]->par()->getChi2Pid() > pim_chi2_cd_cut_lo))
               )
             )
            {              
              double proton_totrecoe = 0., neutron_totrecoe = 0., pim_totrecoe = 0., pip_totrecoe=0.;
              for(int iPr = 0; iPr < protons.size(); iPr++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  proton_totrecoe = proton_totrecoe + pr.E() - mass_p;
                }
              for(int iPim = 0; iPim < pims.size(); iPim++)
                {
                  SetLorentzVector(pim,pims[iPim]);
                  pim_totrecoe = pim_totrecoe + pim.E();
                }
              double totrecoe = proton_totrecoe + neutron_totrecoe + pim_totrecoe + pip_totrecoe + el.E();
              //q2_1pXn1pim_MC->Fill(q2);
              //totrecoe_1pXn1pim_MC->Fill(totrecoe);
              //W_1pXn1pim_MC->Fill(W);
              //Find invariant mass of the system (help to see if it comes from a true event?)
              TLorentzVector inv_mass_1pXn1pim = pr + pip + pim - q;
              nuc_inv_mass_W_1pXn1pim_MC->Fill(W,inv_mass_1pXn1pim.M());
              nuc_inv_mass_omega_1pXn1pim_MC->Fill(q.E(),inv_mass_1pXn1pim.M());
              nuc_inv_mass_q2_1pXn1pim_MC->Fill(q2,inv_mass_1pXn1pim.M());
              nuc_inv_mass_xb_1pXn1pim_MC->Fill(x_b,inv_mass_1pXn1pim.M());
              nuc_inv_mass_y_1pXn1pim_MC->Fill(y,inv_mass_1pXn1pim.M());
              //Fill 2D histograms
              totrecoe_W_1pXn1pim_MC->Fill(W,totrecoe);
              totrecoe_omega_1pXn1pim_MC->Fill(q.E(),totrecoe);
              totrecoe_q2_1pXn1pim_MC->Fill(q2,totrecoe);
              totrecoe_xb_1pXn1pim_MC->Fill(x_b,totrecoe);
              totrecoe_y_1pXn1pim_MC->Fill(y,totrecoe);
              W_q2_1pXn1pim_MC->Fill(q2,W);
              W_xb_1pXn1pim_MC->Fill(x_b,W);
              W_y_1pXn1pim_MC->Fill(y,W);
              q2_xb_1pXn1pim_MC->Fill(x_b,q2);
              q2_y_1pXn1pim_MC->Fill(y,q2);
              y_xb_1pXn1pim_MC->Fill(x_b,y);
            }
          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          //1D Histograms of particular topologies with cuts
          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          //QE single proton, no pions: only one proton and without charged pions, y:(-300,300), x:(0.8,1.2), W:(0.800,1.050), q2:(0.0,3.0)
          //Truth values
          if (protons.size()==1 && pips.size()==0 && pims.size()==0 && processid==1)
            {
              W_qe_1p_truth_MC->Fill(W);
              xb_qe_1p_truth_MC->Fill(x_b);
              counter_qe_1p0pi_truth_MC++;
              for(int iPr = 0; iPr < protons.size(); iPr++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  double p_proton = pr.P();
                  p_proton_qe_1p_truth_MC->Fill(p_proton);
                }
            }
          //Cuts (to be like data)
          if (protons.size()==1 && pips.size()==0 && pims.size()==0 && y>-0.300 && y<0.300 && x_b>0.8 && x_b<1.2 && W>0.8 && W<1.050 && q2>0. && q2<3.0)
            {
              W_qe_1p_cuts_MC->Fill(W);
              xb_qe_1p_cuts_MC->Fill(x_b);
              counter_qe_1p0pi_cuts_MC++;
              for(int iPr = 0; iPr < protons.size(); iPr++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  double p_proton = pr.P();
                  p_proton_qe_1p_cuts_MC->Fill(p_proton);
                }
            }

          //QE neutron: only one neutron, y:(-300,300), x:(0.8,1.2), W:(0.800,1.050), q2:(0.0,3.0)


          //MEC single proton, no charged pions
          //Truth values
          if (protons.size()==1 && pips.size()==0 && pims.size()==0 && processid==2)
            {
              W_mec_1p_truth_MC->Fill(W);
              xb_mec_1p_truth_MC->Fill(x_b);
              counter_mec_1p0pi_truth_MC++;
              for(int iPr = 0; iPr < protons.size(); iPr++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  double p_proton = pr.P();
                  p_proton_mec_1p_truth_MC->Fill(p_proton);
                }
            }
          //Cuts (to be like data)
          if (protons.size()==1 && pips.size()==0 && pims.size()==0 && y>0.05 && y<0.300 && x_b>0.8 && x_b<1.2 && W>0.8 && W<1.5 && q2>0. && q2<3.0)
            {
              W_mec_1p_cuts_MC->Fill(W);
              xb_mec_1p_cuts_MC->Fill(x_b);
              counter_mec_1p0pi_cuts_MC++;
              for(int iPr = 0; iPr < protons.size(); iPr++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  double p_proton = pr.P();
                  p_proton_mec_1p_cuts_MC->Fill(p_proton);
                }
            }

          //MEC two proton (proton-proton), no charged pions
          //Truth values
          if (protons.size()==2 && pips.size()==0 && pims.size()==0 && processid==2)
            {
              W_mec_2p_truth_MC->Fill(W);
              xb_mec_2p_truth_MC->Fill(x_b);
              counter_mec_2p0pi_truth_MC++;
              for(int iPr = 0; iPr < protons.size(); iPr++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  double p_proton = pr.P();
                  p_proton_mec_2p_truth_MC->Fill(p_proton);
                }
            }
          //Cuts (to be like data)
          if (protons.size()==2 && pips.size()==0 && pims.size()==0 && y>0.05 && y<0.300 && x_b>0.8 && x_b<1.2 && W>0.8 && W<1.5 && q2>0. && q2<3.0)
            {
              W_mec_2p_cuts_MC->Fill(W);
              xb_mec_2p_cuts_MC->Fill(x_b);
              counter_mec_2p0pi_cuts_MC++;
              for(int iPr = 0; iPr < protons.size(); iPr++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  double p_proton = pr.P();
                  p_proton_mec_2p_cuts_MC->Fill(p_proton);
                }
            }

          //MEC proton-neutron


          //RES single charged pion (no protons)
          //Truth values
          if (protons.size()==0 && ((pips.size()==1 && pims.size()==0) || (pips.size()==0 && pims.size()==1)) && processid==3 && resid>=0)
            {
              W_res_0p1pi_truth_MC->Fill(W);
              xb_res_0p1pi_truth_MC->Fill(x_b);
              counter_res_0p1pi_truth_MC++;
              for(int iPr = 0; iPr < pims.size()+pips.size(); iPr++)
                {
                  if (pips.size()==1 && pims.size()==0)
                    {
                      SetLorentzVector(pip,pips[iPr]);
                      double p_pip = pip.P();
                      p_pi_res_0p1pi_truth_MC->Fill(p_pip);
                    }
                  if (pips.size()==0 && pims.size()==1)
                    {
                      SetLorentzVector(pim,pims[iPr]);
                      double p_pim = pim.P();
                      p_pi_res_0p1pi_truth_MC->Fill(p_pim);
                    }
                }
            }
          //Cuts (to be like data)
          if (protons.size()==0 && ((pips.size()==1 && pims.size()==0) || (pips.size()==0 && pims.size()==1)) && y>0.2 && y<0.6 && x_b>0.4 && x_b<0.8 && W>1.050 && W<2.0 && q2>0.8 && q2<6.0)
            {
              W_res_0p1pi_cuts_MC->Fill(W);
              xb_res_0p1pi_cuts_MC->Fill(x_b);
              counter_res_0p1pi_cuts_MC++;
              for(int iPr = 0; iPr < pims.size()+pips.size(); iPr++)
                {
                  if (pips.size()==1 && pims.size()==0)
                    {
                      SetLorentzVector(pip,pips[iPr]);
                      double p_pip = pip.P();
                      p_pi_res_0p1pi_cuts_MC->Fill(p_pip);
                    }
                  if (pips.size()==0 && pims.size()==1)
                    {
                      SetLorentzVector(pim,pims[iPr]);
                      double p_pim = pim.P();
                      p_pi_res_0p1pi_cuts_MC->Fill(p_pim);
                    }
                }
            }

          //RES single charged pion with single proton
          //Truth values
          if (protons.size()==1 && ((pips.size()==1 && pims.size()==0) || (pips.size()==0 && pims.size()==1)) && processid==3 && resid>=0)
            {
              W_res_1p1pi_truth_MC->Fill(W);
              xb_res_1p1pi_truth_MC->Fill(x_b);
              counter_res_1p1pi_truth_MC++;
              for(int iPr = 0; iPr < protons.size(); iPr++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  double p_proton = pr.P();
                  p_proton_res_1p1pi_truth_MC->Fill(p_proton);
                }
              for(int iPr = 0; iPr < pims.size()+pips.size(); iPr++)
                {
                  if (pips.size()==1 && pims.size()==0)
                    {
                      SetLorentzVector(pip,pips[iPr]);
                      double p_pip = pip.P();
                      p_pi_res_1p1pi_truth_MC->Fill(p_pip);
                    }
                  if (pips.size()==0 && pims.size()==1)
                    {
                      SetLorentzVector(pim,pims[iPr]);
                      double p_pim = pim.P();
                      p_pi_res_1p1pi_truth_MC->Fill(p_pim);
                    }
                }
            }
          //Cuts (to be like data)
          if (protons.size()==1 && ((pips.size()==1 && pims.size()==0) || (pips.size()==0 && pims.size()==1)) && y>0.2 && y<0.6 && x_b>0.4 && x_b<0.8 && W>1.050 && W<2.0 && q2>0.8 && q2<6.0)
            {
              W_res_1p1pi_cuts_MC->Fill(W);
              xb_res_1p1pi_cuts_MC->Fill(x_b);
              counter_res_1p1pi_cuts_MC++;
              for(int iPr = 0; iPr < protons.size(); iPr++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  double p_proton = pr.P();
                  p_proton_res_1p1pi_cuts_MC->Fill(p_proton);
                }
              for(int iPr = 0; iPr < pims.size()+pips.size(); iPr++)
                {
                  if (pips.size()==1 && pims.size()==0)
                    {
                      SetLorentzVector(pip,pips[iPr]);
                      double p_pip = pip.P();
                      p_pi_res_1p1pi_cuts_MC->Fill(p_pip);
                    }
                  if (pips.size()==0 && pims.size()==1)
                    {
                      SetLorentzVector(pim,pims[iPr]);
                      double p_pim = pim.P();
                      p_pi_res_1p1pi_cuts_MC->Fill(p_pim);
                    }
                }
            }
          /*************  Plot for Larry: MC  *************/ 
          //0p1n0pi
          if ( 
               protons.size()==0
            && neutrons.size()==1
            /*&& (   //Must check that each neutron has a good chi2 value range in either the forward or central detectors
                   ((neutrons[0]->getRegion()==FD) && (neutrons[0]->par()->getChi2Pid() < n_chi2_fd_cut_hi) && (neutrons[0]->par()->getChi2Pid() > n_chi2_fd_cut_lo))
                || ((neutrons[0]->getRegion()==CD) && (neutrons[0]->par()->getChi2Pid() < n_chi2_cd_cut_hi) && (neutrons[0]->par()->getChi2Pid() > n_chi2_cd_cut_lo))
               )*/
            && pips.size()==0
            && pims.size()==0 
             )
            {
              double proton_totrecoe = 0., neutron_totrecoe = 0., pim_totrecoe = 0., pip_totrecoe=0.;
              for(int iNr = 0; iNr < neutrons.size(); iNr++)
                {
                   SetLorentzVector(nr,neutrons[iNr]);
                  neutron_totrecoe = neutron_totrecoe + nr.E() - mass_n;
                }
              double totrecoe = proton_totrecoe + neutron_totrecoe + pim_totrecoe + pip_totrecoe + el.E();
              q2_0p1n0pi_MC->Fill(q2);
              totrecoe_0p1n0pi_MC->Fill(totrecoe-0.08);
              W_0p1n0pi_MC->Fill(W);
            }
          //1pXn0pi--no neutron/neutral particle cuts
          if ( 
               protons.size()==1
            && (   //Must check that each proton has a good chi2 value range in either the forward or central detectors
                   ((protons[0]->getRegion()==FD) && (protons[0]->par()->getChi2Pid() < p_chi2_fd_cut_hi) && (protons[0]->par()->getChi2Pid() > p_chi2_fd_cut_lo))
                || ((protons[0]->getRegion()==CD) && (protons[0]->par()->getChi2Pid() < p_chi2_cd_cut_hi) && (protons[0]->par()->getChi2Pid() > p_chi2_cd_cut_lo))
               )
            && pips.size()==0
            && pims.size()==0 
             )
            {
              double proton_totrecoe = 0., neutron_totrecoe = 0., pim_totrecoe = 0., pip_totrecoe=0.;
              for(int iPr = 0; iPr < protons.size(); iPr++)
                {
                  SetLorentzVector(pr,protons[iPr]);
                  proton_totrecoe = proton_totrecoe + pr.E() - mass_p;
                }
              double totrecoe = proton_totrecoe + neutron_totrecoe + pim_totrecoe + pip_totrecoe + el.E();
              q2_1pXn0pi_MC->Fill(q2);
              totrecoe_1pXn0pi_MC->Fill(totrecoe-0.08);
              W_1pXn0pi_MC->Fill(W);
            }


          //RES pion-neutron


          //RES pion-pion-neutron

          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	        counter_MC++;
	      }
    }

  //Fudge some scaling
  double MC_to_Data_Scaling = ( (double)counter_data ) / ( (double)counter_MC );
  //double MC_to_Data_qe_Scaling = ( (double)counter_qe_MC) * ( () );
  //Scale, color, and place all histograms
  /*************  Plot for Larry: MC  *************/
  q2_0p1n0pi_MC->SetLineColor(2); q2_0p1n0pi_MC->SetMarkerColor(2); q2_0p1n0pi_MC->Sumw2(); q2_0p1n0pi_MC->Scale(MC_to_Data_Scaling); /*q2_byChannel_stacked->Add(q2_0p1n0pi_MC);*/
  q2_1pXn0pi_MC->SetLineColor(3); q2_1pXn0pi_MC->SetMarkerColor(3); q2_1pXn0pi_MC->Sumw2(); q2_1pXn0pi_MC->Scale(MC_to_Data_Scaling); /*q2_byChannel_stacked->Add(q2_1pXn0pi_MC);*/
  totrecoe_0p1n0pi_MC->SetLineColor(2); totrecoe_0p1n0pi_MC->SetMarkerColor(2); totrecoe_0p1n0pi_MC->Sumw2(); totrecoe_0p1n0pi_MC->Scale(MC_to_Data_Scaling); /*totrecoe_byChannel_stacked->Add(totrecoe_0p1n0pi_MC);*/
  totrecoe_1pXn0pi_MC->SetLineColor(3); totrecoe_1pXn0pi_MC->SetMarkerColor(3); totrecoe_1pXn0pi_MC->Sumw2(); totrecoe_1pXn0pi_MC->Scale(MC_to_Data_Scaling); /*totrecoe_byChannel_stacked->Add(totrecoe_1pXn0pi_MC);*/
  W_0p1n0pi_MC->SetLineColor(2); W_0p1n0pi_MC->SetMarkerColor(2); W_0p1n0pi_MC->Sumw2(); W_0p1n0pi_MC->Scale(MC_to_Data_Scaling); /*W_byChannel_stacked->Add(W_0p1n0pi_MC);*/
  W_1pXn0pi_MC->SetLineColor(3); W_1pXn0pi_MC->SetMarkerColor(3); W_1pXn0pi_MC->Sumw2(); W_1pXn0pi_MC->Scale(MC_to_Data_Scaling); /*W_byChannel_stacked->Add(W_1pXn0pi_MC);*/
  //Main histograms
  mult_p_MC->SetMarkerColor(kRed); mult_p_MC->Sumw2(); mult_p_MC->Scale(MC_to_Data_Scaling); mult_p_MC->SetLineColor(kRed); mult_p_stacked->Add(mult_p_MC); //mult_p_MC->Clear();
  //mult_n_MC->SetMarkerColor(kRed); mult_n_MC->Sumw2(); mult_n_MC->Scale(MC_to_Data_Scaling); mult_n_MC->SetLineColor(kRed); mult_n_MC->SetLineStyle(1); mult_n_stacked->Add(mult_n_MC);
  mult_n_nocuts_MC->SetMarkerColor(kRed); mult_n_nocuts_MC->Sumw2(); mult_n_nocuts_MC->Scale(MC_to_Data_Scaling); mult_n_nocuts_MC->SetLineColor(kRed); mult_n_nocuts_MC->SetLineStyle(9); mult_n_stacked->Add(mult_n_nocuts_MC);
  mult_pip_MC->SetMarkerColor(kRed); mult_pip_MC->Sumw2(); mult_pip_MC->Scale(MC_to_Data_Scaling); mult_pip_MC->SetLineColor(kRed); mult_pip_stacked->Add(mult_pip_MC); //mult_pip_MC->Clear();
  mult_pim_MC->SetMarkerColor(kRed); mult_pim_MC->Sumw2(); mult_pim_MC->Scale(MC_to_Data_Scaling); mult_pim_MC->SetLineColor(kRed); mult_pim_stacked->Add(mult_pim_MC); //mult_pim_MC->Clear();
  //q2_mc_h_MC->SetMarkerColor(kRed); q2_mc_h_MC->Sumw2(); q2_mc_h_MC->Scale(MC_to_Data_Scaling); q2_mc_h_MC->SetLineColor(kRed); q2_mc_h_stacked->Add(q2_mc_h_MC); //q2_mc_h_MC->Clear();
  q2_MC->SetMarkerColor(kRed); q2_MC->Sumw2(); q2_MC->Scale(MC_to_Data_Scaling); q2_MC->SetLineColor(kRed); q2_stacked->Add(q2_MC); //q2_MC->Clear();
  q2_qe_MC->SetMarkerColor(kGreen); q2_qe_MC->Sumw2(); q2_qe_MC->Scale(MC_to_Data_Scaling); q2_qe_MC->SetLineColor(kGreen); q2_stacked->Add(q2_qe_MC);
  q2_mec_MC->SetMarkerColor(kOrange); q2_mec_MC->Sumw2(); q2_mec_MC->Scale(MC_to_Data_Scaling); q2_mec_MC->SetLineColor(kOrange); q2_stacked->Add(q2_mec_MC);
  q2_res_MC->SetMarkerColor(kBlue); q2_res_MC->Sumw2(); q2_res_MC->Scale(MC_to_Data_Scaling); q2_res_MC->SetLineColor(kBlue); q2_stacked->Add(q2_res_MC);
  q2_dis_MC->SetMarkerColor(kCyan); q2_dis_MC->Sumw2(); q2_dis_MC->Scale(MC_to_Data_Scaling); q2_dis_MC->SetLineColor(kCyan); q2_stacked->Add(q2_dis_MC);
  p_chi2_cd_MC->SetMarkerColor(kRed); p_chi2_cd_MC->Sumw2(); p_chi2_cd_MC->Scale(MC_to_Data_Scaling); p_chi2_cd_MC->SetLineColor(kRed); p_chi2_cd_MC->SetLineStyle(1); p_chi2_stacked->Add(p_chi2_cd_MC); //p_chi2_cd_MC->Clear();
  n_chi2_cd_MC->SetMarkerColor(kRed); n_chi2_cd_MC->Sumw2(); n_chi2_cd_MC->Scale(MC_to_Data_Scaling); n_chi2_cd_MC->SetLineColor(kRed); n_chi2_cd_MC->SetLineStyle(1); n_chi2_stacked->Add(n_chi2_cd_MC);
  n_chi2_nocuts_MC->SetMarkerColor(kGreen); n_chi2_nocuts_MC->Sumw2(); n_chi2_nocuts_MC->Scale(MC_to_Data_Scaling); n_chi2_nocuts_MC->SetLineColor(kGreen); n_chi2_nocuts_MC->SetLineStyle(1); n_chi2_stacked->Add(n_chi2_nocuts_MC);
  pip_chi2_cd_MC->SetMarkerColor(kRed); pip_chi2_cd_MC->Sumw2(); pip_chi2_cd_MC->Scale(MC_to_Data_Scaling); pip_chi2_cd_MC->SetLineColor(kRed); pip_chi2_cd_MC->SetLineStyle(1); pip_chi2_stacked->Add(pip_chi2_cd_MC); //pip_chi2_cd_MC->Clear();
  pim_chi2_cd_MC->SetMarkerColor(kRed); pim_chi2_cd_MC->Sumw2(); pim_chi2_cd_MC->Scale(MC_to_Data_Scaling); pim_chi2_cd_MC->SetLineColor(kRed); pim_chi2_cd_MC->SetLineStyle(1); pim_chi2_stacked->Add(pim_chi2_cd_MC); //pim_chi2_cd_MC->Clear();
  p_chi2_fd_MC->SetMarkerColor(kRed); p_chi2_fd_MC->Sumw2(); p_chi2_fd_MC->Scale(MC_to_Data_Scaling); p_chi2_fd_MC->SetLineColor(kRed); p_chi2_fd_MC->SetLineStyle(9); p_chi2_stacked->Add(p_chi2_fd_MC); //p_chi2_fd_MC->Clear();
  n_chi2_fd_MC->SetMarkerColor(kRed); n_chi2_fd_MC->Sumw2(); n_chi2_fd_MC->Scale(MC_to_Data_Scaling); n_chi2_fd_MC->SetLineColor(kRed); n_chi2_fd_MC->SetLineStyle(9); n_chi2_stacked->Add(n_chi2_fd_MC);
  pip_chi2_fd_MC->SetMarkerColor(kRed); pip_chi2_fd_MC->Sumw2(); pip_chi2_fd_MC->Scale(MC_to_Data_Scaling); pip_chi2_fd_MC->SetLineColor(kRed); pip_chi2_fd_MC->SetLineStyle(9); pip_chi2_stacked->Add(pip_chi2_fd_MC); //pip_chi2_fd_MC->Clear();
  pim_chi2_fd_MC->SetMarkerColor(kRed); pim_chi2_fd_MC->Sumw2(); pim_chi2_fd_MC->Scale(MC_to_Data_Scaling); pim_chi2_fd_MC->SetLineColor(kRed); pim_chi2_fd_MC->SetLineStyle(9); pim_chi2_stacked->Add(pim_chi2_fd_MC); //pim_chi2_fd_MC->Clear();
  energy_transfer_0_MC->SetMarkerColor(kRed); energy_transfer_0_MC->Sumw2(); energy_transfer_0_MC->Scale(MC_to_Data_Scaling); energy_transfer_0_MC->SetLineColor(kRed); energy_transfer_0_stacked->Add(energy_transfer_0_MC); //energy_transfer_0_MC->Clear();
  energy_transfer_1_MC->SetMarkerColor(kRed); energy_transfer_1_MC->Sumw2(); energy_transfer_1_MC->Scale(MC_to_Data_Scaling); energy_transfer_1_MC->SetLineColor(kRed); energy_transfer_1_stacked->Add(energy_transfer_1_MC); //energy_transfer_1_MC->Clear();
  energy_transfer_2_MC->SetMarkerColor(kRed); energy_transfer_2_MC->Sumw2(); energy_transfer_2_MC->Scale(MC_to_Data_Scaling); energy_transfer_2_MC->SetLineColor(kRed); energy_transfer_2_stacked->Add(energy_transfer_2_MC); //energy_transfer_2_MC->Clear();
  energy_transfer_3_MC->SetMarkerColor(kRed); energy_transfer_3_MC->Sumw2(); energy_transfer_3_MC->Scale(MC_to_Data_Scaling); energy_transfer_3_MC->SetLineColor(kRed); energy_transfer_3_stacked->Add(energy_transfer_3_MC); //energy_transfer_3_MC->Clear();
  energy_transfer_qe_0_MC->SetMarkerColor(kGreen); energy_transfer_qe_0_MC->Sumw2(); energy_transfer_qe_0_MC->Scale(MC_to_Data_Scaling); energy_transfer_qe_0_MC->SetLineColor(kGreen); energy_transfer_0_stacked->Add(energy_transfer_qe_0_MC);
  energy_transfer_qe_1_MC->SetMarkerColor(kGreen); energy_transfer_qe_1_MC->Sumw2(); energy_transfer_qe_1_MC->Scale(MC_to_Data_Scaling); energy_transfer_qe_1_MC->SetLineColor(kGreen); energy_transfer_1_stacked->Add(energy_transfer_qe_1_MC);
  energy_transfer_qe_2_MC->SetMarkerColor(kGreen); energy_transfer_qe_2_MC->Sumw2(); energy_transfer_qe_2_MC->Scale(MC_to_Data_Scaling); energy_transfer_qe_2_MC->SetLineColor(kGreen); energy_transfer_2_stacked->Add(energy_transfer_qe_2_MC);
  energy_transfer_qe_3_MC->SetMarkerColor(kGreen); energy_transfer_qe_3_MC->Sumw2(); energy_transfer_qe_3_MC->Scale(MC_to_Data_Scaling); energy_transfer_qe_3_MC->SetLineColor(kGreen); energy_transfer_3_stacked->Add(energy_transfer_qe_3_MC);
  energy_transfer_mec_0_MC->SetMarkerColor(kOrange); energy_transfer_mec_0_MC->Sumw2(); energy_transfer_mec_0_MC->Scale(MC_to_Data_Scaling); energy_transfer_mec_0_MC->SetLineColor(kOrange); energy_transfer_0_stacked->Add(energy_transfer_mec_0_MC);
  energy_transfer_mec_1_MC->SetMarkerColor(kOrange); energy_transfer_mec_1_MC->Sumw2(); energy_transfer_mec_1_MC->Scale(MC_to_Data_Scaling); energy_transfer_mec_1_MC->SetLineColor(kOrange); energy_transfer_1_stacked->Add(energy_transfer_mec_1_MC);
  energy_transfer_mec_2_MC->SetMarkerColor(kOrange); energy_transfer_mec_2_MC->Sumw2(); energy_transfer_mec_2_MC->Scale(MC_to_Data_Scaling); energy_transfer_mec_2_MC->SetLineColor(kOrange); energy_transfer_2_stacked->Add(energy_transfer_mec_2_MC);
  energy_transfer_mec_3_MC->SetMarkerColor(kOrange); energy_transfer_mec_3_MC->Sumw2(); energy_transfer_mec_3_MC->Scale(MC_to_Data_Scaling); energy_transfer_mec_3_MC->SetLineColor(kOrange); energy_transfer_3_stacked->Add(energy_transfer_mec_3_MC);
  energy_transfer_res_0_MC->SetMarkerColor(kBlue); energy_transfer_res_0_MC->Sumw2(); energy_transfer_res_0_MC->Scale(MC_to_Data_Scaling); energy_transfer_res_0_MC->SetLineColor(kBlue); energy_transfer_0_stacked->Add(energy_transfer_res_0_MC);
  energy_transfer_res_1_MC->SetMarkerColor(kBlue); energy_transfer_res_1_MC->Sumw2(); energy_transfer_res_1_MC->Scale(MC_to_Data_Scaling); energy_transfer_res_1_MC->SetLineColor(kBlue); energy_transfer_1_stacked->Add(energy_transfer_res_1_MC);
  energy_transfer_res_2_MC->SetMarkerColor(kBlue); energy_transfer_res_2_MC->Sumw2(); energy_transfer_res_2_MC->Scale(MC_to_Data_Scaling); energy_transfer_res_2_MC->SetLineColor(kBlue); energy_transfer_2_stacked->Add(energy_transfer_res_2_MC);
  energy_transfer_res_3_MC->SetMarkerColor(kBlue); energy_transfer_res_3_MC->Sumw2(); energy_transfer_res_3_MC->Scale(MC_to_Data_Scaling); energy_transfer_res_3_MC->SetLineColor(kBlue); energy_transfer_3_stacked->Add(energy_transfer_res_3_MC);
  energy_transfer_dis_0_MC->SetMarkerColor(kCyan); energy_transfer_dis_0_MC->Sumw2(); energy_transfer_dis_0_MC->Scale(MC_to_Data_Scaling); energy_transfer_dis_0_MC->SetLineColor(kCyan); energy_transfer_0_stacked->Add(energy_transfer_dis_0_MC);
  energy_transfer_dis_1_MC->SetMarkerColor(kCyan); energy_transfer_dis_1_MC->Sumw2(); energy_transfer_dis_1_MC->Scale(MC_to_Data_Scaling); energy_transfer_dis_1_MC->SetLineColor(kCyan); energy_transfer_1_stacked->Add(energy_transfer_dis_1_MC);
  energy_transfer_dis_2_MC->SetMarkerColor(kCyan); energy_transfer_dis_2_MC->Sumw2(); energy_transfer_dis_2_MC->Scale(MC_to_Data_Scaling); energy_transfer_dis_2_MC->SetLineColor(kCyan); energy_transfer_2_stacked->Add(energy_transfer_dis_2_MC);
  energy_transfer_dis_3_MC->SetMarkerColor(kCyan); energy_transfer_dis_3_MC->Sumw2(); energy_transfer_dis_3_MC->Scale(MC_to_Data_Scaling); energy_transfer_dis_3_MC->SetLineColor(kCyan); energy_transfer_3_stacked->Add(energy_transfer_dis_3_MC);
  invariant_mass_MC->SetMarkerColor(kRed); invariant_mass_MC->Sumw2(); invariant_mass_MC->Scale(MC_to_Data_Scaling); invariant_mass_MC->SetLineColor(kRed); invariant_mass_stacked->Add(invariant_mass_MC); //invariant_mass_MC->Clear();
  el_vz_h_MC->SetMarkerColor(kRed); el_vz_h_MC->Sumw2(); el_vz_h_MC->Scale(MC_to_Data_Scaling); el_vz_h_MC->SetLineColor(kRed); el_vz_h_stacked->Add(el_vz_h_MC); //el_vz_h_MC->Clear();
  W_MC->SetMarkerColor(kRed); W_MC->Sumw2(); W_MC->Scale(MC_to_Data_Scaling); W_MC->SetLineColor(kRed); W_stacked->Add(W_MC); //W_MC->Clear();
  W_qe_MC->SetMarkerColor(kGreen); W_qe_MC->Sumw2(); W_qe_MC->Scale(MC_to_Data_Scaling); W_qe_MC->SetLineColor(kGreen); W_stacked->Add(W_qe_MC);
  W_mec_MC->SetMarkerColor(kOrange); W_mec_MC->Sumw2(); W_mec_MC->Scale(MC_to_Data_Scaling); W_mec_MC->SetLineColor(kOrange); W_stacked->Add(W_mec_MC);
  W_res_MC->SetMarkerColor(kBlue); W_res_MC->Sumw2(); W_res_MC->Scale(MC_to_Data_Scaling); W_res_MC->SetLineColor(kBlue); W_stacked->Add(W_res_MC);
  W_resonances_1p2_to_2p1_MC->SetMarkerColor(kBlue); W_resonances_1p2_to_2p1_MC->Sumw2(); W_resonances_1p2_to_2p1_MC->Scale(MC_to_Data_Scaling); W_resonances_1p2_to_2p1_MC->SetMarkerColor(kBlue);
  W_dis_MC->SetMarkerColor(kCyan); W_dis_MC->Sumw2(); W_dis_MC->Scale(MC_to_Data_Scaling); W_dis_MC->SetLineColor(kCyan); W_stacked->Add(W_dis_MC);
  xb_MC->SetMarkerColor(kRed); xb_MC->Sumw2(); xb_MC->Scale(MC_to_Data_Scaling); xb_MC->SetLineColor(kRed); xb_stacked->Add(xb_MC);
  xb_qe_MC->SetMarkerColor(kGreen); xb_qe_MC->Sumw2(); xb_qe_MC->Scale(MC_to_Data_Scaling); xb_qe_MC->SetLineColor(kGreen); xb_stacked->Add(xb_qe_MC);
  xb_mec_MC->SetMarkerColor(kOrange); xb_mec_MC->Sumw2(); xb_mec_MC->Scale(MC_to_Data_Scaling); xb_mec_MC->SetLineColor(kOrange); xb_stacked->Add(xb_mec_MC);
  xb_res_MC->SetMarkerColor(kBlue); xb_res_MC->Sumw2(); xb_res_MC->Scale(MC_to_Data_Scaling); xb_res_MC->SetLineColor(kBlue); xb_stacked->Add(xb_res_MC);
  xb_dis_MC->SetMarkerColor(kCyan); xb_dis_MC->Sumw2(); xb_dis_MC->Scale(MC_to_Data_Scaling); xb_dis_MC->SetLineColor(kCyan); xb_stacked->Add(xb_dis_MC);
  y_MC->SetMarkerColor(kRed); y_MC->Sumw2(); y_MC->Scale(MC_to_Data_Scaling); y_MC->SetLineColor(kRed); y_stacked->Add(y_MC);
  y_qe_MC->SetMarkerColor(kGreen); y_qe_MC->Sumw2(); y_qe_MC->Scale(MC_to_Data_Scaling); y_qe_MC->SetLineColor(kGreen); y_stacked->Add(y_qe_MC);
  y_mec_MC->SetMarkerColor(kOrange); y_mec_MC->Sumw2(); y_mec_MC->Scale(MC_to_Data_Scaling); y_mec_MC->SetLineColor(kOrange); y_stacked->Add(y_mec_MC);
  y_res_MC->SetMarkerColor(kBlue); y_res_MC->Sumw2(); y_res_MC->Scale(MC_to_Data_Scaling); y_res_MC->SetLineColor(kBlue); y_stacked->Add(y_res_MC);
  y_dis_MC->SetMarkerColor(kCyan); y_dis_MC->Sumw2(); y_dis_MC->Scale(MC_to_Data_Scaling); y_dis_MC->SetLineColor(kCyan); y_stacked->Add(y_dis_MC);
  omega_MC->SetMarkerColor(kRed); omega_MC->Sumw2(); omega_MC->Scale(MC_to_Data_Scaling); omega_MC->SetLineColor(kRed); omega_stacked->Add(omega_MC);
  omega_qe_MC->SetMarkerColor(kGreen); omega_qe_MC->Sumw2(); omega_qe_MC->Scale(MC_to_Data_Scaling); omega_qe_MC->SetLineColor(kGreen); omega_stacked->Add(omega_qe_MC);
  omega_mec_MC->SetMarkerColor(kOrange); omega_mec_MC->Sumw2(); omega_mec_MC->Scale(MC_to_Data_Scaling); omega_mec_MC->SetLineColor(kOrange); omega_stacked->Add(omega_mec_MC);
  omega_res_MC->SetMarkerColor(kBlue); omega_res_MC->Sumw2(); omega_res_MC->Scale(MC_to_Data_Scaling); omega_res_MC->SetLineColor(kBlue); omega_stacked->Add(omega_res_MC);
  omega_dis_MC->SetMarkerColor(kCyan); omega_dis_MC->Sumw2(); omega_dis_MC->Scale(MC_to_Data_Scaling); omega_dis_MC->SetLineColor(kCyan); omega_stacked->Add(omega_dis_MC);
  theta_MC->SetMarkerColor(kRed); theta_MC->Sumw2(); theta_MC->Scale(MC_to_Data_Scaling); theta_MC->SetLineColor(kRed); theta_stacked->Add(theta_MC);
  theta_qe_MC->SetMarkerColor(kGreen); theta_qe_MC->Sumw2(); theta_qe_MC->Scale(MC_to_Data_Scaling); theta_qe_MC->SetLineColor(kGreen); theta_stacked->Add(theta_qe_MC);
  theta_mec_MC->SetMarkerColor(kOrange); theta_mec_MC->Sumw2(); theta_mec_MC->Scale(MC_to_Data_Scaling); theta_mec_MC->SetLineColor(kOrange); theta_stacked->Add(theta_mec_MC);
  theta_res_MC->SetMarkerColor(kBlue); theta_res_MC->Sumw2(); theta_res_MC->Scale(MC_to_Data_Scaling); theta_res_MC->SetLineColor(kBlue); theta_stacked->Add(theta_res_MC);
  theta_dis_MC->SetMarkerColor(kCyan); theta_dis_MC->Sumw2(); theta_dis_MC->Scale(MC_to_Data_Scaling); theta_dis_MC->SetLineColor(kCyan); theta_stacked->Add(theta_dis_MC);
  p_proton_MC->SetMarkerColor(kRed); p_proton_MC->Sumw2(); p_proton_MC->Scale(MC_to_Data_Scaling); p_proton_MC->SetLineColor(kRed); p_proton_stacked->Add(p_proton_MC);
  p_miss_MC->SetMarkerColor(kRed); p_miss_MC->Sumw2(); p_miss_MC->Scale(MC_to_Data_Scaling); p_miss_MC->SetLineColor(kRed); p_miss_stacked->Add(p_miss_MC);
  p_miss_theta_MC->SetMarkerColor(kRed); p_miss_theta_MC->Sumw2(); p_miss_theta_MC->Scale(MC_to_Data_Scaling); p_miss_theta_MC->SetLineColor(kRed); p_miss_theta_stacked->Add(p_miss_theta_MC);
  processid_MC->SetMarkerColor(kRed); processid_MC->Sumw2(); processid_MC->Scale(MC_to_Data_Scaling); processid_MC->SetLineColor(kRed);
  resid_MC->SetMarkerColor(kRed); resid_MC->Sumw2(); /*resid_MC->Scale(MC_to_Data_Scaling);*/ resid_MC->SetLineColor(kRed);
  resonances_1p2_to_2p1_MC->SetMarkerColor(kBlue); resonances_1p2_to_2p1_MC->Sumw2(); resonances_1p2_to_2p1_MC->Scale(MC_to_Data_Scaling); resonances_1p2_to_2p1_MC->SetLineColor(kBlue); 

  //1D Plots by sector
  //W
  W_bySector_MC->SetMarkerColor(1); W_bySector_MC->Sumw2(); W_bySector_MC->Scale(MC_to_Data_Scaling); W_bySector_MC->SetLineColor(1); W_bySector_MC_stacked->Add(W_bySector_MC);
  W_S1_MC->SetLineColor(2); W_S1_MC->Sumw2(); W_S1_MC->Scale(6.*MC_to_Data_Scaling); W_S1_MC->SetLineColor(2); W_bySector_MC_stacked->Add(W_S1_MC);
  W_S2_MC->SetLineColor(3); W_S2_MC->Sumw2(); W_S2_MC->Scale(6.*MC_to_Data_Scaling); W_S2_MC->SetLineColor(3); W_bySector_MC_stacked->Add(W_S2_MC);
  W_S3_MC->SetLineColor(4); W_S3_MC->Sumw2(); W_S3_MC->Scale(6.*MC_to_Data_Scaling); W_S3_MC->SetLineColor(4); W_bySector_MC_stacked->Add(W_S3_MC);
  W_S4_MC->SetLineColor(5); W_S4_MC->Sumw2(); W_S4_MC->Scale(6.*MC_to_Data_Scaling); W_S4_MC->SetLineColor(5); W_bySector_MC_stacked->Add(W_S4_MC);
  W_S5_MC->SetLineColor(6); W_S5_MC->Sumw2(); W_S5_MC->Scale(6.*MC_to_Data_Scaling); W_S5_MC->SetLineColor(6); W_bySector_MC_stacked->Add(W_S5_MC);
  W_S6_MC->SetLineColor(7); W_S6_MC->Sumw2(); W_S6_MC->Scale(6.*MC_to_Data_Scaling); W_S6_MC->SetLineColor(7); W_bySector_MC_stacked->Add(W_S6_MC);
  //Q^2
  q2_bySector_MC->SetMarkerColor(1); q2_bySector_MC->Sumw2(); q2_bySector_MC->Scale(MC_to_Data_Scaling); q2_bySector_MC->SetLineColor(1); q2_bySector_MC_stacked->Add(q2_bySector_MC);
  q2_S1_MC->SetLineColor(2); q2_S1_MC->Sumw2(); q2_S1_MC->Scale(6.*MC_to_Data_Scaling); q2_S1_MC->SetLineColor(2); q2_bySector_MC_stacked->Add(q2_S1_MC);
  q2_S2_MC->SetLineColor(3); q2_S2_MC->Sumw2(); q2_S2_MC->Scale(6.*MC_to_Data_Scaling); q2_S2_MC->SetLineColor(3); q2_bySector_MC_stacked->Add(q2_S2_MC);
  q2_S3_MC->SetLineColor(4); q2_S3_MC->Sumw2(); q2_S3_MC->Scale(6.*MC_to_Data_Scaling); q2_S3_MC->SetLineColor(4); q2_bySector_MC_stacked->Add(q2_S3_MC);
  q2_S4_MC->SetLineColor(5); q2_S4_MC->Sumw2(); q2_S4_MC->Scale(6.*MC_to_Data_Scaling); q2_S4_MC->SetLineColor(5); q2_bySector_MC_stacked->Add(q2_S4_MC);
  q2_S5_MC->SetLineColor(6); q2_S5_MC->Sumw2(); q2_S5_MC->Scale(6.*MC_to_Data_Scaling); q2_S5_MC->SetLineColor(6); q2_bySector_MC_stacked->Add(q2_S5_MC);
  q2_S6_MC->SetLineColor(7); q2_S6_MC->Sumw2(); q2_S6_MC->Scale(6.*MC_to_Data_Scaling); q2_S6_MC->SetLineColor(7); q2_bySector_MC_stacked->Add(q2_S6_MC);
  
  //Histograms with kinematics cuts for particular topologies
  //QE Single Proton
  W_qe_1p_cuts_MC->SetLineColor(kBlue); W_qe_1p_cuts_MC->Sumw2(); W_qe_1p_cuts_MC->Scale(MC_to_Data_Scaling); W_qe_1p_cuts_MC->SetLineColor(kBlue); W_qe_1p_cuts_stacked->Add(W_qe_1p_cuts_MC);
  xb_qe_1p_cuts_MC->SetLineColor(kBlue); xb_qe_1p_cuts_MC->Sumw2(); xb_qe_1p_cuts_MC->Scale(MC_to_Data_Scaling); xb_qe_1p_cuts_MC->SetLineColor(kBlue); xb_qe_1p_cuts_stacked->Add(xb_qe_1p_cuts_MC);
  p_proton_qe_1p_cuts_MC->SetLineColor(kBlue); p_proton_qe_1p_cuts_MC->Sumw2(); p_proton_qe_1p_cuts_MC->Scale(MC_to_Data_Scaling); p_proton_qe_1p_cuts_MC->SetLineColor(kBlue); p_proton_qe_1p_cuts_stacked->Add(p_proton_qe_1p_cuts_MC);
  W_qe_1p_truth_MC->SetLineColor(kGreen); W_qe_1p_truth_MC->Sumw2(); W_qe_1p_truth_MC->Scale(MC_to_Data_Scaling); W_qe_1p_truth_MC->SetLineColor(kGreen); W_qe_1p_cuts_stacked->Add(W_qe_1p_truth_MC);
  xb_qe_1p_truth_MC->SetLineColor(kGreen); xb_qe_1p_truth_MC->Sumw2(); xb_qe_1p_truth_MC->Scale(MC_to_Data_Scaling); xb_qe_1p_truth_MC->SetLineColor(kGreen); xb_qe_1p_cuts_stacked->Add(xb_qe_1p_truth_MC);
  p_proton_qe_1p_truth_MC->SetLineColor(kGreen); p_proton_qe_1p_truth_MC->Sumw2(); p_proton_qe_1p_truth_MC->Scale(MC_to_Data_Scaling); p_proton_qe_1p_truth_MC->SetLineColor(kGreen); p_proton_qe_1p_cuts_stacked->Add(p_proton_qe_1p_truth_MC);
  
  //QE single neutron


  //MEC single proton
  W_mec_1p_cuts_MC->SetLineColor(kBlue); W_mec_1p_cuts_MC->Sumw2(); W_mec_1p_cuts_MC->Scale(MC_to_Data_Scaling); W_mec_1p_cuts_MC->SetLineColor(kBlue); W_mec_1p_cuts_stacked->Add(W_mec_1p_cuts_MC);
  xb_mec_1p_cuts_MC->SetLineColor(kBlue); xb_mec_1p_cuts_MC->Sumw2(); xb_mec_1p_cuts_MC->Scale(MC_to_Data_Scaling); xb_mec_1p_cuts_MC->SetLineColor(kBlue); xb_mec_1p_cuts_stacked->Add(xb_mec_1p_cuts_MC);
  p_proton_mec_1p_cuts_MC->SetLineColor(kBlue); p_proton_mec_1p_cuts_MC->Sumw2(); p_proton_mec_1p_cuts_MC->Scale(MC_to_Data_Scaling); p_proton_mec_1p_cuts_MC->SetLineColor(kBlue); p_proton_mec_1p_cuts_stacked->Add(p_proton_mec_1p_cuts_MC);
  W_mec_1p_truth_MC->SetLineColor(kGreen); W_mec_1p_truth_MC->Sumw2(); W_mec_1p_truth_MC->Scale(MC_to_Data_Scaling); W_mec_1p_truth_MC->SetLineColor(kGreen); W_mec_1p_cuts_stacked->Add(W_mec_1p_truth_MC);
  xb_mec_1p_truth_MC->SetLineColor(kGreen); xb_mec_1p_truth_MC->Sumw2(); xb_mec_1p_truth_MC->Scale(MC_to_Data_Scaling); xb_mec_1p_truth_MC->SetLineColor(kGreen); xb_mec_1p_cuts_stacked->Add(xb_mec_1p_truth_MC);
  p_proton_mec_1p_truth_MC->SetLineColor(kGreen); p_proton_mec_1p_truth_MC->Sumw2(); p_proton_mec_1p_truth_MC->Scale(MC_to_Data_Scaling); p_proton_mec_1p_truth_MC->SetLineColor(kGreen); p_proton_mec_1p_cuts_stacked->Add(p_proton_mec_1p_truth_MC);

  //MEC two proton (proton-proton)
  W_mec_2p_cuts_MC->SetLineColor(kBlue); W_mec_2p_cuts_MC->Sumw2(); W_mec_2p_cuts_MC->Scale(MC_to_Data_Scaling); W_mec_2p_cuts_MC->SetLineColor(kBlue); W_mec_2p_cuts_stacked->Add(W_mec_2p_cuts_MC);
  xb_mec_2p_cuts_MC->SetLineColor(kBlue); xb_mec_2p_cuts_MC->Sumw2(); xb_mec_2p_cuts_MC->Scale(MC_to_Data_Scaling); xb_mec_2p_cuts_MC->SetLineColor(kBlue); xb_mec_2p_cuts_stacked->Add(xb_mec_2p_cuts_MC);
  p_proton_mec_2p_cuts_MC->SetLineColor(kBlue); p_proton_mec_2p_cuts_MC->Sumw2(); p_proton_mec_2p_cuts_MC->Scale(MC_to_Data_Scaling); p_proton_mec_2p_cuts_MC->SetLineColor(kBlue); p_proton_mec_2p_cuts_stacked->Add(p_proton_mec_2p_cuts_MC);
  W_mec_2p_truth_MC->SetLineColor(kGreen); W_mec_2p_truth_MC->Sumw2(); W_mec_2p_truth_MC->Scale(MC_to_Data_Scaling); W_mec_2p_truth_MC->SetLineColor(kGreen); W_mec_2p_cuts_stacked->Add(W_mec_2p_truth_MC);
  xb_mec_2p_truth_MC->SetLineColor(kGreen); xb_mec_2p_truth_MC->Sumw2(); xb_mec_2p_truth_MC->Scale(MC_to_Data_Scaling); xb_mec_2p_truth_MC->SetLineColor(kGreen); xb_mec_2p_cuts_stacked->Add(xb_mec_2p_truth_MC);
  p_proton_mec_2p_truth_MC->SetLineColor(kGreen); p_proton_mec_2p_truth_MC->Sumw2(); p_proton_mec_2p_truth_MC->Scale(MC_to_Data_Scaling); p_proton_mec_2p_truth_MC->SetLineColor(kGreen); p_proton_mec_2p_cuts_stacked->Add(p_proton_mec_2p_truth_MC);

  //MEC proton-neutron (1p1n)


  //RES single charged pion (0p1pi)
  W_res_0p1pi_cuts_MC->SetLineColor(kBlue); W_res_0p1pi_cuts_MC->Sumw2(); W_res_0p1pi_cuts_MC->Scale(MC_to_Data_Scaling); W_res_0p1pi_cuts_MC->SetLineColor(kBlue); W_res_0p1pi_cuts_stacked->Add(W_res_0p1pi_cuts_MC);
  xb_res_0p1pi_cuts_MC->SetLineColor(kBlue); xb_res_0p1pi_cuts_MC->Sumw2(); xb_res_0p1pi_cuts_MC->Scale(MC_to_Data_Scaling); xb_res_0p1pi_cuts_MC->SetLineColor(kBlue); xb_res_0p1pi_cuts_stacked->Add(xb_res_0p1pi_cuts_MC);
  p_pi_res_0p1pi_cuts_MC->SetLineColor(kBlue); p_pi_res_0p1pi_cuts_MC->Sumw2(); p_pi_res_0p1pi_cuts_MC->Scale(MC_to_Data_Scaling); p_pi_res_0p1pi_cuts_MC->SetLineColor(kBlue); p_pi_res_0p1pi_cuts_stacked->Add(p_pi_res_0p1pi_cuts_MC);
  W_res_0p1pi_truth_MC->SetLineColor(kGreen); W_res_0p1pi_truth_MC->Sumw2(); W_res_0p1pi_truth_MC->Scale(MC_to_Data_Scaling); W_res_0p1pi_truth_MC->SetLineColor(kGreen); W_res_0p1pi_cuts_stacked->Add(W_res_0p1pi_truth_MC);
  xb_res_0p1pi_truth_MC->SetLineColor(kGreen); xb_res_0p1pi_truth_MC->Sumw2(); xb_res_0p1pi_truth_MC->Scale(MC_to_Data_Scaling); xb_res_0p1pi_truth_MC->SetLineColor(kGreen); xb_res_0p1pi_cuts_stacked->Add(xb_res_0p1pi_truth_MC);
  p_pi_res_0p1pi_truth_MC->SetLineColor(kGreen); p_pi_res_0p1pi_truth_MC->Sumw2(); p_pi_res_0p1pi_truth_MC->Scale(MC_to_Data_Scaling); p_pi_res_0p1pi_truth_MC->SetLineColor(kGreen); p_pi_res_0p1pi_cuts_stacked->Add(p_pi_res_0p1pi_truth_MC);

  //RES pion-proton (1p1pi)
  W_res_1p1pi_cuts_MC->SetLineColor(kBlue); W_res_1p1pi_cuts_MC->Sumw2(); W_res_1p1pi_cuts_MC->Scale(MC_to_Data_Scaling); W_res_1p1pi_cuts_MC->SetLineColor(kBlue); W_res_1p1pi_cuts_stacked->Add(W_res_1p1pi_cuts_MC);
  xb_res_1p1pi_cuts_MC->SetLineColor(kBlue); xb_res_1p1pi_cuts_MC->Sumw2(); xb_res_1p1pi_cuts_MC->Scale(MC_to_Data_Scaling); xb_res_1p1pi_cuts_MC->SetLineColor(kBlue); xb_res_1p1pi_cuts_stacked->Add(xb_res_1p1pi_cuts_MC);
  p_pi_res_1p1pi_cuts_MC->SetLineColor(kBlue); p_pi_res_1p1pi_cuts_MC->Sumw2(); p_pi_res_1p1pi_cuts_MC->Scale(MC_to_Data_Scaling); p_pi_res_1p1pi_cuts_MC->SetLineColor(kBlue); p_pi_res_1p1pi_cuts_stacked->Add(p_pi_res_1p1pi_cuts_MC);
  p_proton_res_1p1pi_cuts_MC->SetLineColor(kBlue); p_proton_res_1p1pi_cuts_MC->Sumw2(); p_proton_res_1p1pi_cuts_MC->Scale(MC_to_Data_Scaling); p_proton_res_1p1pi_cuts_MC->SetLineColor(kBlue); p_proton_res_1p1pi_cuts_stacked->Add(p_proton_res_1p1pi_cuts_MC);
  W_res_1p1pi_truth_MC->SetLineColor(kGreen); W_res_1p1pi_truth_MC->Sumw2(); W_res_1p1pi_truth_MC->Scale(MC_to_Data_Scaling); W_res_1p1pi_truth_MC->SetLineColor(kGreen); W_res_1p1pi_cuts_stacked->Add(W_res_1p1pi_truth_MC);
  xb_res_1p1pi_truth_MC->SetLineColor(kGreen); xb_res_1p1pi_truth_MC->Sumw2(); xb_res_1p1pi_truth_MC->Scale(MC_to_Data_Scaling); xb_res_1p1pi_truth_MC->SetLineColor(kGreen); xb_res_1p1pi_cuts_stacked->Add(xb_res_1p1pi_truth_MC);
  p_pi_res_1p1pi_truth_MC->SetLineColor(kGreen); p_pi_res_1p1pi_truth_MC->Sumw2(); p_pi_res_1p1pi_truth_MC->Scale(MC_to_Data_Scaling); p_pi_res_1p1pi_truth_MC->SetLineColor(kGreen); p_pi_res_1p1pi_cuts_stacked->Add(p_pi_res_1p1pi_truth_MC);
  p_proton_res_1p1pi_truth_MC->SetLineColor(kGreen); p_proton_res_1p1pi_truth_MC->Sumw2(); p_proton_res_1p1pi_truth_MC->Scale(MC_to_Data_Scaling); p_proton_res_1p1pi_truth_MC->SetLineColor(kGreen); p_proton_res_1p1pi_cuts_stacked->Add(p_proton_res_1p1pi_truth_MC);

  //RES pion-neutron (1n1pi)


  //RES pion-pion-neutron (1n2pi)


  
  //Print out "inclusive" event rates
  outfile_byProcess << "Inclusive event scales and rates" << std::endl;
  outfile_byProcess << "Good data electron events: " << counter_data << std::endl;
  outfile_byProcess << "Good simulation electron events: " << counter_MC << std::endl;
  outfile_byProcess << "MC to Data scaling: " << MC_to_Data_Scaling << std::endl;
  outfile_byProcess << "Actual simulation QE events: " << counter_qe_MC << std::endl;
  outfile_byProcess << "Expected simulation QE events: " << counter_qe_MC*MC_to_Data_Scaling << std::endl;
  outfile_byProcess << "Actual simulation MEC events: " << counter_mec_MC << std::endl;
  outfile_byProcess << "Expected simulation MEC events: " << counter_mec_MC*MC_to_Data_Scaling << std::endl;
  outfile_byProcess << "Actual simulation RES events: " << counter_res_MC << std::endl;
  outfile_byProcess << "Expected simulation RES events: " << counter_res_MC*MC_to_Data_Scaling << std::endl;
  outfile_byProcess << "Actual simulation DIS events: " << counter_dis_MC << std::endl;
  outfile_byProcess << "Expected simulation DIS events: " << counter_dis_MC*MC_to_Data_Scaling << std::endl;
  //Print out specific "exclusive" topology event rates
  //QE single proton
  outfile_byProcess << "QE Single Proton rates" << std::endl;
  outfile_byProcess << "Data with cuts: " << counter_qe_1p0pi_cuts_data << std::endl;
  outfile_byProcess << "GENIE Truth: " << counter_qe_1p0pi_truth_MC << std::endl;
  outfile_byProcess << "GENIE Truth, Scaled: " << counter_qe_1p0pi_truth_MC*MC_to_Data_Scaling << std::endl;
  outfile_byProcess << "GENIE with cuts: " << counter_qe_1p0pi_cuts_MC << std::endl;
  outfile_byProcess << "GENIE with cuts, Scaled: " << counter_qe_1p0pi_cuts_MC*MC_to_Data_Scaling << std::endl;
/*  //QE single neutron
  outfile_byProcess << "QE Single Neutron rates" << std::endl;
  outfile_byProcess << "Data with cuts: " << counter_qe_1p0pi_cuts_data << std::endl;
  outfile_byProcess << "GENIE Truth: " << counter_qe_1p0pi_truth_MC << std::endl;
  outfile_byProcess << "GENIE Truth, Scaled: " << counter_qe_1p0pi_truth_MC*MC_to_Data_Scaling << std::endl;
  outfile_byProcess << "GENIE with cuts: " << counter_qe_1p0pi_cuts_MC << std::endl;
  outfile_byProcess << "GENIE with cuts, Scaled: " << counter_qe_1p0pi_cuts_MC*MC_to_Data_Scaling << std::endl;
*/
  //MEC single proton
  outfile_byProcess << "MEC Single Proton rates" << std::endl;
  outfile_byProcess << "Data with cuts: " << counter_mec_1p0pi_cuts_data << std::endl;
  outfile_byProcess << "GENIE Truth: " << counter_mec_1p0pi_truth_MC << std::endl;
  outfile_byProcess << "GENIE Truth, Scaled: " << counter_mec_1p0pi_truth_MC*MC_to_Data_Scaling << std::endl;
  outfile_byProcess << "GENIE with cuts: " << counter_mec_1p0pi_cuts_MC << std::endl;
  outfile_byProcess << "GENIE with cuts, Scaled: " << counter_mec_1p0pi_cuts_MC*MC_to_Data_Scaling << std::endl;
  //MEC two proton (proton-proton) 
  outfile_byProcess << "MEC Two Proton rates" << std::endl;
  outfile_byProcess << "Data with cuts: " << counter_mec_2p0pi_cuts_data << std::endl;
  outfile_byProcess << "GENIE Truth: " << counter_mec_2p0pi_truth_MC << std::endl;
  outfile_byProcess << "GENIE Truth, Scaled: " << counter_mec_2p0pi_truth_MC*MC_to_Data_Scaling << std::endl;
  outfile_byProcess << "GENIE with cuts: " << counter_mec_2p0pi_cuts_MC << std::endl;
  outfile_byProcess << "GENIE with cuts, Scaled: " << counter_mec_2p0pi_cuts_MC*MC_to_Data_Scaling << std::endl;
/*  //MEC proton-neutron (1p1n)
  outfile_byProcess << "MEC Proton-Neutron rates" << std::endl;
  outfile_byProcess << "Data with cuts: " << counter_qe_1p0pi_cuts_data << std::endl;
  outfile_byProcess << "GENIE Truth: " << counter_qe_1p0pi_truth_MC << std::endl;
  outfile_byProcess << "GENIE Truth, Scaled: " << counter_qe_1p0pi_truth_MC*MC_to_Data_Scaling << std::endl;
  outfile_byProcess << "GENIE with cuts: " << counter_qe_1p0pi_cuts_MC << std::endl;
  outfile_byProcess << "GENIE with cuts, Scaled: " << counter_qe_1p0pi_cuts_MC*MC_to_Data_Scaling << std::endl;
*/  //RES single charged pion (0p1pi)
  outfile_byProcess << "MEC Single Charged Pion (0p1pi) rates" << std::endl;
  outfile_byProcess << "Data with cuts: " << counter_res_0p1pi_cuts_data << std::endl;
  outfile_byProcess << "GENIE Truth: " << counter_res_0p1pi_truth_MC << std::endl;
  outfile_byProcess << "GENIE Truth, Scaled: " << counter_res_0p1pi_truth_MC*MC_to_Data_Scaling << std::endl;
  outfile_byProcess << "GENIE with cuts: " << counter_res_0p1pi_cuts_MC << std::endl;
  outfile_byProcess << "GENIE with cuts, Scaled: " << counter_res_0p1pi_cuts_MC*MC_to_Data_Scaling << std::endl;
  //RES single charged pion with single proton (1p1pi)
  outfile_byProcess << "MEC Proton-Pion rates" << std::endl;
  outfile_byProcess << "Data with cuts: " << counter_res_1p1pi_cuts_data << std::endl;
  outfile_byProcess << "GENIE Truth: " << counter_res_1p1pi_truth_MC << std::endl;
  outfile_byProcess << "GENIE Truth, Scaled: " << counter_res_1p1pi_truth_MC*MC_to_Data_Scaling << std::endl;
  outfile_byProcess << "GENIE with cuts: " << counter_res_1p1pi_cuts_MC << std::endl;
  outfile_byProcess << "GENIE with cuts, Scaled: " << counter_res_1p1pi_cuts_MC*MC_to_Data_Scaling << std::endl;
/*  //RES pion-neutron (1n1pi)
  outfile_byProcess << "MEC Single Neutron-Pion rates" << std::endl;
  outfile_byProcess << "Data with cuts: " << counter_qe_1p0pi_cuts_data << std::endl;
  outfile_byProcess << "GENIE Truth: " << counter_qe_1p0pi_truth_MC << std::endl;
  outfile_byProcess << "GENIE Truth, Scaled: " << counter_qe_1p0pi_truth_MC*MC_to_Data_Scaling << std::endl;
  outfile_byProcess << "GENIE with cuts: " << counter_qe_1p0pi_cuts_MC << std::endl;
  outfile_byProcess << "GENIE with cuts, Scaled: " << counter_qe_1p0pi_cuts_MC*MC_to_Data_Scaling << std::endl;
  //RES pion-pion-neutron (1n2pi)
  outfile_byProcess << "MEC Sinle Neutron Two Pion rates" << std::endl;
  outfile_byProcess << "Data with cuts: " << counter_qe_1p0pi_cuts_data << std::endl;
  outfile_byProcess << "GENIE Truth: " << counter_qe_1p0pi_truth_MC << std::endl;
  outfile_byProcess << "GENIE Truth, Scaled: " << counter_qe_1p0pi_truth_MC*MC_to_Data_Scaling << std::endl;
  outfile_byProcess << "GENIE with cuts: " << counter_qe_1p0pi_cuts_MC << std::endl;
  outfile_byProcess << "GENIE with cuts, Scaled: " << counter_qe_1p0pi_cuts_MC*MC_to_Data_Scaling << std::endl;
*/


  TFile *outFile = new TFile(outputFile + ".root","RECREATE");

  TCanvas *c1 = new TCanvas("c1","c1",6000,4000);
  TString fileName = outputFile + ".pdf";

  //Set some style options
  gStyle->SetStatY(0.9);
  gStyle->SetStatX(0.9);
  //gStyle->SetOptStat(0);
  gStyle->SetTitleX(0.5);
  gStyle->SetTitleAlign(23);
  

  //First page of the PDF
  //Divide page into 2x3 grid
  c1->Divide(2,3);
  //Top left
  c1->cd(1);
  el_vz_h_stacked->Draw("hist && nostack");
  c1->cd(1)->BuildLegend(0.7,0.7,0.9,0.9);
  //Top right
  c1->cd(2);
  c1->cd(2)->SetLogx();
  c1->cd(2)->SetLogy();
  q2_stacked->Draw("hist && nostack");
  c1->cd(2)->BuildLegend(0.7,0.7,0.9,0.9);
  //Middle left
  c1->cd(3);
  c1->cd(3)->SetLogx();
  q2_theta_data->Draw("colz"); //q2_theta_data->Write();
  //Middle right
  c1->cd(4);
  c1->cd(4)->SetLogx();
  q2_theta_MC->Draw("colz"); //q2_theta_MC->Write();
  //Bottom left
  c1->cd(5);
  theta_mom_data->Draw("colz"); //theta_mom_data->Write();
  //Bottom right
  c1->cd(6);
  theta_mom_MC->Draw("colz"); //theta_mom_MC->Write();
  //Start print of page to file
  c1->Print(fileName+"(");
  c1->Clear();
  

  //Second page
  //Divide page into 2x3 grid
  c1->Divide(2,3);
  //Top left
  c1->cd(1);
  W_stacked->Draw("hist && nostack");
  c1->cd(1)->BuildLegend(0.7,0.7,0.9,0.9);
  //Top right
  c1->cd(2);
  c1->cd(2)->SetLogx();
  c1->cd(2)->SetLogy();
  q2_stacked->Draw("hist && nostack");
  c1->cd(2)->BuildLegend(0.7,0.7,0.9,0.9);
  //Middle right
  c1->cd(3);
  c1->cd(3)->SetLogy();
  q2_W_data->Draw("colz");
  //Middle left
  c1->cd(4);
  c1->cd(4)->SetLogy();
  q2_W_MC->Draw("colz");
  //Bottom right
  c1->cd(5);
  W_theta_data->Draw("colz");
  //Bottom left
  c1->cd(6);
  W_theta_MC->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //Divide page into 2x3 grid
  c1->Divide(2,3);
  //Top Left
  c1->cd(1);
  c1->cd(1)->SetLogz();
  q0_vs_q3_data->Draw("colz");
  //Top Right
  c1->cd(2);
  c1->cd(2)->SetLogz();
  q0_vs_q3_MC->Draw("colz");
  //Middle Left
  c1->cd(3);
  c1->cd(3)->SetLogz();
  q0_vs_q3_1pXn0pi_data->Draw("colz");
  //Middle Right
  c1->cd(4);
  c1->cd(4)->SetLogz();
  q0_vs_q3_1pXn0pi_MC->Draw("colz");
  //Bottom Left
  c1->cd(5);
  c1->cd(5)->SetLogz();
  q0_vs_q3_2pXn0pi_data->Draw("colz");
  //Bottom Right
  c1->cd(6);
  c1->cd(6)->SetLogz();
  q0_vs_q3_2pXn0pi_MC->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();


  //Divide page into 2x2 grid
  c1->Divide(2,2);
  //Top left
  c1->cd(1);
  W_bySector_data_stacked->Draw("hist && nostack");
  c1->cd(1)->BuildLegend(0.7,0.7,0.9,0.9);
  //Top right
  c1->cd(2);
  W_bySector_MC_stacked->Draw("hist && nostack");
  c1->cd(2)->BuildLegend(0.7,0.7,0.9,0.9);
  //Bottom right
  c1->cd(3);
  c1->cd(3)->SetLogx();
  c1->cd(3)->SetLogy();
  q2_bySector_data_stacked->Draw("hist && nostack");
  c1->cd(2)->BuildLegend(0.7,0.7,0.9,0.9);
  //Bottom left
  c1->cd(4);
  c1->cd(4)->SetLogx();
  c1->cd(4)->SetLogy();
  q2_bySector_MC_stacked->Draw("hist && nostack");
  c1->cd(2)->BuildLegend(0.7,0.7,0.9,0.9);
  //Print page to file
  c1->Print(fileName);
  c1->Clear();


  //Divide page into 2x6 grid
  c1->Divide(2,6);
  c1->cd(1);
  W_theta_S1_data->Draw("colz");
  c1->cd(2);
  W_theta_S1_MC->Draw("colz");
  c1->cd(3);
  W_theta_S2_data->Draw("colz");
  c1->cd(4);
  W_theta_S2_MC->Draw("colz");
  c1->cd(5);
  W_theta_S3_data->Draw("colz");
  c1->cd(6);
  W_theta_S3_MC->Draw("colz");
  c1->cd(7);
  W_theta_S4_data->Draw("colz");
  c1->cd(8);
  W_theta_S4_MC->Draw("colz");
  c1->cd(9);
  W_theta_S5_data->Draw("colz");
  c1->cd(10);
  W_theta_S5_MC->Draw("colz");
  c1->cd(11);
  W_theta_S6_data->Draw("colz");
  c1->cd(12);
  W_theta_S6_MC->Draw("colz");  
  //Print page to file
  c1->Print(fileName);
  c1->Clear();


  //Third page
  //Divide page into 2x2 grid
  c1->Divide(2,3);
  //Top left
  c1->cd(1);
  omega_stacked->Draw("hist && nostack");
  c1->cd(1)->BuildLegend(0.7,0.7,0.9,0.9);
  //Middle left
  c1->cd(3);
  energy_transfer_0_stacked->Draw("hist && nostack");
  c1->cd(3)->BuildLegend(0.7,0.7,0.9,0.9);
  //Middle right
  c1->cd(4);
  energy_transfer_1_stacked->Draw("hist && nostack");
  c1->cd(4)->BuildLegend(0.7,0.7,0.9,0.9);
  //Middle right
  c1->cd(5);
  energy_transfer_2_stacked->Draw("hist && nostack");
  c1->cd(5)->BuildLegend(0.7,0.7,0.9,0.9);
  //Middle left
  c1->cd(6);
  energy_transfer_3_stacked->Draw("hist && nostack");
  c1->cd(6)->BuildLegend(0.1,0.7,0.25,0.9);
  //Bottom left
  //Print page to file
  c1->Print(fileName);
  c1->Clear();


  //Fourth page
  //Divide page into 2x2 grid
  c1->Divide(2,3);
  //Top left
  c1->cd(1);
  omega_theta_data->Draw("colz");
  //Top right
  c1->cd(2);
  omega_theta_MC->Draw("colz");
  //Middle left
  c1->cd(3);
  energy_transfer_0_stacked->Draw("hist && nostack");
  c1->cd(3)->BuildLegend(0.7,0.7,0.9,0.9);
  //Middle right
  c1->cd(4);
  energy_transfer_1_stacked->Draw("hist && nostack");
  c1->cd(4)->BuildLegend(0.7,0.7,0.9,0.9);
  //Middle right
  c1->cd(5);
  energy_transfer_2_stacked->Draw("hist && nostack");
  c1->cd(5)->BuildLegend(0.7,0.7,0.9,0.9);
  //Middle left
  c1->cd(6);
  energy_transfer_3_stacked->Draw("hist && nostack");
  c1->cd(6)->BuildLegend(0.1,0.7,0.25,0.9);
  //Bottom left
  //Print page to file
  c1->Print(fileName);
  c1->Clear();


  //Fifth page
  gStyle->SetStatY(0.9);
  gStyle->SetStatX(0.9);
  //Divide page into 2x3 grid
  c1->Divide(2,3);
  //Top left
  c1->cd(1);
  //c1->cd(1)->SetLogz();
  p_miss_theta_vs_omega_data->Draw("colz");
  //Top right
  c1->cd(2);
  //c1->cd(2)->SetLogz();
  p_miss_theta_vs_omega_MC->Draw("colz");
  gStyle->SetStatY(0.9);
  gStyle->SetStatX(0.9);
  //Middle left
  c1->cd(3);
  p_proton_stacked->Draw("hist && nostack");
  c1->cd(3)->BuildLegend(0.7,0.7,0.9,0.9);
  //Middle right
  c1->cd(4);
  p_miss_stacked->Draw("hist && nostack");
  c1->cd(4)->BuildLegend(0.7,0.7,0.9,0.9);
  //Bottom left
  gStyle->SetStatY(0.9);
  gStyle->SetStatX(0.1);
  c1->cd(5);
  p_miss_theta_stacked->Draw("hist && nostack");
  c1->cd(5)->BuildLegend(0.3,0.7,0.6,0.9);
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //Divide page into 2x2 grid
  c1->Divide(2,2);
  //Top left
  c1->cd(1);
  //c1->cd(1)->SetLogz();
  p_chi2_stacked->Draw("hist && nostack");
  c1->cd(1)->BuildLegend(0.7,0.7,0.9,0.9);
  //Top right
  c1->cd(2);
  //c1->cd(2)->SetLogz();
  n_chi2_stacked->Draw("hist && nostack");
  c1->cd(2)->BuildLegend(0.7,0.7,0.9,0.9);
  //Bottom left
  c1->cd(3);
  pip_chi2_stacked->Draw("hist && nostack");
  c1->cd(3)->BuildLegend(0.7,0.7,0.9,0.9);
  //Bottom right
  c1->cd(4);
  pim_chi2_stacked->Draw("hist && nostack");
  c1->cd(4)->BuildLegend(0.7,0.7,0.9,0.9);
  //Print page to file
  c1->Print(fileName);
  c1->Clear();


  //Sixth page
  //Divide page into 2x4 grid
  gStyle->SetStatY(0.9);
  gStyle->SetStatX(0.9);
  c1->Divide(2,4);
  //Top left
  c1->cd(1);
  c1->cd(1)->SetLogy();//Make the multiplcity plots logscale
  mult_p_stacked->Draw("hist && nostack");
  c1->cd(1)->BuildLegend(0.6,0.6,0.9,0.9);
  //Top right
  c1->cd(2);
  c1->cd(2)->SetLogy();//Make the multiplcity plots logscale
  mult_n_stacked->Draw("hist && nostack");
  c1->cd(2)->BuildLegend(0.6,0.6,0.9,0.9);
  //Middle left
  c1->cd(3);
  c1->cd(3)->SetLogy();
  mult_pip_stacked->Draw("hist && nostack");
  c1->cd(3)->BuildLegend(0.6,0.6,0.9,0.9);
  //Middle right
  c1->cd(4);
  c1->cd(4)->SetLogy();
  mult_pim_stacked->Draw("hist && nostack");
  c1->cd(4)->BuildLegend(0.6,0.6,0.9,0.9);
  //Bottom left
  c1->cd(5);
  mult_p_pis_data->Draw("colz && text");
  c1->cd(5)->SetLogz();
  //Bottom right
  c1->cd(6);
  mult_p_pis_MC->Draw("colz && text");
  c1->cd(6)->SetLogz();
  //Bottom left
  c1->cd(7);
  mult_p_n_data->Draw("colz && text");
  c1->cd(7)->SetLogz();
  //Bottom right
  c1->cd(8);
  mult_p_n_MC->Draw("colz && text");
  c1->cd(8)->SetLogz();
  //Print page to file
  c1->Print(fileName);
  c1->Clear();


  //Seventh page
  //Divide page into 2x2 grid
  c1->Divide(2,2);
  //Top left
  c1->cd(1);
  q2_stacked->Draw("hist && nostack");
  //c1->cd(1)->SetLogy();
  c1->cd(1)->BuildLegend(0.6,0.6,0.9,0.9);
  //Top right
  c1->cd(2);
  W_stacked->Draw("hist && nostack");
  //c1->cd(2)->SetLogy();
  c1->cd(2)->BuildLegend(0.1,0.6,0.35,0.9);
  //Bottom left
  c1->cd(3);
  xb_stacked->Draw("hist && nostack");
  //c1->cd(3)->SetLogy();
  c1->cd(3)->BuildLegend(0.6,0.6,0.9,0.9);
  //Bottom right
  c1->cd(4);
  y_stacked->Draw("hist && nostack");
  //c1->cd(4)->SetLogy();
  c1->cd(4)->BuildLegend(0.1,0.6,0.35,0.9);
  //Print page to file
  c1->Print(fileName);
  c1->Clear();


  //Seventh page (part deux, making everything log)
  //Divide page into 2x2 grid
  c1->Divide(2,2);
  //Top left
  c1->cd(1);
  q2_stacked->Draw("hist && nostack");
  c1->cd(1)->SetLogy();
  c1->cd(1)->BuildLegend(0.6,0.6,0.9,0.9);
  //Top right
  c1->cd(2);
  W_stacked->Draw("hist && nostack");
  c1->cd(2)->SetLogy();
  c1->cd(2)->BuildLegend(0.1,0.6,0.35,0.9);
  //Bottom left
  c1->cd(3);
  xb_stacked->Draw("hist && nostack");
  c1->cd(3)->SetLogy();
  c1->cd(3)->BuildLegend(0.6,0.6,0.9,0.9);
  //Bottom right
  c1->cd(4);
  y_stacked->Draw("hist && nostack");
  c1->cd(4)->SetLogy();
  c1->cd(4)->BuildLegend(0.1,0.6,0.35,0.9);
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //Divide page into 2x2 grid
  c1->Divide(2,2);
  //Top left
  c1->cd(1);
  q2_byChannel_stacked->Draw("hist && nostack");
  //c1->cd(1)->SetLogy();
  c1->cd(1)->BuildLegend(0.6,0.6,0.9,0.9);
  //Top right
  c1->cd(2);
  W_byChannel_stacked->Draw("hist && nostack");
  //c1->cd(2)->SetLogy();
  c1->cd(2)->BuildLegend(0.6,0.6,0.9,0.9);
  //Bottom left
  c1->cd(3);
  totrecoe_byChannel_stacked->Draw("hist && nostack");
  //c1->cd(3)->SetLogy();
  c1->cd(3)->BuildLegend(0.6,0.6,0.9,0.9);
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //Divide page into 2x2 grid
  c1->Divide(2,2);
  //Top left
  c1->cd(1);
  q2_byChannel_stacked->Draw("hist && nostack");
  c1->cd(1)->SetLogy();
  c1->cd(1)->BuildLegend(0.7,0.6,0.9,0.9);
  //Top right
  c1->cd(2);
  W_byChannel_stacked->Draw("hist && nostack");
  c1->cd(2)->SetLogy();
  c1->cd(2)->BuildLegend(0.7,0.6,0.9,0.9);
  //Bottom left
  c1->cd(3);
  totrecoe_byChannel_stacked->Draw("hist && nostack");
  c1->cd(3)->SetLogy();
  c1->cd(3)->BuildLegend(0.7,0.6,0.9,0.9);
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  
  //Eighth page
  //Divide page into 2x3 grid
  c1->Divide(2,3);
  //Top left
  c1->cd(1);
  W_stacked->Draw("hist && nostack");
  c1->cd(1)->BuildLegend(0.6,0.6,0.9,0.9);
  //Top right
  c1->cd(2);
  resid_MC->Draw("hist");
  c1->cd(2)->BuildLegend(0.6,0.6,0.9,0.9);
  //Middle left
  c1->cd(3);
  W_resonances_1p2_to_2p1_MC->Draw("hist");
  c1->cd(3)->BuildLegend(0.6,0.6,0.9,0.9);
  //Middle right
  c1->cd(4);
  resonances_1p2_to_2p1_MC->Draw("hist");
  c1->cd(4)->BuildLegend(0.6,0.6,0.9,0.9);
  //Bottom left
  c1->cd(5);
  W_resid_1p2_to_2p1_MC->Draw("colz");
  c1->cd(5)->BuildLegend(0.6,0.6,0.9,0.9); 
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  
  //Ninth page
  //Divide page into 2x2 grid
  c1->Divide(2,2);
  //Top left
  c1->cd(1);
  W_qe_1p_cuts_stacked->Draw("hist && nostack");
  c1->cd(1)->BuildLegend(0.1,0.6,0.35,0.9);
  //Top right
  c1->cd(2);
  xb_qe_1p_cuts_stacked->Draw("hist && nostack");
  c1->cd(2)->BuildLegend(0.1,0.6,0.35,0.9);
  //Bottom left
  c1->cd(3);
  p_proton_qe_1p_cuts_stacked->Draw("hist && nostack");
  c1->cd(3)->BuildLegend(0.65,0.65,0.9,0.9);
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  
  //Tenth page
  //Divide page into 2x2 grid
  c1->Divide(2,2);
  //Top left
  c1->cd(1);
  W_mec_1p_cuts_stacked->Draw("hist && nostack");
  c1->cd(1)->BuildLegend(0.1,0.6,0.35,0.9);
  //Top right
  c1->cd(2);
  xb_mec_1p_cuts_stacked->Draw("hist && nostack");
  c1->cd(2)->BuildLegend(0.1,0.6,0.35,0.9);
  //Bottom left
  c1->cd(3);
  p_proton_mec_1p_cuts_stacked->Draw("hist && nostack");
  c1->cd(3)->BuildLegend(0.65,0.65,0.9,0.9);
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  
  //Eleventh page
  //Divide page into 2x2 grid
  c1->Divide(2,2);
  //Top left
  c1->cd(1);
  W_mec_2p_cuts_stacked->Draw("hist && nostack");
  c1->cd(1)->BuildLegend(0.1,0.6,0.35,0.9);
  //Top right
  c1->cd(2);
  xb_mec_2p_cuts_stacked->Draw("hist && nostack");
  c1->cd(2)->BuildLegend(0.1,0.6,0.35,0.9);
  //Bottom left
  c1->cd(3);
  p_proton_mec_2p_cuts_stacked->Draw("hist && nostack");
  c1->cd(3)->BuildLegend(0.65,0.65,0.9,0.9);
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  
  //Twelth page
  //Divide page into 2x2 grid
  c1->Divide(2,2);
  //Top left
  c1->cd(1);
  W_res_0p1pi_cuts_stacked->Draw("hist && nostack");
  c1->cd(1)->BuildLegend(0.65,0.65,0.9,0.9);
  //Top right
  c1->cd(2);
  xb_res_0p1pi_cuts_stacked->Draw("hist && nostack");
  c1->cd(2)->BuildLegend(0.65,0.65,0.9,0.9);
  //Bottom left
  c1->cd(3);
  p_pi_res_0p1pi_cuts_stacked->Draw("hist && nostack");
  c1->cd(3)->BuildLegend(0.65,0.65,0.9,0.9);
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  
  //Thirteenth page
  //Divide page into 2x2 grid
  c1->Divide(2,2);
  //Top left
  c1->cd(1);
  W_res_1p1pi_cuts_stacked->Draw("hist && nostack");
  c1->cd(1)->BuildLegend(0.65,0.65,0.9,0.9);
  //Top right
  c1->cd(2);
  xb_res_1p1pi_cuts_stacked->Draw("hist && nostack");
  c1->cd(2)->BuildLegend(0.65,0.65,0.9,0.9);
  //Bottom left
  c1->cd(3);
  p_pi_res_1p1pi_cuts_stacked->Draw("hist && nostack");
  c1->cd(3)->BuildLegend(0.65,0.65,0.9,0.9);
  //Bottom right
  c1->cd(4);
  p_proton_res_1p1pi_cuts_stacked->Draw("hist && nostack");
  c1->cd(3)->BuildLegend(0.65,0.65,0.9,0.9);
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  /////////////////////////////////////////////////////////////////////////////////
  //2D Histograms of topologies
  /////////////////////////////////////////////////////////////////////////////////
  //Fourteenth page (0p1pip totrecoe vs. kinematic variables)
  //Divide page into 2x5 grid
  c1->Divide(2,5);
  //Top left
  c1->cd(1);
  totrecoe_W_0p_1pip_data->Draw("colz");
  //Top right
  c1->cd(2);
  totrecoe_W_0p_1pip_MC->Draw("colz");
  //Middle top left
  c1->cd(3);
  totrecoe_q2_0p_1pip_data->Draw("colz");
  //Middle top right
  c1->cd(4);
  totrecoe_q2_0p_1pip_MC->Draw("colz");
  //Middle bottom left
  c1->cd(5);
  totrecoe_xb_0p_1pip_data->Draw("colz");
  //Middle bottom right
  c1->cd(6);
  totrecoe_xb_0p_1pip_MC->Draw("colz");
  //Bottom left
  c1->cd(7);
  totrecoe_y_0p_1pip_data->Draw("colz");
  //Bottom right
  c1->cd(8);
  totrecoe_y_0p_1pip_MC->Draw("colz");
  //Very bottom left
  c1->cd(9);
  totrecoe_omega_0p_1pip_data->Draw("colz");
  //Very bottom right
  c1->cd(10);
  totrecoe_omega_0p_1pip_MC->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //Divide page into 2x5 grid
  c1->Divide(2,5);
  //Top left
  c1->cd(1);
  nuc_inv_mass_W_0p_1pip_data->Draw("colz");
  //Top right
  c1->cd(2);
  nuc_inv_mass_W_0p_1pip_MC->Draw("colz");
  //Middle top left
  c1->cd(3);
  nuc_inv_mass_q2_0p_1pip_data->Draw("colz");
  //Middle top right
  c1->cd(4);
  nuc_inv_mass_q2_0p_1pip_MC->Draw("colz");
  //Middle bottom left
  c1->cd(5);
  nuc_inv_mass_xb_0p_1pip_data->Draw("colz");
  //Middle bottom right
  c1->cd(6);
  nuc_inv_mass_xb_0p_1pip_MC->Draw("colz");
  //Bottom left
  c1->cd(7);
  nuc_inv_mass_y_0p_1pip_data->Draw("colz");
  //Bottom right
  c1->cd(8);
  nuc_inv_mass_y_0p_1pip_MC->Draw("colz");
  //Very bottom left
  c1->cd(9);
  nuc_inv_mass_omega_0p_1pip_data->Draw("colz");
  //Very bottom right
  c1->cd(10);
  nuc_inv_mass_omega_0p_1pip_MC->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //Fifteenth page (0p1pip W vs. remaining kinematic variables)
  //Divide page into 2x3 grid
  c1->Divide(2,3);
  //Top left
  c1->cd(1);
  W_q2_0p_1pip_data->Draw("colz");
  //Top right
  c1->cd(2);
  W_q2_0p_1pip_MC->Draw("colz");
  //Middle left
  c1->cd(3);
  W_xb_0p_1pip_data->Draw("colz");
  //Middle right
  c1->cd(4);
  W_xb_0p_1pip_MC->Draw("colz");
  //Bottom left
  c1->cd(5);
  W_y_0p_1pip_data->Draw("colz");
  //Bottom right
  c1->cd(6);
  W_y_0p_1pip_MC->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //Sixteenth page (0p1pip totrecoe vs. kinematic variables)
  //Divide page into 2x2 grid
  c1->Divide(2,2);
  //Top left
  c1->cd(1);
  q2_xb_0p_1pip_data->Draw("colz");
  //Top right
  c1->cd(2);
  q2_xb_0p_1pip_MC->Draw("colz");
  //Bottom left
  c1->cd(3);
  q2_y_0p_1pip_data->Draw("colz");
  //Bottom right
  c1->cd(4);
  q2_y_0p_1pip_MC->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();
  
  //Seventeenth page (0p1pip xb vs y)
  //Divide page into 2x2 grid
  c1->Divide(2,2);
  //Left
  c1->cd(1);
  y_xb_0p_1pip_data->Draw("colz");
  //Right
  c1->cd(2);
  y_xb_0p_1pip_MC->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //Eighteenth page (0p1pip totrecoe vs. kinematic variables)
  //Divide page into 2x5 grid
  c1->Divide(2,5);
  //Top left
  c1->cd(1);
  totrecoe_W_0p_1pim_data->Draw("colz");
  //Top right
  c1->cd(2);
  totrecoe_W_0p_1pim_MC->Draw("colz");
  //Middle top left
  c1->cd(3);
  totrecoe_q2_0p_1pim_data->Draw("colz");
  //Middle top right
  c1->cd(4);
  totrecoe_q2_0p_1pim_MC->Draw("colz");
  //Middle bottom left
  c1->cd(5);
  totrecoe_xb_0p_1pim_data->Draw("colz");
  //Middle bottom right
  c1->cd(6);
  totrecoe_xb_0p_1pim_MC->Draw("colz");
  //Bottom left
  c1->cd(7);
  totrecoe_y_0p_1pim_data->Draw("colz");
  //Bottom right
  c1->cd(8);
  totrecoe_y_0p_1pim_MC->Draw("colz");
  //Very bottom left
  c1->cd(9);
  totrecoe_omega_0p_1pim_data->Draw("colz");
  //Very bottom right
  c1->cd(10);
  totrecoe_omega_0p_1pim_MC->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //Divide page into 2x5 grid
  c1->Divide(2,5);
  //Top left
  c1->cd(1);
  nuc_inv_mass_W_0p_1pim_data->Draw("colz");
  //Top right
  c1->cd(2);
  nuc_inv_mass_W_0p_1pim_MC->Draw("colz");
  //Middle top left
  c1->cd(3);
  nuc_inv_mass_q2_0p_1pim_data->Draw("colz");
  //Middle top right
  c1->cd(4);
  nuc_inv_mass_q2_0p_1pim_MC->Draw("colz");
  //Middle bottom left
  c1->cd(5);
  nuc_inv_mass_xb_0p_1pim_data->Draw("colz");
  //Middle bottom right
  c1->cd(6);
  nuc_inv_mass_xb_0p_1pim_MC->Draw("colz");
  //Bottom left
  c1->cd(7);
  nuc_inv_mass_y_0p_1pim_data->Draw("colz");
  //Bottom right
  c1->cd(8);
  nuc_inv_mass_y_0p_1pim_MC->Draw("colz");
  //Very bottom left
  c1->cd(9);
  nuc_inv_mass_omega_0p_1pim_data->Draw("colz");
  //Very bottom right
  c1->cd(10);
  nuc_inv_mass_omega_0p_1pim_MC->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //Nineteenth page (0p1pip W vs. remaining kinematic variables)
  //Divide page into 2x3 grid
  c1->Divide(2,3);
  //Top left
  c1->cd(1);
  W_q2_0p_1pim_data->Draw("colz");
  //Top right
  c1->cd(2);
  W_q2_0p_1pim_MC->Draw("colz");
  //Middle left
  c1->cd(3);
  W_xb_0p_1pim_data->Draw("colz");
  //Middle right
  c1->cd(4);
  W_xb_0p_1pim_MC->Draw("colz");
  //Bottom left
  c1->cd(5);
  W_y_0p_1pim_data->Draw("colz");
  //Bottom right
  c1->cd(6);
  W_y_0p_1pim_MC->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //Twentieth page (0p1pip totrecoe vs. kinematic variables)
  //Divide page into 2x2 grid
  c1->Divide(2,2);
  //Top left
  c1->cd(1);
  q2_xb_0p_1pim_data->Draw("colz");
  //Top right
  c1->cd(2);
  q2_xb_0p_1pim_MC->Draw("colz");
  //Bottom left
  c1->cd(3);
  q2_y_0p_1pim_data->Draw("colz");
  //Bottom right
  c1->cd(4);
  q2_y_0p_1pim_MC->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();
  
  //Twenty-first page (0p1pip xb vs y)
  //Divide page into 2x2 grid
  c1->Divide(2,2);
  //Left
  c1->cd(1);
  y_xb_0p_1pim_data->Draw("colz");
  //Right
  c1->cd(2);
  y_xb_0p_1pim_MC->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //Twenty-second page (0p1pip totrecoe vs. kinematic variables)
  //Divide page into 2x5 grid
  c1->Divide(2,5);
  //Top left
  c1->cd(1);
  totrecoe_W_1pXn0pi_data->Draw("colz");
  //Top right
  c1->cd(2);
  totrecoe_W_1pXn0pi_MC->Draw("colz");
  //Middle top left
  c1->cd(3);
  totrecoe_q2_1pXn0pi_data->Draw("colz");
  //Middle top right
  c1->cd(4);
  totrecoe_q2_1pXn0pi_MC->Draw("colz");
  //Middle bottom left
  c1->cd(5);
  totrecoe_xb_1pXn0pi_data->Draw("colz");
  //Middle bottom right
  c1->cd(6);
  totrecoe_xb_1pXn0pi_MC->Draw("colz");
  //Bottom left
  c1->cd(7);
  totrecoe_y_1pXn0pi_data->Draw("colz");
  //Bottom right
  c1->cd(8);
  totrecoe_y_1pXn0pi_MC->Draw("colz");
  //Very bottom left
  c1->cd(9);
  totrecoe_omega_1pXn0pi_data->Draw("colz");
  //Very bottom right
  c1->cd(10);
  totrecoe_omega_1pXn0pi_MC->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //Divide page into 2x5 grid
  c1->Divide(2,5);
  //Top left
  c1->cd(1);
  nuc_inv_mass_W_1pXn0pi_data->Draw("colz");
  //Top right
  c1->cd(2);
  nuc_inv_mass_W_1pXn0pi_MC->Draw("colz");
  //Middle top left
  c1->cd(3);
  nuc_inv_mass_q2_1pXn0pi_data->Draw("colz");
  //Middle top right
  c1->cd(4);
  nuc_inv_mass_q2_1pXn0pi_MC->Draw("colz");
  //Middle bottom left
  c1->cd(5);
  nuc_inv_mass_xb_1pXn0pi_data->Draw("colz");
  //Middle bottom right
  c1->cd(6);
  nuc_inv_mass_xb_1pXn0pi_MC->Draw("colz");
  //Bottom left
  c1->cd(7);
  nuc_inv_mass_y_1pXn0pi_data->Draw("colz");
  //Bottom right
  c1->cd(8);
  nuc_inv_mass_y_1pXn0pi_MC->Draw("colz");
  //Very bottom left
  c1->cd(9);
  nuc_inv_mass_omega_1pXn0pi_data->Draw("colz");
  //Very bottom right
  c1->cd(10);
  nuc_inv_mass_omega_1pXn0pi_MC->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //Twenty-third page (0p1pip W vs. remaining kinematic variables)
  //Divide page into 2x3 grid
  c1->Divide(2,3);
  //Top left
  c1->cd(1);
  W_q2_1pXn0pi_data->Draw("colz");
  //Top right
  c1->cd(2);
  W_q2_1pXn0pi_MC->Draw("colz");
  //Middle left
  c1->cd(3);
  W_xb_1pXn0pi_data->Draw("colz");
  //Middle right
  c1->cd(4);
  W_xb_1pXn0pi_MC->Draw("colz");
  //Bottom left
  c1->cd(5);
  W_y_1pXn0pi_data->Draw("colz");
  //Bottom right
  c1->cd(6);
  W_y_1pXn0pi_MC->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //Twenty-fourth page (0p1pip totrecoe vs. kinematic variables)
  //Divide page into 2x2 grid
  c1->Divide(2,2);
  //Top left
  c1->cd(1);
  q2_xb_1pXn0pi_data->Draw("colz");
  //Top right
  c1->cd(2);
  q2_xb_1pXn0pi_MC->Draw("colz");
  //Bottom left
  c1->cd(3);
  q2_y_1pXn0pi_data->Draw("colz");
  //Bottom right
  c1->cd(4);
  q2_y_1pXn0pi_MC->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();
  
  //Twenty-fifth page (0p1pip xb vs y)
  //Divide page into 2x2 grid
  c1->Divide(2,2);
  //Left
  c1->cd(1);
  y_xb_1pXn0pi_data->Draw("colz");
  //Right
  c1->cd(2);
  y_xb_1pXn0pi_MC->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //Twenty-sixth page (0p1pip totrecoe vs. kinematic variables)
  //Divide page into 2x5 grid
  c1->Divide(2,5);
  //Top left
  c1->cd(1);
  totrecoe_W_1pXn1pim_data->Draw("colz");
  //Top right
  c1->cd(2);
  totrecoe_W_1pXn1pim_MC->Draw("colz");
  //Middle top left
  c1->cd(3);
  totrecoe_q2_1pXn1pim_data->Draw("colz");
  //Middle top right
  c1->cd(4);
  totrecoe_q2_1pXn1pim_MC->Draw("colz");
  //Middle bottom left
  c1->cd(5);
  totrecoe_xb_1pXn1pim_data->Draw("colz");
  //Middle bottom right
  c1->cd(6);
  totrecoe_xb_1pXn1pim_MC->Draw("colz");
  //Bottom left
  c1->cd(7);
  totrecoe_y_1pXn1pim_data->Draw("colz");
  //Bottom right
  c1->cd(8);
  totrecoe_y_1pXn1pim_MC->Draw("colz");
  //Very bottom left
  c1->cd(9);
  totrecoe_omega_1pXn1pim_data->Draw("colz");
  //Very bottom right
  c1->cd(10);
  totrecoe_omega_1pXn1pim_MC->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //Divide page into 2x5 grid
  c1->Divide(2,5);
  //Top left
  c1->cd(1);
  nuc_inv_mass_W_1pXn1pim_data->Draw("colz");
  //Top right
  c1->cd(2);
  nuc_inv_mass_W_1pXn1pim_MC->Draw("colz");
  //Middle top left
  c1->cd(3);
  nuc_inv_mass_q2_1pXn1pim_data->Draw("colz");
  //Middle top right
  c1->cd(4);
  nuc_inv_mass_q2_1pXn1pim_MC->Draw("colz");
  //Middle bottom left
  c1->cd(5);
  nuc_inv_mass_xb_1pXn1pim_data->Draw("colz");
  //Middle bottom right
  c1->cd(6);
  nuc_inv_mass_xb_1pXn1pim_MC->Draw("colz");
  //Bottom left
  c1->cd(7);
  nuc_inv_mass_y_1pXn1pim_data->Draw("colz");
  //Bottom right
  c1->cd(8);
  nuc_inv_mass_y_1pXn1pim_MC->Draw("colz");
  //Very bottom left
  c1->cd(9);
  nuc_inv_mass_omega_1pXn1pim_data->Draw("colz");
  //Very bottom right
  c1->cd(10);
  nuc_inv_mass_omega_1pXn1pim_MC->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //Twenty-seventh page (0p1pip W vs. remaining kinematic variables)
  //Divide page into 2x3 grid
  c1->Divide(2,3);
  //Top left
  c1->cd(1);
  W_q2_1pXn1pim_data->Draw("colz");
  //Top right
  c1->cd(2);
  W_q2_1pXn1pim_MC->Draw("colz");
  //Middle left
  c1->cd(3);
  W_xb_1pXn1pim_data->Draw("colz");
  //Middle right
  c1->cd(4);
  W_xb_1pXn1pim_MC->Draw("colz");
  //Bottom left
  c1->cd(5);
  W_y_1pXn1pim_data->Draw("colz");
  //Bottom right
  c1->cd(6);
  W_y_1pXn1pim_MC->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();

  //Twenty-eighth page (0p1pip totrecoe vs. kinematic variables)
  //Divide page into 2x2 grid
  c1->Divide(2,2);
  //Top left
  c1->cd(1);
  q2_xb_1pXn1pim_data->Draw("colz");
  //Top right
  c1->cd(2);
  q2_xb_1pXn1pim_MC->Draw("colz");
  //Bottom left
  c1->cd(3);
  q2_y_1pXn1pim_data->Draw("colz");
  //Bottom right
  c1->cd(4);
  q2_y_1pXn1pim_MC->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();
  
  //Twenty-ninth page (0p1pip xb vs y)
  //Divide page into 2x2 grid
  c1->Divide(2,2);
  //Left
  c1->cd(1);
  y_xb_1pXn1pim_data->Draw("colz");
  //Right
  c1->cd(2);
  y_xb_1pXn1pim_MC->Draw("colz");
  //Print page to file
  c1->Print(fileName);
  c1->Clear();
  
  //Primakoff plot
  c1->Divide(1,1);
  c1->cd(1);
  invmass_vs_theta_0pXn1pip1pim_data->Draw("colz");


  //Print final page to file and close pdf
  c1->Print(fileName+")");
  c1->Clear();

  //Write out histograms to file
  q0_vs_q3_MC->Write();
  q0_vs_q3_data->Write();
  q0_vs_q3_1pXn0pi_MC->Write();
  q0_vs_q3_1pXn0pi_data->Write();
  q0_vs_q3_2pXn0pi_MC->Write();
  q0_vs_q3_2pXn0pi_data->Write();
  mult_p_stacked->Write();
  mult_n_stacked->Write();
  mult_p_data->Write();
  mult_p_MC->Write();
  mult_n_data->Write();
  mult_n_MC->Write();
  mult_pip_stacked->Write();
  mult_pip_data->Write();
  mult_pip_MC->Write();
  mult_pim_stacked->Write();
  mult_pim_data->Write();
  mult_pim_MC->Write();
  mult_p_pis_data->Write();
  mult_p_pis_MC->Write();
  q2_theta_data->Write();
  q2_theta_MC->Write();
  q2_stacked->Write();
  q2_data->Write();
  q2_MC->Write();
  q2_W_data->Write();
  q2_W_MC->Write();
  theta_mom_data->Write();
  theta_mom_MC->Write();
  p_chi2_stacked->Write();
  n_chi2_stacked->Write();
  pip_chi2_stacked->Write();
  pim_chi2_stacked->Write();
  el_vz_h_stacked->Write();
  el_vz_h_data->Write();
  el_vz_h_MC->Write();
  W_stacked->Write();
  W_data->Write();
  W_MC->Write();
  energy_transfer_0_stacked->Write();
  energy_transfer_0_data->Write();
  energy_transfer_0_MC->Write();
  energy_transfer_1_stacked->Write();
  energy_transfer_1_data->Write();
  energy_transfer_1_MC->Write();
  energy_transfer_2_stacked->Write();
  energy_transfer_2_data->Write();
  energy_transfer_2_MC->Write();
  energy_transfer_3_stacked->Write();
  energy_transfer_3_data->Write();
  energy_transfer_3_MC->Write();
  invariant_mass_stacked->Write();
  invariant_mass_data->Write();
  invariant_mass_MC->Write();
  ecal_data->Write();
  ecal_MC->Write();
  p_proton_stacked->Write();
  p_proton_data->Write();
  p_proton_MC->Write();
  p_miss_stacked->Write();
  p_miss_theta_stacked->Write();
  xb_stacked->Write();
  y_stacked->Write();

  //Specific topologies write out
  W_qe_1p_cuts_stacked->Write();
  xb_qe_1p_cuts_stacked->Write();
  p_proton_qe_1p_cuts_stacked->Write();
  W_mec_1p_cuts_stacked->Write();
  xb_mec_1p_cuts_stacked->Write();
  p_proton_mec_1p_cuts_stacked->Write();
  W_mec_2p_cuts_stacked->Write();
  xb_mec_2p_cuts_stacked->Write();
  p_proton_mec_2p_cuts_stacked->Write();
  W_res_0p1pi_cuts_stacked->Write();
  xb_res_0p1pi_cuts_stacked->Write();
  p_pi_res_0p1pi_cuts_stacked->Write();
  W_res_1p1pi_cuts_stacked->Write();
  xb_res_1p1pi_cuts_stacked->Write();
  p_pi_res_1p1pi_cuts_stacked->Write();
  p_proton_res_1p1pi_cuts_stacked->Write();
  //For Larry plot
  q2_0p1n0pi_MC->Write(); q2_0p1n0pi_data->Write();
  q2_1pXn0pi_MC->Write(); q2_1pXn0pi_data->Write();
  totrecoe_0p1n0pi_MC->Write(); totrecoe_0p1n0pi_data->Write();
  totrecoe_1pXn0pi_MC->Write(); totrecoe_1pXn0pi_data->Write();
  W_0p1n0pi_MC->Write(); W_0p1n0pi_data->Write();
  W_1pXn0pi_MC->Write(); W_1pXn0pi_data->Write();
  
  gBenchmark->Stop("timer");
  gBenchmark->Print("timer");
  
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count() << std:: endl;
  std::cout << "Data events = " << counter_data << ", and MC events = " << counter_MC << std::endl;
}
