#include <iostream>
#include <vector>
#include <algorithm>

#include "TH1F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TH2F.h"

#include "ntpleUtils.h"
#include "ConfigParser.h"

int main (int argc, char** argv) {

  if(argc!=2){ std::cout<<" Not correct number of input parameter --> Need Just one cfg file exit "<<std::endl; return -1; }

  // Input Information 
  parseConfigFile(argv[1]);

  // Load TTree Lybrary
  gSystem->Load("libTree.so");


  std::string InputFileName  = gConfigParser -> readStringOption("Input::InputRootFile");
  std::string InputTreeName  = gConfigParser -> readStringOption("Input::InputTreeName");

  std::string OutputFileName = gConfigParser -> readStringOption("Output::OutputRootFileName");
  std::string OutputPlotDir  = gConfigParser -> readStringOption("Output::OutputPlotDirectory");;
  std::string OutputRootDir  = gConfigParser -> readStringOption("Output::OutputRootDirectory");

  int  NJet         = gConfigParser -> readIntOption("Options::NumberOfJets");
  int  Nbin         = gConfigParser -> readIntOption("Options::Nbin");
  int  Nbin2D       = gConfigParser -> readIntOption("Options::Nbin2D");

  double EtaMinMax  = gConfigParser -> readDoubleOption("Options::EtaMinMax");
  double PhiMinMax  = gConfigParser -> readDoubleOption("Options::PhiMinMax");
  double PtEMassMin = gConfigParser -> readDoubleOption("Options::PtEMassMin");
  double PtEMassMax = gConfigParser -> readDoubleOption("Options::PtEMassMax");
  double DRMin      = gConfigParser -> readDoubleOption("Options::DRMin");
  double DRMax      = gConfigParser -> readDoubleOption("Options::DRMax");

  double PtEMassMin2D = gConfigParser -> readDoubleOption("Options::PtEMassMin2D");
  double PtEMassMax2D = gConfigParser -> readDoubleOption("Options::PtEMassMax2D");
  double DRMin2D      = gConfigParser -> readDoubleOption("Options::DRMin2D");
  double DRMax2D      = gConfigParser -> readDoubleOption("Options::DRMax2D");

  double DRMatchingCut      = gConfigParser -> readDoubleOption("Options::DRMatchingCut");
  double PtMatchingFraction = gConfigParser -> readDoubleOption("Options::PtMatchingFraction");
  double PtCut              = gConfigParser -> readDoubleOption("Options::PtCut");

  // Set Root style from global enviroment path
  std::string ROOTStyle =  getenv ("ROOTStyle");

  gROOT->ProcessLine((".x "+ROOTStyle+"/rootLogon.C").c_str());
  gROOT->ProcessLine((".x "+ROOTStyle+"/rootPalette.C").c_str());
  gROOT->ProcessLine((".x "+ROOTStyle+"/rootColors.C").c_str());
  gROOT->ProcessLine((".x "+ROOTStyle+"/setTDRStyle.C").c_str());
  
  // Set more Style Optionns
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);

  // Cancel output png in the directory 

  std::string command =  " if [ -e "+OutputPlotDir+"/ ] ; then rm "+OutputPlotDir+"/*.png ;  else  mkdir "+OutputPlotDir+" ; fi " ;
  system(command.c_str()); 

  command = " if [ ! -e "+OutputRootDir+"/ ] ; then mkdir "+OutputRootDir+" ; fi " ;
  system(command.c_str()); 


  // Open input File and Tree
  TFile * inputFile  = new TFile(InputFileName.c_str(),"READ");
  TFile * outputFile = new TFile((OutputRootDir+"/"+OutputFileName).c_str(),"RECREATE");
  TTree * fTree      = (TTree*) inputFile -> Get(InputTreeName.c_str());
  
  if(fTree==0) return -1;

  // Declaration of Branch Variables for SetBranchAddress

  float W_Parton_px [2],W_Parton_py[2],W_Parton_pz[2],W_Parton_e[2] ;

  float GroomedJet_CA8_pt[6],GroomedJet_CA8_eta[6],GroomedJet_CA8_phi[6],GroomedJet_CA8_e[6];
 
  float GroomedJet_CA8_pt_pr[6],GroomedJet_CA8_eta_pr[6],GroomedJet_CA8_phi_pr[6],GroomedJet_CA8_e_pr[6];
  float GroomedJet_CA8_pt_tr[6],GroomedJet_CA8_eta_tr[6],GroomedJet_CA8_phi_tr[6],GroomedJet_CA8_e_tr[6];
  float GroomedJet_CA8_pt_ft[6],GroomedJet_CA8_eta_ft[6],GroomedJet_CA8_phi_ft[6],GroomedJet_CA8_e_ft[6];

  float GroomedJet_CA8_prsubjet1_px[6],GroomedJet_CA8_prsubjet1_py[6],GroomedJet_CA8_prsubjet1_pz[6],GroomedJet_CA8_prsubjet1_e[6];
  float GroomedJet_CA8_prsubjet2_px[6],GroomedJet_CA8_prsubjet2_py[6],GroomedJet_CA8_prsubjet2_pz[6],GroomedJet_CA8_prsubjet2_e[6];

  float GroomedJet_AK5_pt[6],GroomedJet_AK5_eta[6],GroomedJet_AK5_phi[6],GroomedJet_AK5_e[6];

  float GenGroomedJet_CA8_pt[6],GenGroomedJet_CA8_eta[6],GenGroomedJet_CA8_phi[6],GenGroomedJet_CA8_e[6];
  float GenGroomedJet_AK5_pt[6],GenGroomedJet_AK5_eta[6],GenGroomedJet_AK5_phi[6],GenGroomedJet_AK5_e[6];

  float GenGroomedJet_CA8_pt_pr[6],GenGroomedJet_CA8_eta_pr[6],GenGroomedJet_CA8_phi_pr[6],GenGroomedJet_CA8_e_pr[6];
  float GenGroomedJet_CA8_pt_tr[6],GenGroomedJet_CA8_eta_tr[6],GenGroomedJet_CA8_phi_tr[6],GenGroomedJet_CA8_e_tr[6];
  float GenGroomedJet_CA8_pt_ft[6],GenGroomedJet_CA8_eta_ft[6],GenGroomedJet_CA8_phi_ft[6],GenGroomedJet_CA8_e_ft[6];

  float GenGroomedJet_CA8_prsubjet1_px[6],GenGroomedJet_CA8_prsubjet1_py[6],GenGroomedJet_CA8_prsubjet1_pz[6],GenGroomedJet_CA8_prsubjet1_e[6];
  float GenGroomedJet_CA8_prsubjet2_px[6],GenGroomedJet_CA8_prsubjet2_py[6],GenGroomedJet_CA8_prsubjet2_pz[6],GenGroomedJet_CA8_prsubjet2_e[6];

  float W_Lepton_px,W_Lepton_py,W_Lepton_pz,W_Lepton_e ;
  float W_Met_px,W_Met_py,W_Met_pz,W_Met_e ;
  float W_px,W_py,W_pz,W_e;

  int   ggdboostedWevt ;

  float W_pt, event_met_pfmet ;
 

  // Set Branch Address 
  fTree->SetBranchStatus("*",0);
  fTree->SetBranchStatus("W_Parton*",1);
  fTree->SetBranchStatus("GroomedJet_CA8*",1);
  fTree->SetBranchStatus("GroomedJet_AK5*",1);
  fTree->SetBranchStatus("GenGroomedJet_CA8*",1);
  fTree->SetBranchStatus("GenGroomedJet_AK5*",1);
  fTree->SetBranchStatus("W_Lepton*",1);
  fTree->SetBranchStatus("W_Met*",1);

  fTree->SetBranchStatus("ggdboostedWevt",1);
  fTree->SetBranchStatus("event_met_pfmet",1);

  fTree->SetBranchStatus("W_Lepton*",1);
  fTree->SetBranchStatus("W_p*",1);
  fTree->SetBranchStatus("W_e",1);

  fTree->SetBranchAddress("W_Parton_px[2]",&W_Parton_px);
  fTree->SetBranchAddress("W_Parton_py[2]",&W_Parton_py);
  fTree->SetBranchAddress("W_Parton_pz[2]",&W_Parton_pz);
  fTree->SetBranchAddress("W_Parton_E[2]",&W_Parton_e);

  fTree->SetBranchAddress("W_Lepton_px",&W_Lepton_px);
  fTree->SetBranchAddress("W_Lepton_py",&W_Lepton_py);
  fTree->SetBranchAddress("W_Lepton_pz",&W_Lepton_pz);
  fTree->SetBranchAddress("W_Lepton_E",&W_Lepton_e);

  fTree->SetBranchAddress("W_px",&W_px);
  fTree->SetBranchAddress("W_py",&W_py);
  fTree->SetBranchAddress("W_pz",&W_pz);
  fTree->SetBranchAddress("W_e",&W_e);
  fTree->SetBranchAddress("W_pt",&W_pt);

  fTree->SetBranchAddress("W_Met_px",&W_Met_px);
  fTree->SetBranchAddress("W_Met_py",&W_Met_py);
  fTree->SetBranchAddress("W_Met_pz",&W_Met_pz);
  fTree->SetBranchAddress("W_Met_E",&W_Met_e);

  fTree->SetBranchAddress("ggdboostedWevt",&ggdboostedWevt);
  fTree->SetBranchAddress("event_met_pfmet",&event_met_pfmet);
  
  fTree->SetBranchAddress("GroomedJet_CA8_pt",&GroomedJet_CA8_pt);
  fTree->SetBranchAddress("GroomedJet_CA8_eta",&GroomedJet_CA8_eta);
  fTree->SetBranchAddress("GroomedJet_CA8_phi",&GroomedJet_CA8_phi);
  fTree->SetBranchAddress("GroomedJet_CA8_e",&GroomedJet_CA8_e);

  fTree->SetBranchAddress("GroomedJet_CA8_pt_pr",&GroomedJet_CA8_pt_pr);
  fTree->SetBranchAddress("GroomedJet_CA8_eta_pr",&GroomedJet_CA8_eta_pr);
  fTree->SetBranchAddress("GroomedJet_CA8_phi_pr",&GroomedJet_CA8_phi_pr);
  fTree->SetBranchAddress("GroomedJet_CA8_e_pr",&GroomedJet_CA8_e_pr);

  fTree->SetBranchAddress("GroomedJet_CA8_pt_tr",&GroomedJet_CA8_pt_tr);
  fTree->SetBranchAddress("GroomedJet_CA8_eta_tr",&GroomedJet_CA8_eta_tr);
  fTree->SetBranchAddress("GroomedJet_CA8_phi_tr",&GroomedJet_CA8_phi_tr);
  fTree->SetBranchAddress("GroomedJet_CA8_e_tr",&GroomedJet_CA8_e_tr);

  fTree->SetBranchAddress("GroomedJet_CA8_pt_ft",&GroomedJet_CA8_pt_ft);
  fTree->SetBranchAddress("GroomedJet_CA8_eta_ft",&GroomedJet_CA8_eta_ft);
  fTree->SetBranchAddress("GroomedJet_CA8_phi_ft",&GroomedJet_CA8_phi_ft);
  fTree->SetBranchAddress("GroomedJet_CA8_e_ft",&GroomedJet_CA8_e_ft);

  fTree->SetBranchAddress("GroomedJet_CA8_prsubjet1_px",&GroomedJet_CA8_prsubjet1_px);
  fTree->SetBranchAddress("GroomedJet_CA8_prsubjet1_py",&GroomedJet_CA8_prsubjet1_py);
  fTree->SetBranchAddress("GroomedJet_CA8_prsubjet1_pz",&GroomedJet_CA8_prsubjet1_pz);
  fTree->SetBranchAddress("GroomedJet_CA8_prsubjet1_e",&GroomedJet_CA8_prsubjet1_e);

  fTree->SetBranchAddress("GroomedJet_CA8_prsubjet2_px",&GroomedJet_CA8_prsubjet2_px);
  fTree->SetBranchAddress("GroomedJet_CA8_prsubjet2_py",&GroomedJet_CA8_prsubjet2_py);
  fTree->SetBranchAddress("GroomedJet_CA8_prsubjet2_pz",&GroomedJet_CA8_prsubjet2_pz);
  fTree->SetBranchAddress("GroomedJet_CA8_prsubjet2_e",&GroomedJet_CA8_prsubjet2_e);

  fTree->SetBranchAddress("GroomedJet_AK5_pt",&GroomedJet_AK5_pt);
  fTree->SetBranchAddress("GroomedJet_AK5_eta",&GroomedJet_AK5_eta);
  fTree->SetBranchAddress("GroomedJet_AK5_phi",&GroomedJet_AK5_phi);
  fTree->SetBranchAddress("GroomedJet_AK5_e",&GroomedJet_AK5_e);

  fTree->SetBranchAddress("GenGroomedJet_CA8_pt",&GenGroomedJet_CA8_pt);
  fTree->SetBranchAddress("GenGroomedJet_CA8_eta",&GenGroomedJet_CA8_eta);
  fTree->SetBranchAddress("GenGroomedJet_CA8_phi",&GenGroomedJet_CA8_phi);
  fTree->SetBranchAddress("GenGroomedJet_CA8_e",&GenGroomedJet_CA8_e);

  fTree->SetBranchAddress("GenGroomedJet_CA8_prsubjet1_px",&GenGroomedJet_CA8_prsubjet1_px);
  fTree->SetBranchAddress("GenGroomedJet_CA8_prsubjet1_py",&GenGroomedJet_CA8_prsubjet1_py);
  fTree->SetBranchAddress("GenGroomedJet_CA8_prsubjet1_pz",&GenGroomedJet_CA8_prsubjet1_pz);
  fTree->SetBranchAddress("GenGroomedJet_CA8_prsubjet1_e",&GenGroomedJet_CA8_prsubjet1_e);

  fTree->SetBranchAddress("GenGroomedJet_CA8_prsubjet2_px",&GenGroomedJet_CA8_prsubjet2_px);
  fTree->SetBranchAddress("GenGroomedJet_CA8_prsubjet2_py",&GenGroomedJet_CA8_prsubjet2_py);
  fTree->SetBranchAddress("GenGroomedJet_CA8_prsubjet2_pz",&GenGroomedJet_CA8_prsubjet2_pz);
  fTree->SetBranchAddress("GenGroomedJet_CA8_prsubjet2_e",&GenGroomedJet_CA8_prsubjet2_e);

  fTree->SetBranchAddress("GenGroomedJet_AK5_pt",&GenGroomedJet_AK5_pt);
  fTree->SetBranchAddress("GenGroomedJet_AK5_eta",&GenGroomedJet_AK5_eta);
  fTree->SetBranchAddress("GenGroomedJet_AK5_phi",&GenGroomedJet_AK5_phi);
  fTree->SetBranchAddress("GenGroomedJet_AK5_e",&GenGroomedJet_AK5_e);

  fTree->SetBranchAddress("GenGroomedJet_CA8_pt_pr",&GenGroomedJet_CA8_pt_pr);
  fTree->SetBranchAddress("GenGroomedJet_CA8_eta_pr",&GenGroomedJet_CA8_eta_pr);
  fTree->SetBranchAddress("GenGroomedJet_CA8_phi_pr",&GenGroomedJet_CA8_phi_pr);
  fTree->SetBranchAddress("GenGroomedJet_CA8_e_pr",&GenGroomedJet_CA8_e_pr);

  fTree->SetBranchAddress("GenGroomedJet_CA8_pt_tr",&GenGroomedJet_CA8_pt_tr);
  fTree->SetBranchAddress("GenGroomedJet_CA8_eta_tr",&GenGroomedJet_CA8_eta_tr);
  fTree->SetBranchAddress("GenGroomedJet_CA8_phi_tr",&GenGroomedJet_CA8_phi_tr);
  fTree->SetBranchAddress("GenGroomedJet_CA8_e_tr",&GenGroomedJet_CA8_e_tr);

  fTree->SetBranchAddress("GenGroomedJet_CA8_pt_ft",&GenGroomedJet_CA8_pt_ft);
  fTree->SetBranchAddress("GenGroomedJet_CA8_eta_ft",&GenGroomedJet_CA8_eta_ft);
  fTree->SetBranchAddress("GenGroomedJet_CA8_phi_ft",&GenGroomedJet_CA8_phi_ft);
  fTree->SetBranchAddress("GenGroomedJet_CA8_e_ft",&GenGroomedJet_CA8_e_ft);

  // Gen Jet Histogram 

  TH1F ** DPt_GenGroomedJet_AK5_GenWhad  = new TH1F* [NJet];
  TH1F ** DEta_GenGroomedJet_AK5_GenWhad = new TH1F* [NJet];
  TH1F ** DPhi_GenGroomedJet_AK5_GenWhad = new TH1F* [NJet];
  TH1F ** DR_GenGroomedJet_AK5_GenWhad   = new TH1F* [NJet];
  TH1F ** DE_GenGroomedJet_AK5_GenWhad   = new TH1F* [NJet];
  TH1F ** DMass_GenGroomedJet_AK5_GenWhad   = new TH1F* [NJet];

  TH2F ** DPt_DR_GenGroomedJet_AK5_GenWhad = new TH2F* [NJet];
  TH2F ** DMass_DR_GenGroomedJet_AK5_GenWhad = new TH2F* [NJet];

  TH1F ** DPt_GenGroomedJet_CA8_GenWhad  = new TH1F* [NJet];
  TH1F ** DEta_GenGroomedJet_CA8_GenWhad = new TH1F* [NJet];
  TH1F ** DPhi_GenGroomedJet_CA8_GenWhad = new TH1F* [NJet];
  TH1F ** DR_GenGroomedJet_CA8_GenWhad   = new TH1F* [NJet];
  TH1F ** DE_GenGroomedJet_CA8_GenWhad   = new TH1F* [NJet];
  TH1F ** DMass_GenGroomedJet_CA8_GenWhad   = new TH1F* [NJet];

  TH2F ** DPt_DR_GenGroomedJet_CA8_GenWhad = new TH2F* [NJet];
  TH2F ** DMass_DR_GenGroomedJet_CA8_GenWhad = new TH2F* [NJet];

  TH1F ** DPt_GenGroomedJet_CA8_subjet1_WParton  = new TH1F* [NJet];
  TH1F ** DEta_GenGroomedJet_CA8_subjet1_WParton = new TH1F* [NJet];
  TH1F ** DPhi_GenGroomedJet_CA8_subjet1_WParton = new TH1F* [NJet];
  TH1F ** DR_GenGroomedJet_CA8_subjet1_WParton   = new TH1F* [NJet];
  TH1F ** DE_GenGroomedJet_CA8_subjet1_WParton   = new TH1F* [NJet];
  TH1F ** DMass_GenGroomedJet_CA8_subjet1_WParton   = new TH1F* [NJet];

  TH2F ** DPt_DR_GenGroomedJet_CA8_subjet1_WParton = new TH2F* [NJet];
  TH2F ** DMass_DR_GenGroomedJet_CA8_subjet1_WParton = new TH2F* [NJet];

  TH1F ** DPt_GenGroomedJet_CA8_subjet2_WParton  = new TH1F* [NJet];
  TH1F ** DEta_GenGroomedJet_CA8_subjet2_WParton = new TH1F* [NJet];
  TH1F ** DPhi_GenGroomedJet_CA8_subjet2_WParton = new TH1F* [NJet];
  TH1F ** DR_GenGroomedJet_CA8_subjet2_WParton   = new TH1F* [NJet];
  TH1F ** DE_GenGroomedJet_CA8_subjet2_WParton   = new TH1F* [NJet];
  TH1F ** DMass_GenGroomedJet_CA8_subjet2_WParton   = new TH1F* [NJet];

  TH2F ** DPt_DR_GenGroomedJet_CA8_subjet2_WParton = new TH2F* [NJet];
  TH2F ** DMass_DR_GenGroomedJet_CA8_subjet2_WParton = new TH2F* [NJet];

  TH1F ** DPt_GenGroomedJet_CA8_GenWhad_pr  = new TH1F* [NJet];
  TH1F ** DEta_GenGroomedJet_CA8_GenWhad_pr = new TH1F* [NJet];
  TH1F ** DPhi_GenGroomedJet_CA8_GenWhad_pr = new TH1F* [NJet];
  TH1F ** DR_GenGroomedJet_CA8_GenWhad_pr   = new TH1F* [NJet];
  TH1F ** DE_GenGroomedJet_CA8_GenWhad_pr   = new TH1F* [NJet];
  TH1F ** DMass_GenGroomedJet_CA8_GenWhad_pr   = new TH1F* [NJet];

  TH2F ** DPt_DR_GenGroomedJet_CA8_GenWhad_pr = new TH2F* [NJet];
  TH2F ** DMass_DR_GenGroomedJet_CA8_GenWhad_pr = new TH2F* [NJet];

  TH1F ** DPt_GenGroomedJet_CA8_GenWhad_tr  = new TH1F* [NJet];
  TH1F ** DEta_GenGroomedJet_CA8_GenWhad_tr = new TH1F* [NJet];
  TH1F ** DPhi_GenGroomedJet_CA8_GenWhad_tr = new TH1F* [NJet];
  TH1F ** DR_GenGroomedJet_CA8_GenWhad_tr   = new TH1F* [NJet];
  TH1F ** DE_GenGroomedJet_CA8_GenWhad_tr   = new TH1F* [NJet];
  TH1F ** DMass_GenGroomedJet_CA8_GenWhad_tr   = new TH1F* [NJet];

  TH2F ** DPt_DR_GenGroomedJet_CA8_GenWhad_tr = new TH2F* [NJet];
  TH2F ** DMass_DR_GenGroomedJet_CA8_GenWhad_tr = new TH2F* [NJet];

  TH1F ** DPt_GenGroomedJet_CA8_GenWhad_ft  = new TH1F* [NJet];
  TH1F ** DEta_GenGroomedJet_CA8_GenWhad_ft = new TH1F* [NJet];
  TH1F ** DPhi_GenGroomedJet_CA8_GenWhad_ft = new TH1F* [NJet];
  TH1F ** DR_GenGroomedJet_CA8_GenWhad_ft   = new TH1F* [NJet];
  TH1F ** DE_GenGroomedJet_CA8_GenWhad_ft   = new TH1F* [NJet];
  TH1F ** DMass_GenGroomedJet_CA8_GenWhad_ft   = new TH1F* [NJet];

  TH2F ** DPt_DR_GenGroomedJet_CA8_GenWhad_ft = new TH2F* [NJet];
  TH2F ** DMass_DR_GenGroomedJet_CA8_GenWhad_ft = new TH2F* [NJet];
 
  // reco jet Histogramm

  TH1F ** DPt_GroomedJet_AK5_GenWhad  = new TH1F* [NJet];
  TH1F ** DEta_GroomedJet_AK5_GenWhad = new TH1F* [NJet];
  TH1F ** DPhi_GroomedJet_AK5_GenWhad = new TH1F* [NJet];
  TH1F ** DR_GroomedJet_AK5_GenWhad   = new TH1F* [NJet];
  TH1F ** DE_GroomedJet_AK5_GenWhad   = new TH1F* [NJet];
  TH1F ** DMass_GroomedJet_AK5_GenWhad   = new TH1F* [NJet];

  TH2F ** DPt_DR_GroomedJet_AK5_GenWhad = new TH2F* [NJet];
  TH2F ** DMass_DR_GroomedJet_AK5_GenWhad = new TH2F* [NJet];

  TH1F ** DPt_GroomedJet_CA8_GenWhad  = new TH1F* [NJet];
  TH1F ** DEta_GroomedJet_CA8_GenWhad = new TH1F* [NJet];
  TH1F ** DPhi_GroomedJet_CA8_GenWhad = new TH1F* [NJet];
  TH1F ** DR_GroomedJet_CA8_GenWhad   = new TH1F* [NJet];
  TH1F ** DE_GroomedJet_CA8_GenWhad   = new TH1F* [NJet];
  TH1F ** DMass_GroomedJet_CA8_GenWhad   = new TH1F* [NJet];

  TH2F ** DPt_DR_GroomedJet_CA8_GenWhad = new TH2F* [NJet];
  TH2F ** DMass_DR_GroomedJet_CA8_GenWhad = new TH2F* [NJet];

  TH1F ** DPt_GroomedJet_CA8_subjet1_WParton  = new TH1F* [NJet];
  TH1F ** DEta_GroomedJet_CA8_subjet1_WParton = new TH1F* [NJet];
  TH1F ** DPhi_GroomedJet_CA8_subjet1_WParton = new TH1F* [NJet];
  TH1F ** DR_GroomedJet_CA8_subjet1_WParton   = new TH1F* [NJet];
  TH1F ** DE_GroomedJet_CA8_subjet1_WParton   = new TH1F* [NJet];
  TH1F ** DMass_GroomedJet_CA8_subjet1_WParton   = new TH1F* [NJet];

  TH2F ** DPt_DR_GroomedJet_CA8_subjet1_WParton = new TH2F* [NJet];
  TH2F ** DMass_DR_GroomedJet_CA8_subjet1_WParton = new TH2F* [NJet];

  TH1F ** DPt_GroomedJet_CA8_subjet2_WParton  = new TH1F* [NJet];
  TH1F ** DEta_GroomedJet_CA8_subjet2_WParton = new TH1F* [NJet];
  TH1F ** DPhi_GroomedJet_CA8_subjet2_WParton = new TH1F* [NJet];
  TH1F ** DR_GroomedJet_CA8_subjet2_WParton   = new TH1F* [NJet];
  TH1F ** DE_GroomedJet_CA8_subjet2_WParton   = new TH1F* [NJet];
  TH1F ** DMass_GroomedJet_CA8_subjet2_WParton   = new TH1F* [NJet];

  TH2F ** DPt_DR_GroomedJet_CA8_subjet2_WParton = new TH2F* [NJet];
  TH2F ** DMass_DR_GroomedJet_CA8_subjet2_WParton = new TH2F* [NJet];

  TH1F ** DPt_GroomedJet_CA8_GenWhad_pr  = new TH1F* [NJet];
  TH1F ** DEta_GroomedJet_CA8_GenWhad_pr = new TH1F* [NJet];
  TH1F ** DPhi_GroomedJet_CA8_GenWhad_pr = new TH1F* [NJet];
  TH1F ** DR_GroomedJet_CA8_GenWhad_pr   = new TH1F* [NJet];
  TH1F ** DE_GroomedJet_CA8_GenWhad_pr   = new TH1F* [NJet];
  TH1F ** DMass_GroomedJet_CA8_GenWhad_pr   = new TH1F* [NJet];

  TH2F ** DPt_DR_GroomedJet_CA8_GenWhad_pr = new TH2F* [NJet];
  TH2F ** DMass_DR_GroomedJet_CA8_GenWhad_pr = new TH2F* [NJet];

  TH1F ** DPt_GroomedJet_CA8_GenWhad_tr  = new TH1F* [NJet];
  TH1F ** DEta_GroomedJet_CA8_GenWhad_tr = new TH1F* [NJet];
  TH1F ** DPhi_GroomedJet_CA8_GenWhad_tr = new TH1F* [NJet];
  TH1F ** DR_GroomedJet_CA8_GenWhad_tr   = new TH1F* [NJet];
  TH1F ** DE_GroomedJet_CA8_GenWhad_tr   = new TH1F* [NJet];
  TH1F ** DMass_GroomedJet_CA8_GenWhad_tr   = new TH1F* [NJet];

  TH2F ** DPt_DR_GroomedJet_CA8_GenWhad_tr = new TH2F* [NJet];
  TH2F ** DMass_DR_GroomedJet_CA8_GenWhad_tr = new TH2F* [NJet];

  TH1F ** DPt_GroomedJet_CA8_GenWhad_ft  = new TH1F* [NJet];
  TH1F ** DEta_GroomedJet_CA8_GenWhad_ft = new TH1F* [NJet];
  TH1F ** DPhi_GroomedJet_CA8_GenWhad_ft = new TH1F* [NJet];
  TH1F ** DR_GroomedJet_CA8_GenWhad_ft   = new TH1F* [NJet];
  TH1F ** DE_GroomedJet_CA8_GenWhad_ft   = new TH1F* [NJet];
  TH1F ** DMass_GroomedJet_CA8_GenWhad_ft   = new TH1F* [NJet];

  TH2F ** DPt_DR_GroomedJet_CA8_GenWhad_ft = new TH2F* [NJet];
  TH2F ** DMass_DR_GroomedJet_CA8_GenWhad_ft = new TH2F* [NJet];

  TH1F ** DPt_GroomedJet_CA8_WLep  = new TH1F* [NJet];
  TH1F ** DEta_GroomedJet_CA8_WLep = new TH1F* [NJet];
  TH1F ** DPhi_GroomedJet_CA8_WLep = new TH1F* [NJet];
  TH1F ** DR_GroomedJet_CA8_WLep   = new TH1F* [NJet];
  TH1F ** DE_GroomedJet_CA8_WLep = new TH1F* [NJet];
  TH1F ** DMass_GroomedJet_CA8_WLep   = new TH1F* [NJet];

  TH2F ** DPt_DR_GroomedJet_CA8_WLep = new TH2F* [NJet];
  TH2F ** DMass_DR_GroomedJet_CA8_WLep = new TH2F* [NJet];
  

  // Allocate Histogramm

  for( int i =0 ; i<NJet ; i++){

   TString name = Form("DPt_GenGroomedJet_AK5_[%i]_GenWhad",i);
   std::string Name = name.Data();
   // Gen Jet Histogramm Allocation

   DPt_GenGroomedJet_AK5_GenWhad[i]  = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DPt_GenGroomedJet_AK5_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DEta_GenGroomedJet_AK5_[%i]_GenWhad",i);
   DEta_GenGroomedJet_AK5_GenWhad[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-EtaMinMax,EtaMinMax);
   DEta_GenGroomedJet_AK5_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DPhi_GenGroomedJet_AK5_[%i]_GenWhad",i);
   DPhi_GenGroomedJet_AK5_GenWhad[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-PhiMinMax,PhiMinMax);
   DPhi_GenGroomedJet_AK5_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DR_GenGroomedJet_AK5_[%i]_GenWhad",i);
   DR_GenGroomedJet_AK5_GenWhad[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,DRMin,DRMax);
   DR_GenGroomedJet_AK5_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DE_GenGroomedJet_AK5_[%i]_GenWhad",i);
   DE_GenGroomedJet_AK5_GenWhad[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DE_GenGroomedJet_AK5_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DMass_GenGroomedJet_AK5_[%i]_GenWhad",i);
   DMass_GenGroomedJet_AK5_GenWhad[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DMass_GenGroomedJet_AK5_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());

   Name = Form("DPt_DR_GenGroomedJet_AK5_[%i]_GenWhad",i);
   DPt_DR_GenGroomedJet_AK5_GenWhad[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DPt_DR_GenGroomedJet_AK5_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DMass_DR_GenGroomedJet_AK5_[%i]_GenWhad",i);
   DMass_DR_GenGroomedJet_AK5_GenWhad[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DMass_DR_GenGroomedJet_AK5_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());


   Name = Form("DPt_GenGroomedJet_CA8_[%i]_GenWhad",i);
   DPt_GenGroomedJet_CA8_GenWhad[i]  = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DPt_GenGroomedJet_CA8_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DEta_GenGroomedJet_CA8_[%i]_GenWhad",i);
   DEta_GenGroomedJet_CA8_GenWhad[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-EtaMinMax,EtaMinMax);
   DEta_GenGroomedJet_CA8_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DPhi_GenGroomedJet_CA8_[%i]_GenWhad",i);
   DPhi_GenGroomedJet_CA8_GenWhad[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-PhiMinMax,PhiMinMax);
   DPhi_GenGroomedJet_CA8_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DR_GenGroomedJet_CA8_[%i]_GenWhad",i);
   DR_GenGroomedJet_CA8_GenWhad[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,DRMin,DRMax);
   DR_GenGroomedJet_CA8_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DE_GenGroomedJet_CA8_[%i]_GenWhad",i);
   DE_GenGroomedJet_CA8_GenWhad[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DE_GenGroomedJet_CA8_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DMass_GenGroomedJet_CA8_[%i]_GenWhad",i);
   DMass_GenGroomedJet_CA8_GenWhad[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DMass_GenGroomedJet_CA8_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());

   Name = Form("DPt_DR_GenGroomedJet_CA8_[%i]_GenWhad",i);
   DPt_DR_GenGroomedJet_CA8_GenWhad[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DPt_DR_GenGroomedJet_CA8_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DMass_DR_GenGroomedJet_CA8_[%i]_GenWhad",i);
   DMass_DR_GenGroomedJet_CA8_GenWhad[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DMass_DR_GenGroomedJet_CA8_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());

   Name = Form("DPt_GenGroomedJet_CA8_[%i]_subjet1_WParton",i);
   DPt_GenGroomedJet_CA8_subjet1_WParton[i]  = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DPt_GenGroomedJet_CA8_subjet1_WParton[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DEta_GenGroomedJet_CA8_[%i]_subjet1_WParton",i);
   DEta_GenGroomedJet_CA8_subjet1_WParton[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-EtaMinMax,EtaMinMax);
   DEta_GenGroomedJet_CA8_subjet1_WParton[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DPhi_GenGroomedJet_CA8_[%i]_subjet1_WParton",i);
   DPhi_GenGroomedJet_CA8_subjet1_WParton[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-PhiMinMax,PhiMinMax);
   DPhi_GenGroomedJet_CA8_subjet1_WParton[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DR_GenGroomedJet_CA8_[%i]_subjet1_WParton",i);
   DR_GenGroomedJet_CA8_subjet1_WParton[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,DRMin,DRMax);
   DR_GenGroomedJet_CA8_subjet1_WParton[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DE_GenGroomedJet_CA8_[%i]_subjet1_WParton",i);
   DE_GenGroomedJet_CA8_subjet1_WParton[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DE_GenGroomedJet_CA8_subjet1_WParton[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DMass_GenGroomedJet_CA8_[%i]_subjet1_WParton",i);
   DMass_GenGroomedJet_CA8_subjet1_WParton[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DMass_GenGroomedJet_CA8_subjet1_WParton[i]->GetXaxis()->SetTitle(Name.c_str());

   Name = Form("DPt_DR_GenGroomedJet_CA8_[%i]_subjet1_WParton",i);
   DPt_DR_GenGroomedJet_CA8_subjet1_WParton[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DPt_DR_GenGroomedJet_CA8_subjet1_WParton[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DMass_DR_GenGroomedJet_CA8_[%i]_subjet1_WParton",i);
   DMass_DR_GenGroomedJet_CA8_subjet1_WParton[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DMass_DR_GenGroomedJet_CA8_subjet1_WParton[i]->GetXaxis()->SetTitle(Name.c_str());

   Name = Form("DPt_GenGroomedJet_CA8_[%i]_subjet2_WParton",i);
   DPt_GenGroomedJet_CA8_subjet2_WParton[i]  = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DPt_GenGroomedJet_CA8_subjet2_WParton[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DEta_GenGroomedJet_CA8_[%i]_subjet2_WParton",i);
   DEta_GenGroomedJet_CA8_subjet2_WParton[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-EtaMinMax,EtaMinMax);
   DEta_GenGroomedJet_CA8_subjet2_WParton[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DPhi_GenGroomedJet_CA8_[%i]_subjet2_WParton",i);
   DPhi_GenGroomedJet_CA8_subjet2_WParton[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-PhiMinMax,PhiMinMax);
   DPhi_GenGroomedJet_CA8_subjet2_WParton[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DR_GenGroomedJet_CA8_[%i]_subjet2_WParton",i);
   DR_GenGroomedJet_CA8_subjet2_WParton[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,DRMin,DRMax);
   DR_GenGroomedJet_CA8_subjet2_WParton[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DE_GenGroomedJet_CA8_[%i]_subjet2_WParton",i);
   DE_GenGroomedJet_CA8_subjet2_WParton[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DE_GenGroomedJet_CA8_subjet2_WParton[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DMass_GenGroomedJet_CA8_[%i]_subjet2_WParton",i);
   DMass_GenGroomedJet_CA8_subjet2_WParton[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DMass_GenGroomedJet_CA8_subjet2_WParton[i]->GetXaxis()->SetTitle(Name.c_str());

   Name = Form("DPt_DR_GenGroomedJet_CA8_[%i]_subjet2_WParton",i);
   DPt_DR_GenGroomedJet_CA8_subjet2_WParton[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DPt_DR_GenGroomedJet_CA8_subjet2_WParton[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DMass_DR_GenGroomedJet_CA8_[%i]_subjet2_WParton",i);
   DMass_DR_GenGroomedJet_CA8_subjet2_WParton[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DMass_DR_GenGroomedJet_CA8_subjet2_WParton[i]->GetXaxis()->SetTitle(Name.c_str());

   Name = Form("DPt_GenGroomedJet_CA8_[%i]_GenWhad_pr",i);
   DPt_GenGroomedJet_CA8_GenWhad_pr[i]  = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DPt_GenGroomedJet_CA8_GenWhad_pr[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DEta_GenGroomedJet_CA8_[%i]_GenWhad_pr",i);
   DEta_GenGroomedJet_CA8_GenWhad_pr[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-EtaMinMax,EtaMinMax);
   DEta_GenGroomedJet_CA8_GenWhad_pr[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DPhi_GenGroomedJet_CA8_[%i]_GenWhad_pr",i);
   DPhi_GenGroomedJet_CA8_GenWhad_pr[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-PhiMinMax,PhiMinMax);
   DPhi_GenGroomedJet_CA8_GenWhad_pr[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DR_GenGroomedJet_CA8_[%i]_GenWhad_pr",i);
   DR_GenGroomedJet_CA8_GenWhad_pr[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,DRMin,DRMax);
   DR_GenGroomedJet_CA8_GenWhad_pr[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DE_GenGroomedJet_CA8_[%i]_GenWhad_pr",i);
   DE_GenGroomedJet_CA8_GenWhad_pr[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DE_GenGroomedJet_CA8_GenWhad_pr[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DMass_GenGroomedJet_CA8_[%i]_GenWhad_pr",i);
   DMass_GenGroomedJet_CA8_GenWhad_pr[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DMass_GenGroomedJet_CA8_GenWhad_pr[i]->GetXaxis()->SetTitle(Name.c_str());

   Name = Form("DPt_DR_GenGroomedJet_CA8_[%i]_GenWhad_pr",i);
   DPt_DR_GenGroomedJet_CA8_GenWhad_pr[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DPt_DR_GenGroomedJet_CA8_GenWhad_pr[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DMass_DR_GenGroomedJet_CA8_[%i]_GenWhad_pr",i);
   DMass_DR_GenGroomedJet_CA8_GenWhad_pr[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DMass_DR_GenGroomedJet_CA8_GenWhad_pr[i]->GetXaxis()->SetTitle(Name.c_str());

   Name = Form("DPt_GenGroomedJet_CA8_[%i]_GenWhad_tr",i);
   DPt_GenGroomedJet_CA8_GenWhad_tr[i]  = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DPt_GenGroomedJet_CA8_GenWhad_tr[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DEta_GenGroomedJet_CA8_[%i]_GenWhad_tr",i);
   DEta_GenGroomedJet_CA8_GenWhad_tr[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-EtaMinMax,EtaMinMax);
   DEta_GenGroomedJet_CA8_GenWhad_tr[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DPhi_GenGroomedJet_CA8_[%i]_GenWhad_tr",i);
   DPhi_GenGroomedJet_CA8_GenWhad_tr[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-PhiMinMax,PhiMinMax);
   DPhi_GenGroomedJet_CA8_GenWhad_tr[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DR_GenGroomedJet_CA8_[%i]_GenWhad_tr",i);
   DR_GenGroomedJet_CA8_GenWhad_tr[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,DRMin,DRMax);
   DR_GenGroomedJet_CA8_GenWhad_tr[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DE_GenGroomedJet_CA8_[%i]_GenWhad_tr",i);
   DE_GenGroomedJet_CA8_GenWhad_tr[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DE_GenGroomedJet_CA8_GenWhad_tr[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DMass_GenGroomedJet_CA8_[%i]_GenWhad_tr",i);
   DMass_GenGroomedJet_CA8_GenWhad_tr[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DMass_GenGroomedJet_CA8_GenWhad_tr[i]->GetXaxis()->SetTitle(Name.c_str());

   Name = Form("DPt_DR_GenGroomedJet_CA8_[%i]_GenWhad_tr",i);
   DPt_DR_GenGroomedJet_CA8_GenWhad_tr[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DPt_DR_GenGroomedJet_CA8_GenWhad_tr[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DMass_DR_GenGroomedJet_CA8_[%i]_GenWhad_tr",i);
   DMass_DR_GenGroomedJet_CA8_GenWhad_tr[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DMass_DR_GenGroomedJet_CA8_GenWhad_tr[i]->GetXaxis()->SetTitle(Name.c_str());

   Name = Form("DPt_GenGroomedJet_CA8_[%i]_GenWhad_ft",i);
   DPt_GenGroomedJet_CA8_GenWhad_ft[i]  = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DPt_GenGroomedJet_CA8_GenWhad_ft[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DEta_GenGroomedJet_CA8_[%i]_GenWhad_ft",i);
   DEta_GenGroomedJet_CA8_GenWhad_ft[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-EtaMinMax,EtaMinMax);
   DEta_GenGroomedJet_CA8_GenWhad_ft[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DPhi_GenGroomedJet_CA8_[%i]_GenWhad_ft",i);
   DPhi_GenGroomedJet_CA8_GenWhad_ft[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-PhiMinMax,PhiMinMax);
   DPhi_GenGroomedJet_CA8_GenWhad_ft[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DR_GenGroomedJet_CA8_[%i]_GenWhad_ft",i);
   DR_GenGroomedJet_CA8_GenWhad_ft[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,DRMin,DRMax);
   DR_GenGroomedJet_CA8_GenWhad_ft[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DE_GenGroomedJet_CA8_[%i]_GenWhad_ft",i);
   DE_GenGroomedJet_CA8_GenWhad_ft[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DE_GenGroomedJet_CA8_GenWhad_ft[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DMass_GenGroomedJet_CA8_[%i]_GenWhad_ft",i);
   DMass_GenGroomedJet_CA8_GenWhad_ft[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DMass_GenGroomedJet_CA8_GenWhad_ft[i]->GetXaxis()->SetTitle(Name.c_str());

   Name = Form("DPt_DR_GenGroomedJet_CA8_[%i]_GenWhad_ft",i);
   DPt_DR_GenGroomedJet_CA8_GenWhad_ft[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DPt_DR_GenGroomedJet_CA8_GenWhad_ft[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DMass_DR_GenGroomedJet_CA8_[%i]_GenWhad_ft",i);
   DMass_DR_GenGroomedJet_CA8_GenWhad_ft[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DMass_DR_GenGroomedJet_CA8_GenWhad_ft[i]->GetXaxis()->SetTitle(Name.c_str());

   // Reco Jet Histogramm Allocation

   Name = Form("DPt_GroomedJet_AK5_[%i]_GenWhad",i);
   DPt_GroomedJet_AK5_GenWhad[i]  = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DPt_GroomedJet_AK5_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DEta_GroomedJet_AK5_[%i]_GenWhad",i);
   DEta_GroomedJet_AK5_GenWhad[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-EtaMinMax,EtaMinMax);
   DEta_GroomedJet_AK5_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DPhi_GroomedJet_AK5_[%i]_GenWhad",i);
   DPhi_GroomedJet_AK5_GenWhad[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-PhiMinMax,PhiMinMax);
   DPhi_GroomedJet_AK5_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DR_GroomedJet_AK5_[%i]_GenWhad",i);
   DR_GroomedJet_AK5_GenWhad[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,DRMin,DRMax);
   DR_GroomedJet_AK5_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DE_GroomedJet_AK5_[%i]_GenWhad",i);
   DE_GroomedJet_AK5_GenWhad[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DE_GroomedJet_AK5_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DMass_GroomedJet_AK5_[%i]_GenWhad",i);
   DMass_GroomedJet_AK5_GenWhad[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DMass_GroomedJet_AK5_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());

   Name = Form("DPt_DR_GroomedJet_AK5_[%i]_GenWhad",i);
   DPt_DR_GroomedJet_AK5_GenWhad[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DPt_DR_GroomedJet_AK5_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DMass_DR_GroomedJet_AK5_[%i]_GenWhad",i);
   DMass_DR_GroomedJet_AK5_GenWhad[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DMass_DR_GroomedJet_AK5_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());


   Name = Form("DPt_GroomedJet_CA8_[%i]_GenWhad",i);
   DPt_GroomedJet_CA8_GenWhad[i]  = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DPt_GroomedJet_CA8_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DEta_GroomedJet_CA8_[%i]_GenWhad",i);
   DEta_GroomedJet_CA8_GenWhad[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-EtaMinMax,EtaMinMax);
   DEta_GroomedJet_CA8_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DPhi_GroomedJet_CA8_[%i]_GenWhad",i);
   DPhi_GroomedJet_CA8_GenWhad[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-PhiMinMax,PhiMinMax);
   DPhi_GroomedJet_CA8_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DR_GroomedJet_CA8_[%i]_GenWhad",i);
   DR_GroomedJet_CA8_GenWhad[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,DRMin,DRMax);
   DR_GroomedJet_CA8_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DE_GroomedJet_CA8_[%i]_GenWhad",i);
   DE_GroomedJet_CA8_GenWhad[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DE_GroomedJet_CA8_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str()); 
   Name = Form("DMass_GroomedJet_CA8_[%i]_GenWhad",i);
   DMass_GroomedJet_CA8_GenWhad[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DMass_GroomedJet_CA8_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());


   Name = Form("DPt_DR_GroomedJet_CA8_[%i]_GenWhad",i);
   DPt_DR_GroomedJet_CA8_GenWhad[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DPt_DR_GroomedJet_CA8_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DMass_DR_GroomedJet_CA8_[%i]_GenWhad",i);
   DMass_DR_GroomedJet_CA8_GenWhad[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DMass_DR_GroomedJet_CA8_GenWhad[i]->GetXaxis()->SetTitle(Name.c_str());
 

   Name = Form("DPt_GroomedJet_CA8_[%i]_subjet1_WParton",i);
   DPt_GroomedJet_CA8_subjet1_WParton[i]  = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DPt_GroomedJet_CA8_subjet1_WParton[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DEta_GroomedJet_CA8_[%i]_subjet1_WParton",i);
   DEta_GroomedJet_CA8_subjet1_WParton[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-EtaMinMax,EtaMinMax);
   DEta_GroomedJet_CA8_subjet1_WParton[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DPhi_GroomedJet_CA8_[%i]_subjet1_WParton",i);
   DPhi_GroomedJet_CA8_subjet1_WParton[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-PhiMinMax,PhiMinMax);
   DPhi_GroomedJet_CA8_subjet1_WParton[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DR_GroomedJet_CA8_[%i]_subjet1_WParton",i);
   DR_GroomedJet_CA8_subjet1_WParton[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,DRMin,DRMax);
   DR_GroomedJet_CA8_subjet1_WParton[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DE_GroomedJet_CA8_[%i]_subjet1_WParton",i);
   DE_GroomedJet_CA8_subjet1_WParton[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DE_GroomedJet_CA8_subjet1_WParton[i]->GetXaxis()->SetTitle(Name.c_str()); 
   Name = Form("DMass_GroomedJet_CA8_[%i]_subjet1_WParton",i);
   DMass_GroomedJet_CA8_subjet1_WParton[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DMass_GroomedJet_CA8_subjet1_WParton[i]->GetXaxis()->SetTitle(Name.c_str());


   Name = Form("DPt_DR_GroomedJet_CA8_[%i]_subjet1_WParton",i);
   DPt_DR_GroomedJet_CA8_subjet1_WParton[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DPt_DR_GroomedJet_CA8_subjet1_WParton[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DMass_DR_GroomedJet_CA8_[%i]_subjet1_WParton",i);
   DMass_DR_GroomedJet_CA8_subjet1_WParton[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DMass_DR_GroomedJet_CA8_subjet1_WParton[i]->GetXaxis()->SetTitle(Name.c_str());


   Name = Form("DPt_GroomedJet_CA8_[%i]_subjet2_WParton",i);
   DPt_GroomedJet_CA8_subjet2_WParton[i]  = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DPt_GroomedJet_CA8_subjet2_WParton[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DEta_GroomedJet_CA8_[%i]_subjet2_WParton",i);
   DEta_GroomedJet_CA8_subjet2_WParton[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-EtaMinMax,EtaMinMax);
   DEta_GroomedJet_CA8_subjet2_WParton[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DPhi_GroomedJet_CA8_[%i]_subjet2_WParton",i);
   DPhi_GroomedJet_CA8_subjet2_WParton[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-PhiMinMax,PhiMinMax);
   DPhi_GroomedJet_CA8_subjet2_WParton[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DR_GroomedJet_CA8_[%i]_subjet2_WParton",i);
   DR_GroomedJet_CA8_subjet2_WParton[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,DRMin,DRMax);
   DR_GroomedJet_CA8_subjet2_WParton[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DE_GroomedJet_CA8_[%i]_subjet2_WParton",i);
   DE_GroomedJet_CA8_subjet2_WParton[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DE_GroomedJet_CA8_subjet2_WParton[i]->GetXaxis()->SetTitle(Name.c_str()); 
   Name = Form("DMass_GroomedJet_CA8_[%i]_subjet2_WParton",i);
   DMass_GroomedJet_CA8_subjet2_WParton[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DMass_GroomedJet_CA8_subjet2_WParton[i]->GetXaxis()->SetTitle(Name.c_str());

   Name = Form("DPt_DR_GroomedJet_CA8_[%i]_subjet2_WParton",i);
   DPt_DR_GroomedJet_CA8_subjet2_WParton[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DPt_DR_GroomedJet_CA8_subjet2_WParton[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DMass_DR_GroomedJet_CA8_[%i]_subjet2_WParton",i);
   DMass_DR_GroomedJet_CA8_subjet2_WParton[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DMass_DR_GroomedJet_CA8_subjet2_WParton[i]->GetXaxis()->SetTitle(Name.c_str());

   Name = Form("DPt_GroomedJet_CA8_[%i]_GenWhad_pr",i);
   DPt_GroomedJet_CA8_GenWhad_pr[i]  = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DPt_GroomedJet_CA8_GenWhad_pr[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DEta_GroomedJet_CA8_[%i]_GenWhad_pr",i);
   DEta_GroomedJet_CA8_GenWhad_pr[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-EtaMinMax,EtaMinMax);
   DEta_GroomedJet_CA8_GenWhad_pr[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DPhi_GroomedJet_CA8_[%i]_GenWhad_pr",i);
   DPhi_GroomedJet_CA8_GenWhad_pr[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-PhiMinMax,PhiMinMax);
   DPhi_GroomedJet_CA8_GenWhad_pr[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DR_GroomedJet_CA8_[%i]_GenWhad_pr",i);
   DR_GroomedJet_CA8_GenWhad_pr[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,DRMin,DRMax);
   DR_GroomedJet_CA8_GenWhad_pr[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DE_GroomedJet_CA8_[%i]_GenWhad_pr",i);
   DE_GroomedJet_CA8_GenWhad_pr[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DE_GroomedJet_CA8_GenWhad_pr[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DMass_GroomedJet_CA8_[%i]_GenWhad_pr",i);
   DMass_GroomedJet_CA8_GenWhad_pr[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DMass_GroomedJet_CA8_GenWhad_pr[i]->GetXaxis()->SetTitle(Name.c_str());

   Name = Form("DPt_DR_GroomedJet_CA8_[%i]_GenWhad_pr",i);
   DPt_DR_GroomedJet_CA8_GenWhad_pr[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DPt_DR_GroomedJet_CA8_GenWhad_pr[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DMass_DR_GroomedJet_CA8_[%i]_GenWhad_pr",i);
   DMass_DR_GroomedJet_CA8_GenWhad_pr[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DMass_DR_GroomedJet_CA8_GenWhad_pr[i]->GetXaxis()->SetTitle(Name.c_str());
 
   Name = Form("DPt_GroomedJet_CA8_[%i]_GenWhad_tr",i);
   DPt_GroomedJet_CA8_GenWhad_tr[i]  = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DPt_GroomedJet_CA8_GenWhad_tr[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DEta_GroomedJet_CA8_[%i]_GenWhad_tr",i);
   DEta_GroomedJet_CA8_GenWhad_tr[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-EtaMinMax,EtaMinMax);
   DEta_GroomedJet_CA8_GenWhad_tr[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DPhi_GroomedJet_CA8_[%i]_GenWhad_tr",i);
   DPhi_GroomedJet_CA8_GenWhad_tr[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-PhiMinMax,PhiMinMax);
   DPhi_GroomedJet_CA8_GenWhad_tr[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DR_GroomedJet_CA8_[%i]_GenWhad_tr",i);
   DR_GroomedJet_CA8_GenWhad_tr[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,DRMin,DRMax);
   DR_GroomedJet_CA8_GenWhad_tr[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DE_GroomedJet_CA8_[%i]_GenWhad_tr",i);
   DE_GroomedJet_CA8_GenWhad_tr[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DE_GroomedJet_CA8_GenWhad_tr[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DMass_GroomedJet_CA8_[%i]_GenWhad_tr",i);
   DMass_GroomedJet_CA8_GenWhad_tr[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DMass_GroomedJet_CA8_GenWhad_tr[i]->GetXaxis()->SetTitle(Name.c_str());
  
   Name = Form("DPt_DR_GroomedJet_CA8_[%i]_GenWhad_tr",i);
   DPt_DR_GroomedJet_CA8_GenWhad_tr[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DPt_DR_GroomedJet_CA8_GenWhad_tr[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DMass_DR_GroomedJet_CA8_[%i]_GenWhad_tr",i);
   DMass_DR_GroomedJet_CA8_GenWhad_tr[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DMass_DR_GroomedJet_CA8_GenWhad_tr[i]->GetXaxis()->SetTitle(Name.c_str());
 
   Name = Form("DPt_GroomedJet_CA8_[%i]_GenWhad_ft",i);
   DPt_GroomedJet_CA8_GenWhad_ft[i]  = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DPt_GroomedJet_CA8_GenWhad_ft[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DEta_GroomedJet_CA8_[%i]_GenWhad_ft",i);
   DEta_GroomedJet_CA8_GenWhad_ft[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-EtaMinMax,EtaMinMax);
   DEta_GroomedJet_CA8_GenWhad_ft[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DPhi_GroomedJet_CA8_[%i]_GenWhad_ft",i);
   DPhi_GroomedJet_CA8_GenWhad_ft[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-PhiMinMax,PhiMinMax);
   DPhi_GroomedJet_CA8_GenWhad_ft[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DR_GroomedJet_CA8_[%i]_GenWhad_ft",i);
   DR_GroomedJet_CA8_GenWhad_ft[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,DRMin,DRMax);
   DR_GroomedJet_CA8_GenWhad_ft[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DE_GroomedJet_CA8_[%i]_GenWhad_ft",i);
   DE_GroomedJet_CA8_GenWhad_ft[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DE_GroomedJet_CA8_GenWhad_ft[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DMass_GroomedJet_CA8_[%i]_GenWhad_ft",i);
   DMass_GroomedJet_CA8_GenWhad_ft[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DMass_GroomedJet_CA8_GenWhad_ft[i]->GetXaxis()->SetTitle(Name.c_str());

   Name = Form("DPt_DR_GroomedJet_CA8_[%i]_GenWhad_ft",i);
   DPt_DR_GroomedJet_CA8_GenWhad_ft[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DPt_DR_GroomedJet_CA8_GenWhad_ft[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DMass_DR_GroomedJet_CA8_[%i]_GenWhad_ft",i);
   DMass_DR_GroomedJet_CA8_GenWhad_ft[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DMass_DR_GroomedJet_CA8_GenWhad_ft[i]->GetXaxis()->SetTitle(Name.c_str());
 
   Name = Form("DPt_GroomedJet_CA8_[%i]_WLep",i);
   DPt_GroomedJet_CA8_WLep[i]  = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DPt_GroomedJet_CA8_WLep[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DEta_GroomedJet_CA8_[%i]_WLep",i);
   DEta_GroomedJet_CA8_WLep[i]  = new TH1F (Name.c_str(),Name.c_str(),Nbin,-EtaMinMax,EtaMinMax);
   DEta_GroomedJet_CA8_WLep[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DPhi_GroomedJet_CA8_[%i]_WLep",i);
   DPhi_GroomedJet_CA8_WLep[i]  = new TH1F (Name.c_str(),Name.c_str(),Nbin,-PhiMinMax,PhiMinMax);
   DPhi_GroomedJet_CA8_WLep[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DR_GroomedJet_CA8_[%i]_WLep",i);
   DR_GroomedJet_CA8_WLep[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,DRMin,DRMax);
   DR_GroomedJet_CA8_WLep[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DE_GroomedJet_CA8_[%i]_WLep",i);
   DE_GroomedJet_CA8_WLep[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DE_GroomedJet_CA8_WLep[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DMass_GroomedJet_CA8_[%i]_WLep",i);
   DMass_GroomedJet_CA8_WLep[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
   DMass_GroomedJet_CA8_WLep[i]->GetXaxis()->SetTitle(Name.c_str());

   Name = Form("DPt_DR_GroomedJet_CA8_[%i]_WLep",i);
   DPt_DR_GroomedJet_CA8_WLep[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DPt_DR_GroomedJet_CA8_WLep[i]->GetXaxis()->SetTitle(Name.c_str());
   Name = Form("DMass_DR_GroomedJet_CA8_[%i]_WLep",i);
   DMass_DR_GroomedJet_CA8_WLep[i] = new TH2F (Name.c_str(),Name.c_str(),Nbin2D,PtEMassMin2D,PtEMassMax2D,Nbin2D,DRMin2D,DRMax2D);
   DMass_DR_GroomedJet_CA8_WLep[i]->GetXaxis()->SetTitle(Name.c_str());

  }

  TH1F*  DPt_WLep_GenWhad     = new TH1F ("DPt_GenWLep_GenWhad","DPt_GenWLep_GenWhad",Nbin,PtEMassMin,PtEMassMax);
  TH1F*  DEta_WLep_GenWhad    = new TH1F ("DEta_GenWLep_GenWhad","DEta_GenWLep_GenWhad",Nbin,-EtaMinMax,EtaMinMax);
  TH1F*  DPhi_WLep_GenWhad    = new TH1F ("DPhi_GenWLep_GenWhad","DPhi_GenWLep_GenWhad",Nbin,-PhiMinMax,PhiMinMax);
  TH1F*  DR_WLep_GenWhad      = new TH1F ("DR_GenWLep_GenWhad","DR_GenWLep_GenWhad",Nbin,DRMin,DRMax);
  TH1F*  DE_WLep_GenWhad      = new TH1F ("DE_GenWLep_GenWhad","DPhi_GenWLep_GenWhad",Nbin,PtEMassMin,PtEMassMax);
  TH1F*  DMass_WLep_GenWhad   = new TH1F ("DMass_GenWLep_GenWhad","DR_GenWLep_GenWhad",Nbin,PtEMassMin,PtEMassMax);

  // Histogramm after matching selections : inclusive categories

  TH1F**  DPt_GroomedCA8_GenWhad_DRCut      = new TH1F* [NJet];
  TH1F**  DEta_GroomedCA8_GenWhad_DRCut     = new TH1F* [NJet];
  TH1F**  DPhi_GroomedCA8_GenWhad_DRCut     = new TH1F* [NJet];
  TH1F**  DR_GroomedCA8_GenWhad_DRCut       = new TH1F* [NJet];
  TH1F**  DMass_GroomedCA8_GenWhad_DRCut    = new TH1F* [NJet];

  TH1F**  DPt_GroomedAK5_GenWhad_DRCut      = new TH1F* [NJet];
  TH1F**  DEta_GroomedAK5_GenWhad_DRCut     = new TH1F* [NJet];
  TH1F**  DPhi_GroomedAK5_GenWhad_DRCut     = new TH1F* [NJet];
  TH1F**  DR_GroomedAK5_GenWhad_DRCut       = new TH1F* [NJet];
  TH1F**  DMass_GroomedAK5_GenWhad_DRCut    = new TH1F* [NJet];

  TH1F**  DPt_GroomedCA8_GenWhad_DRCut_DPtCut      = new TH1F* [NJet];
  TH1F**  DEta_GroomedCA8_GenWhad_DRCut_DPtCut     = new TH1F* [NJet];
  TH1F**  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut     = new TH1F* [NJet];
  TH1F**  DR_GroomedCA8_GenWhad_DRCut_DPtCut       = new TH1F* [NJet];
  TH1F**  DMass_GroomedCA8_GenWhad_DRCut_DPtCut    = new TH1F* [NJet];

  TH1F**  DPt_GroomedAK5_GenWhad_DRCut_DPtCut      = new TH1F* [NJet];
  TH1F**  DEta_GroomedAK5_GenWhad_DRCut_DPtCut     = new TH1F* [NJet];
  TH1F**  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut     = new TH1F* [NJet];
  TH1F**  DR_GroomedAK5_GenWhad_DRCut_DPtCut       = new TH1F* [NJet];
  TH1F**  DMass_GroomedAK5_GenWhad_DRCut_DPtCut    = new TH1F* [NJet];


  for(int i = 0; i<NJet ; i++){

    TString name = Form("DPt_GroomedCA8_[%i]_GenWhad_DRCut",i);
    std::string Name = name.Data();

    DPt_GroomedCA8_GenWhad_DRCut[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
    DPt_GroomedCA8_GenWhad_DRCut[i]->GetXaxis()->SetTitle(Name.c_str());
    Name = Form("DEta_GroomedCA8_[%i]_GenWhad_DRCut",i);
    DEta_GroomedCA8_GenWhad_DRCut[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-EtaMinMax,EtaMinMax);
    DEta_GroomedCA8_GenWhad_DRCut[i]->GetXaxis()->SetTitle(Name.c_str());
    Name = Form("DPhi_GroomedCA8_[%i]_GenWhad_DRCut",i);
    DPhi_GroomedCA8_GenWhad_DRCut[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-PhiMinMax,PhiMinMax);
    DPhi_GroomedCA8_GenWhad_DRCut[i]->GetXaxis()->SetTitle(Name.c_str());
    Name = Form("DR_GroomedCA8_[%i]_GenWhad_DRCut",i);
    DR_GroomedCA8_GenWhad_DRCut[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,DRMin,DRMatchingCut);
    DR_GroomedCA8_GenWhad_DRCut[i]->GetXaxis()->SetTitle(Name.c_str());
    Name = Form("DMass_GroomedCA8_[%i]_GenWhad_DRCut",i);
    DMass_GroomedCA8_GenWhad_DRCut[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
    DMass_GroomedCA8_GenWhad_DRCut[i]->GetXaxis()->SetTitle(Name.c_str());

    Name = Form("DPt_GroomedAK5_[%i]_GenWhad_DRCut",i);
    DPt_GroomedAK5_GenWhad_DRCut[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
    DPt_GroomedAK5_GenWhad_DRCut[i]->GetXaxis()->SetTitle(Name.c_str());
    Name = Form("DEta_GroomedAK5_[%i]_GenWhad_DRCut",i);
    DEta_GroomedAK5_GenWhad_DRCut[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-EtaMinMax,EtaMinMax);
    DEta_GroomedAK5_GenWhad_DRCut[i]->GetXaxis()->SetTitle(Name.c_str());
    Name = Form("DPhi_GroomedAK5_[%i]_GenWhad_DRCut",i);
    DPhi_GroomedAK5_GenWhad_DRCut[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-PhiMinMax,PhiMinMax);
    DPhi_GroomedAK5_GenWhad_DRCut[i]->GetXaxis()->SetTitle(Name.c_str());
    Name = Form("DR_GroomedAK5_[%i]_GenWhad_DRCut",i);
    DR_GroomedAK5_GenWhad_DRCut[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,DRMin,DRMatchingCut);
    DR_GroomedAK5_GenWhad_DRCut[i]->GetXaxis()->SetTitle(Name.c_str());
    Name = Form("DMass_GroomedAK5_[%i]_GenWhad_DRCut",i);
    DMass_GroomedAK5_GenWhad_DRCut[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
    DMass_GroomedAK5_GenWhad_DRCut[i]->GetXaxis()->SetTitle(Name.c_str());


    Name = Form("DPt_GroomedCA8_[%i]_GenWhad_DRCut_DPtCut",i); 
    DPt_GroomedCA8_GenWhad_DRCut_DPtCut[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
    DPt_GroomedCA8_GenWhad_DRCut_DPtCut[i]->GetXaxis()->SetTitle(Name.c_str());
    Name = Form("DEta_GroomedCA8_[%i]_GenWhad_DRCut_DPtCut",i);
    DEta_GroomedCA8_GenWhad_DRCut_DPtCut[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-EtaMinMax,EtaMinMax);
    DEta_GroomedCA8_GenWhad_DRCut_DPtCut[i]->GetXaxis()->SetTitle(Name.c_str());
    Name = Form("DPhi_GroomedCA8_[%i]_GenWhad_DRCut_DPtCut",i);
    DPhi_GroomedCA8_GenWhad_DRCut_DPtCut[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-PhiMinMax,PhiMinMax);
    DPhi_GroomedCA8_GenWhad_DRCut_DPtCut[i]->GetXaxis()->SetTitle(Name.c_str());
    Name = Form("DR_GroomedCA8_[%i]_GenWhad_DRCut_DPtCut",i);
    DR_GroomedCA8_GenWhad_DRCut_DPtCut[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,DRMin,DRMatchingCut);
    DR_GroomedCA8_GenWhad_DRCut_DPtCut[i]->GetXaxis()->SetTitle(Name.c_str());
    Name = Form("DMass_GroomedCA8_[%i]_GenWhad_DRCut_DPtCut",i);
    DMass_GroomedCA8_GenWhad_DRCut_DPtCut[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
    DMass_GroomedCA8_GenWhad_DRCut_DPtCut[i]->GetXaxis()->SetTitle(Name.c_str());

    Name = Form("DPt_GroomedAK5_[%i]_GenWhad_DRCut_DPtCut",i);
    DPt_GroomedAK5_GenWhad_DRCut_DPtCut[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
    DPt_GroomedAK5_GenWhad_DRCut_DPtCut[i]->GetXaxis()->SetTitle(Name.c_str());
    Name = Form("DEta_GroomedAK5_[%i]_GenWhad_DRCut_DPtCut",i);
    DEta_GroomedAK5_GenWhad_DRCut_DPtCut[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-EtaMinMax,EtaMinMax);
    DEta_GroomedAK5_GenWhad_DRCut_DPtCut[i]->GetXaxis()->SetTitle(Name.c_str());
    Name = Form("DPhi_GroomedAK5_[%i]_GenWhad_DRCut_DPtCut",i);
    DPhi_GroomedAK5_GenWhad_DRCut_DPtCut[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,-PhiMinMax,PhiMinMax);
    DPhi_GroomedAK5_GenWhad_DRCut_DPtCut[i]->GetXaxis()->SetTitle(Name.c_str());
    Name = Form("DR_GroomedAK5_[%i]_GenWhad_DRCut_DPtCut",i);
    DR_GroomedAK5_GenWhad_DRCut_DPtCut[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,DRMin,DRMatchingCut);
    DR_GroomedAK5_GenWhad_DRCut_DPtCut[i]->GetXaxis()->SetTitle(Name.c_str());
    Name = Form("DMass_GroomedAK5_[%i]_GenWhad_DRCut_DPtCut",i);
    DMass_GroomedAK5_GenWhad_DRCut_DPtCut[i] = new TH1F (Name.c_str(),Name.c_str(),Nbin,PtEMassMin,PtEMassMax);
    DMass_GroomedAK5_GenWhad_DRCut_DPtCut[i]->GetXaxis()->SetTitle(Name.c_str());


  }

  //   Histogramm after matching selections : Exclusive 1 CA8 Jet 

  TH1F*  DPt_GroomedCA8_GenWhad_DRCut_1Jet    = new TH1F ("DPt_GroomedCA8_GenWhad_DRCut_1Jet","DPt_GroomedCA8_GenWhad_DRCut_1Jet",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedCA8_GenWhad_DRCut_1Jet->GetXaxis()->SetTitle("DPt_GroomedCA8_GenWhad_DRCut_1Jet");
  TH1F*  DEta_GroomedCA8_GenWhad_DRCut_1Jet   = new TH1F ("DEta_GroomedCA8_GenWhad_DRCut_1Jet","DEta_GroomedCA8_GenWhad_DRCut_1Jet",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedCA8_GenWhad_DRCut_1Jet->GetXaxis()->SetTitle("DEta_GroomedCA8_GenWhad_DRCut_1Jet");
  TH1F*  DPhi_GroomedCA8_GenWhad_DRCut_1Jet   = new TH1F ("DPhi_GroomedCA8_GenWhad_DRCut_1Jet","DPhi_GroomedCA8_GenWhad_DRCut_1Jet",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedCA8_GenWhad_DRCut_1Jet->GetXaxis()->SetTitle("DPhi_GroomedCA8_GenWhad_DRCut_1Jet");
  TH1F*  DR_GroomedCA8_GenWhad_DRCut_1Jet     = new TH1F ("DR_GroomedCA8_GenWhad_DRCut_1Jet","DR_GroomedCA8_GenWhad_DRCut_1Jet",Nbin,DRMin,DRMatchingCut);
  DR_GroomedCA8_GenWhad_DRCut_1Jet->GetXaxis()->SetTitle("DR_GroomedCA8_GenWhad_DRCut_1Jet");
  TH1F*  DMass_GroomedCA8_GenWhad_DRCut_1Jet  = new TH1F ("DMass_GroomedCA8_GenWhad_DRCut_1Jet","DMass_GroomedCA8_GenWhad_DRCut_1Jet",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedCA8_GenWhad_DRCut_1Jet->GetXaxis()->SetTitle("DMass_GroomedCA8_GenWhad_DRCut_1Jet");

  TH1F*  DPt_GroomedAK5_GenWhad_DRCut_1Jet    = new TH1F ("DPt_GroomedAK5_GenWhad_DRCut_1Jet","DPt_GroomedAK5_GenWhad_DRCut_1Jet",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedAK5_GenWhad_DRCut_1Jet->GetXaxis()->SetTitle("DPt_GroomedAK5_GenWhad_DRCut_1Jet");
  TH1F*  DEta_GroomedAK5_GenWhad_DRCut_1Jet   = new TH1F ("DEta_GroomedAK5_GenWhad_DRCut_1Jet","DEta_GroomedAK5_GenWhad_DRCut_1Jet",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedAK5_GenWhad_DRCut_1Jet->GetXaxis()->SetTitle("DEta_GroomedAK5_GenWhad_DRCut_1Jet");
  TH1F*  DPhi_GroomedAK5_GenWhad_DRCut_1Jet   = new TH1F ("DPhi_GroomedAK5_GenWhad_DRCut_1Jet","DPhi_GroomedAK5_GenWhad_DRCut_1Jet",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedAK5_GenWhad_DRCut_1Jet->GetXaxis()->SetTitle("DPhi_GroomedAK5_GenWhad_DRCut_1Jet");
  TH1F*  DR_GroomedAK5_GenWhad_DRCut_1Jet     = new TH1F ("DR_GroomedAK5_GenWhad_DRCut_1Jet","DR_GroomedAK5_GenWhad_DRCut_1Jet",Nbin,DRMin,DRMatchingCut);
  DR_GroomedAK5_GenWhad_DRCut_1Jet->GetXaxis()->SetTitle("DR_GroomedAK5_GenWhad_DRCut_1Jet");
  TH1F*  DMass_GroomedAK5_GenWhad_DRCut_1Jet  = new TH1F ("DMass_GroomedAK5_GenWhad_DRCut_1Jet","DMass_GroomedAK5_GenWhad_DRCut_1Jet",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedAK5_GenWhad_DRCut_1Jet->GetXaxis()->SetTitle("DMass_GroomedAK5_GenWhad_DRCut_1Jet");

  TH1F*  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet    = new TH1F ("DPt_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet","DPt_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet->GetXaxis()->SetTitle("DPt_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet");
  TH1F*  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet   = new TH1F ("DEta_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet","DEta_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet->GetXaxis()->SetTitle("DEta_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet");
  TH1F*  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet   = new TH1F ("DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet","DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet->GetXaxis()->SetTitle("DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet");
  TH1F*  DR_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet     = new TH1F ("DR_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet","DR_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet",Nbin,DRMin,DRMatchingCut);
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet->GetXaxis()->SetTitle("DR_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet");
  TH1F*  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet  = new TH1F ("DMass_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet","DMass_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet->GetXaxis()->SetTitle("DMass_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet");

  TH1F*  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet    = new TH1F ("DPt_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet","DPt_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->GetXaxis()->SetTitle("DPt_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet");
  TH1F*  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet   = new TH1F ("DEta_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet","DEta_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->GetXaxis()->SetTitle("DEta_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet");
  TH1F*  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet   = new TH1F ("DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet","DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->GetXaxis()->SetTitle("DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet");
  TH1F*  DR_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet     = new TH1F ("DR_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet","DR_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet",Nbin,DRMin,DRMatchingCut);
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->GetXaxis()->SetTitle("DR_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet");
  TH1F*  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet  = new TH1F ("DMass_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet","DMass_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->GetXaxis()->SetTitle("DMass_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet");

  //   Histogramm after matching selections : Exclusive 2 CA8 Jet --> hard Matched, II outside

  TH1F*  DPt_GroomedCA8_GenWhad_DRCut_2Jet    = new TH1F ("DPt_GroomedCA8_GenWhad_DRCut_2Jet","DPt_GroomedCA8_GenWhad_DRCut_2Jet",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedCA8_GenWhad_DRCut_2Jet->GetXaxis()->SetTitle("DPt_GroomedCA8_GenWhad_DRCut_2Jet");
  TH1F*  DEta_GroomedCA8_GenWhad_DRCut_2Jet   = new TH1F ("DEta_GroomedCA8_GenWhad_DRCut_2Jet","DEta_GroomedCA8_GenWhad_DRCut_2Jet",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedCA8_GenWhad_DRCut_2Jet->GetXaxis()->SetTitle("DEta_GroomedCA8_GenWhad_DRCut_2Jet");
  TH1F*  DPhi_GroomedCA8_GenWhad_DRCut_2Jet   = new TH1F ("DPhi_GroomedCA8_GenWhad_DRCut_2Jet","DPhi_GroomedCA8_GenWhad_DRCut_2Jet",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedCA8_GenWhad_DRCut_2Jet->GetXaxis()->SetTitle("DPhi_GroomedCA8_GenWhad_DRCut_2Jet");
  TH1F*  DR_GroomedCA8_GenWhad_DRCut_2Jet     = new TH1F ("DR_GroomedCA8_GenWhad_DRCut_2Jet","DR_GroomedCA8_GenWhad_DRCut_2Jet",Nbin,DRMin,DRMatchingCut);
  DR_GroomedCA8_GenWhad_DRCut_2Jet->GetXaxis()->SetTitle("DR_GroomedCA8_GenWhad_DRCut_2Jet");
  TH1F*  DMass_GroomedCA8_GenWhad_DRCut_2Jet  = new TH1F ("DMass_GroomedCA8_GenWhad_DRCut_2Jet","DMass_GroomedCA8_GenWhad_DRCut_2Jet",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedCA8_GenWhad_DRCut_2Jet->GetXaxis()->SetTitle("DMass_GroomedCA8_GenWhad_DRCut_2Jet");

  TH1F*  DPt_GroomedAK5_GenWhad_DRCut_2Jet    = new TH1F ("DPt_GroomedAK5_GenWhad_DRCut_2Jet","DPt_GroomedAK5_GenWhad_DRCut_2Jet",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedAK5_GenWhad_DRCut_2Jet->GetXaxis()->SetTitle("DPt_GroomedAK5_GenWhad_DRCut_2Jet");
  TH1F*  DEta_GroomedAK5_GenWhad_DRCut_2Jet   = new TH1F ("DEta_GroomedAK5_GenWhad_DRCut_2Jet","DEta_GroomedAK5_GenWhad_DRCut_2Jet",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedAK5_GenWhad_DRCut_2Jet->GetXaxis()->SetTitle("DEta_GroomedAK5_GenWhad_DRCut_2Jet");
  TH1F*  DPhi_GroomedAK5_GenWhad_DRCut_2Jet   = new TH1F ("DPhi_GroomedAK5_GenWhad_DRCut_2Jet","DPhi_GroomedAK5_GenWhad_DRCut_2Jet",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedAK5_GenWhad_DRCut_2Jet->GetXaxis()->SetTitle("DPhi_GroomedAK5_GenWhad_DRCut_2Jet");
  TH1F*  DR_GroomedAK5_GenWhad_DRCut_2Jet     = new TH1F ("DR_GroomedAK5_GenWhad_DRCut_2Jet","DR_GroomedAK5_GenWhad_DRCut_2Jet",Nbin,DRMin,DRMatchingCut);
  DR_GroomedAK5_GenWhad_DRCut_2Jet->GetXaxis()->SetTitle("DR_GroomedAK5_GenWhad_DRCut_2Jet");
  TH1F*  DMass_GroomedAK5_GenWhad_DRCut_2Jet  = new TH1F ("DMass_GroomedAK5_GenWhad_DRCut_2Jet","DMass_GroomedAK5_GenWhad_DRCut_2Jet",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedAK5_GenWhad_DRCut_2Jet->GetXaxis()->SetTitle("DMass_GroomedAK5_GenWhad_DRCut_2Jet");

  TH1F*  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet    = new TH1F ("DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet","DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->GetXaxis()->SetTitle("DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet");
  TH1F*  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet   = new TH1F ("DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet","DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->GetXaxis()->SetTitle("DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet");
  TH1F*  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet   = new TH1F ("DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet","DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->GetXaxis()->SetTitle("DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet");
  TH1F*  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet     = new TH1F ("DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet","DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet",Nbin,DRMin,DRMatchingCut);
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->GetXaxis()->SetTitle("DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet");
  TH1F*  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet  = new TH1F ("DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet","DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->GetXaxis()->SetTitle("DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet");

  TH1F*  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet    = new TH1F ("DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet","DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->GetXaxis()->SetTitle("DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet");
  TH1F*  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet   = new TH1F ("DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet","DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->GetXaxis()->SetTitle("DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet");
  TH1F*  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet   = new TH1F ("DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet","DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->GetXaxis()->SetTitle("DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet");
  TH1F*  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet     = new TH1F ("DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet","DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet",Nbin,DRMin,DRMatchingCut);
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->GetXaxis()->SetTitle("DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet");
  TH1F*  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet  = new TH1F ("DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet","DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->GetXaxis()->SetTitle("DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet");

  // Histogramm for CA8 jet matching simultaneusly with another CA8 jet in exclusive 2 jet events

  TH1F*  DPt_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching    = new TH1F ("DPt_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching",
                                                                         "DPt_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->GetXaxis()->SetTitle("DPt_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching");
  TH1F*  DEta_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching   = new TH1F ("DEta_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching",
                                                                          "DEta_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->GetXaxis()->SetTitle("DEta_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching");
  TH1F*  DPhi_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching   = new TH1F ("DPhi_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching",
                                                                         "DPhi_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->GetXaxis()->SetTitle("DPhi_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching");
  TH1F*  DR_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching     = new TH1F ("DR_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching",
                                                                         "DR_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching",Nbin,DRMin,DRMatchingCut);
  DR_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->GetXaxis()->SetTitle("DR_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching");
  TH1F*  DMass_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching  = new TH1F ("DMass_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching",
                                                                         "DMass_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->GetXaxis()->SetTitle("DMass_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching");

  TH1F*  DPt_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching    = new TH1F ("DPt_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching",
                                                                         "DPt_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->GetXaxis()->SetTitle("DPt_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching");
  TH1F*  DEta_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching   = new TH1F ("DEta_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching",
                                                                         "DEta_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->GetXaxis()->SetTitle("DEta_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching");
  TH1F*  DPhi_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching   = new TH1F ("DPhi_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching",
                                                                         "DPhi_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->GetXaxis()->SetTitle("DPhi_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching");
  TH1F*  DR_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching     = new TH1F ("DR_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching",
                                                                         "DR_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching",Nbin,DRMin,DRMatchingCut);
  DR_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->GetXaxis()->SetTitle("DR_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching");
  TH1F*  DMass_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching  = new TH1F ("DMass_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching",
                                                                         "DMass_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->GetXaxis()->SetTitle("DMass_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching");

  TH1F*  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching   = new TH1F ("DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching",
                                                                               "DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetXaxis()->SetTitle("DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching");
  TH1F*  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching   = new TH1F ("DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching",
                                                                                "DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetXaxis()->SetTitle("DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching");
  TH1F*  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching   = new TH1F ("DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching",
                                                                                "DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetXaxis()->SetTitle("DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching");
  TH1F*  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching     = new TH1F ("DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching",
                                                                                "DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching",Nbin,DRMin,DRMatchingCut);
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetXaxis()->SetTitle("DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching");
  TH1F*  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching  = new TH1F ("DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching",
                                                                                "DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetXaxis()->SetTitle("DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching");

  TH1F*  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching    = new TH1F ("DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching",
                                                                                "DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetXaxis()->SetTitle("DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching");
  TH1F*  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching   = new TH1F ("DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching",
                                                                                "DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetXaxis()->SetTitle("DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching");
  TH1F*  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching   = new TH1F ("DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching",
                                                                                "DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetXaxis()->SetTitle("DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching");
  TH1F*  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching     = new TH1F ("DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching",
                                                                                "DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching",Nbin,DRMin,DRMatchingCut);
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetXaxis()->SetTitle("DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching");
  TH1F*  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching  = new TH1F ("DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching",
                                                                                "DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetXaxis()->SetTitle("DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching");

  // Histogramm for CA8 jet matching inverted in exclusive 2 jet events

  TH1F*  DPt_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching    = new TH1F ("DPt_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching",
                                                                         "DPt_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->GetXaxis()->SetTitle("DPt_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching");
  TH1F*  DEta_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching   = new TH1F ("DEta_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching",
                                                                          "DEta_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->GetXaxis()->SetTitle("DEta_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching");
  TH1F*  DPhi_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching   = new TH1F ("DPhi_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching",
                                                                         "DPhi_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->GetXaxis()->SetTitle("DPhi_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching");
  TH1F*  DR_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching     = new TH1F ("DR_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching",
                                                                         "DR_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching",Nbin,DRMin,DRMatchingCut);
  DR_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->GetXaxis()->SetTitle("DR_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching");
  TH1F*  DMass_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching  = new TH1F ("DMass_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching",
                                                                         "DMass_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->GetXaxis()->SetTitle("DMass_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching");

  TH1F*  DPt_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching    = new TH1F ("DPt_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching",
                                                                         "DPt_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->GetXaxis()->SetTitle("DPt_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching");
  TH1F*  DEta_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching   = new TH1F ("DEta_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching",
                                                                         "DEta_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->GetXaxis()->SetTitle("DEta_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching");
  TH1F*  DPhi_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching   = new TH1F ("DPhi_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching",
                                                                         "DPhi_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->GetXaxis()->SetTitle("DPhi_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching");
  TH1F*  DR_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching     = new TH1F ("DR_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching",
                                                                         "DR_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching",Nbin,DRMin,DRMatchingCut);
  DR_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->GetXaxis()->SetTitle("DR_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching");
  TH1F*  DMass_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching  = new TH1F ("DMass_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching",
                                                                         "DMass_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->GetXaxis()->SetTitle("DMass_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching");

  TH1F*  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching   = new TH1F ("DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching",
                                                                               "DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetXaxis()->SetTitle("DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching");
  TH1F*  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching   = new TH1F ("DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching",
                                                                                "DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetXaxis()->SetTitle("DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching");
  TH1F*  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching   = new TH1F ("DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching",
                                                                                "DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetXaxis()->SetTitle("DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching");
  TH1F*  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching     = new TH1F ("DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching",
                                                                                "DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching",Nbin,DRMin,DRMatchingCut);
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetXaxis()->SetTitle("DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching");
  TH1F*  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching  = new TH1F ("DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching",
                                                                                "DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetXaxis()->SetTitle("DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching");

  TH1F*  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching    = new TH1F ("DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching",
                                                                                "DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetXaxis()->SetTitle("DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching");
  TH1F*  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching   = new TH1F ("DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching",
                                                                                "DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetXaxis()->SetTitle("DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching");
  TH1F*  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching   = new TH1F ("DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching",
                                                                                "DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetXaxis()->SetTitle("DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching");
  TH1F*  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching     = new TH1F ("DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching",
                                                                                "DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching",Nbin,DRMin,DRMatchingCut);
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetXaxis()->SetTitle("DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching");
  TH1F*  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching  = new TH1F ("DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching",
                                                                                "DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetXaxis()->SetTitle("DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching");


  // Histogramm for CA8 jet no matching in exclusive 2 jet events

  TH1F*  DPt_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching    = new TH1F ("DPt_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching",
                                                                         "DPt_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->GetXaxis()->SetTitle("DPt_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching");
  TH1F*  DEta_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching   = new TH1F ("DEta_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching",
                                                                          "DEta_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->GetXaxis()->SetTitle("DEta_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching");
  TH1F*  DPhi_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching   = new TH1F ("DPhi_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching",
                                                                         "DPhi_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->GetXaxis()->SetTitle("DPhi_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching");
  TH1F*  DR_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching     = new TH1F ("DR_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching",
                                                                         "DR_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching",Nbin,DRMin,DRMax);
  DR_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->GetXaxis()->SetTitle("DR_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching");
  TH1F*  DMass_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching  = new TH1F ("DMass_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching",
                                                                         "DMass_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->GetXaxis()->SetTitle("DMass_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching");

  TH1F*  DPt_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching    = new TH1F ("DPt_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching",
                                                                         "DPt_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->GetXaxis()->SetTitle("DPt_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching");
  TH1F*  DEta_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching   = new TH1F ("DEta_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching",
                                                                         "DEta_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->GetXaxis()->SetTitle("DEta_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching");
  TH1F*  DPhi_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching   = new TH1F ("DPhi_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching",
                                                                         "DPhi_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->GetXaxis()->SetTitle("DPhi_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching");
  TH1F*  DR_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching     = new TH1F ("DR_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching",
                                                                         "DR_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching",Nbin,DRMin,DRMax);
  DR_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->GetXaxis()->SetTitle("DR_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching");
  TH1F*  DMass_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching  = new TH1F ("DMass_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching",
                                                                         "DMass_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->GetXaxis()->SetTitle("DMass_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching");

  TH1F*  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching   = new TH1F ("DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching",
                                                                               "DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetXaxis()->SetTitle("DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching");
  TH1F*  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching   = new TH1F ("DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching",
                                                                                "DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetXaxis()->SetTitle("DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching");
  TH1F*  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching   = new TH1F ("DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching",
                                                                                "DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetXaxis()->SetTitle("DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching");
  TH1F*  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching     = new TH1F ("DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching",
                                                                                "DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching",Nbin,DRMin,DRMax);
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetXaxis()->SetTitle("DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching");
  TH1F*  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching  = new TH1F ("DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching",
                                                                                "DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetXaxis()->SetTitle("DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching");

  TH1F*  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching    = new TH1F ("DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching",
                                                                                "DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetXaxis()->SetTitle("DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching");
  TH1F*  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching   = new TH1F ("DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching",
                                                                                "DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetXaxis()->SetTitle("DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching");
  TH1F*  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching   = new TH1F ("DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching",
                                                                                "DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetXaxis()->SetTitle("DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching");
  TH1F*  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching     = new TH1F ("DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching",
                                                                                "DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching",Nbin,DRMin,DRMax);
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetXaxis()->SetTitle("DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching");
  TH1F*  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching  = new TH1F ("DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching",
                                                                                "DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetXaxis()->SetTitle("DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching");


  //   Histogramm after matching selections : Exclusive 3 CA8 Jet 

  TH1F*  DPt_GroomedCA8_GenWhad_DRCut_3Jet    = new TH1F ("DPt_GroomedCA8_GenWhad_DRCut_3Jet","DPt_GroomedCA8_GenWhad_DRCut_3Jet",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedCA8_GenWhad_DRCut_3Jet->GetXaxis()->SetTitle("DPt_GroomedCA8_GenWhad_DRCut_3Jet");
  TH1F*  DEta_GroomedCA8_GenWhad_DRCut_3Jet   = new TH1F ("DEta_GroomedCA8_GenWhad_DRCut_3Jet","DEta_GroomedCA8_GenWhad_DRCut_3Jet",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedCA8_GenWhad_DRCut_3Jet->GetXaxis()->SetTitle("DEta_GroomedCA8_GenWhad_DRCut_3Jet");
  TH1F*  DPhi_GroomedCA8_GenWhad_DRCut_3Jet   = new TH1F ("DPhi_GroomedCA8_GenWhad_DRCut_3Jet","DPhi_GroomedCA8_GenWhad_DRCut_3Jet",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedCA8_GenWhad_DRCut_3Jet->GetXaxis()->SetTitle("DPhi_GroomedCA8_GenWhad_DRCut_3Jet");
  TH1F*  DR_GroomedCA8_GenWhad_DRCut_3Jet     = new TH1F ("DR_GroomedCA8_GenWhad_DRCut_3Jet","DR_GroomedCA8_GenWhad_DRCut_3Jet",Nbin,DRMin,DRMatchingCut);
  DR_GroomedCA8_GenWhad_DRCut_3Jet->GetXaxis()->SetTitle("DR_GroomedCA8_GenWhad_DRCut_3Jet");
  TH1F*  DMass_GroomedCA8_GenWhad_DRCut_3Jet  = new TH1F ("DMass_GroomedCA8_GenWhad_DRCut_3Jet","DMass_GroomedCA8_GenWhad_DRCut_3Jet",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedCA8_GenWhad_DRCut_3Jet->GetXaxis()->SetTitle("DMass_GroomedCA8_GenWhad_DRCut_3Jet");

  TH1F*  DPt_GroomedAK5_GenWhad_DRCut_3Jet    = new TH1F ("DPt_GroomedAK5_GenWhad_DRCut_3Jet","DPt_GroomedAK5_GenWhad_DRCut_3Jet",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedAK5_GenWhad_DRCut_3Jet->GetXaxis()->SetTitle("DPt_GroomedAK5_GenWhad_DRCut_3Jet");
  TH1F*  DEta_GroomedAK5_GenWhad_DRCut_3Jet   = new TH1F ("DEta_GroomedAK5_GenWhad_DRCut_3Jet","DEta_GroomedAK5_GenWhad_DRCut_3Jet",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedAK5_GenWhad_DRCut_3Jet->GetXaxis()->SetTitle("DEta_GroomedAK5_GenWhad_DRCut_3Jet");
  TH1F*  DPhi_GroomedAK5_GenWhad_DRCut_3Jet   = new TH1F ("DPhi_GroomedAK5_GenWhad_DRCut_3Jet","DPhi_GroomedAK5_GenWhad_DRCut_3Jet",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedAK5_GenWhad_DRCut_3Jet->GetXaxis()->SetTitle("DPhi_GroomedAK5_GenWhad_DRCut_3Jet");
  TH1F*  DR_GroomedAK5_GenWhad_DRCut_3Jet     = new TH1F ("DR_GroomedAK5_GenWhad_DRCut_3Jet","DR_GroomedAK5_GenWhad_DRCut_3Jet",Nbin,DRMin,DRMatchingCut);
  DR_GroomedAK5_GenWhad_DRCut_3Jet->GetXaxis()->SetTitle("DR_GroomedAK5_GenWhad_DRCut_3Jet");
  TH1F*  DMass_GroomedAK5_GenWhad_DRCut_3Jet  = new TH1F ("DMass_GroomedAK5_GenWhad_DRCut_3Jet","DMass_GroomedAK5_GenWhad_DRCut_3Jet",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedAK5_GenWhad_DRCut_3Jet->GetXaxis()->SetTitle("DMass_GroomedAK5_GenWhad_DRCut_3Jet");

  TH1F*  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet    = new TH1F ("DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet","DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet->GetXaxis()->SetTitle("DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet");
  TH1F*  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet   = new TH1F ("DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet","DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet->GetXaxis()->SetTitle("DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet");
  TH1F*  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet   = new TH1F ("DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet","DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet->GetXaxis()->SetTitle("DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet");
  TH1F*  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet     = new TH1F ("DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet","DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet",Nbin,DRMin,DRMatchingCut);
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet->GetXaxis()->SetTitle("DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet");
  TH1F*  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet  = new TH1F ("DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet","DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet->GetXaxis()->SetTitle("DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet");

  TH1F*  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet    = new TH1F ("DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet","DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet->GetXaxis()->SetTitle("DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet");
  TH1F*  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet   = new TH1F ("DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet","DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet->GetXaxis()->SetTitle("DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet");
  TH1F*  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet   = new TH1F ("DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet","DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet->GetXaxis()->SetTitle("DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet");
  TH1F*  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet     = new TH1F ("DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet","DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet",Nbin,DRMin,DRMatchingCut);
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet->GetXaxis()->SetTitle("DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet");
  TH1F*  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet  = new TH1F ("DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet","DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet->GetXaxis()->SetTitle("DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet");


  //   Histogramm after matching selections : Double Matching in event with three jets

  TH1F*  DPt_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching    = new TH1F ("DPt_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching",
                                                                         "DPt_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->GetXaxis()->SetTitle("DPt_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching");
  TH1F*  DEta_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching   = new TH1F ("DEta_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching",
                                                                         "DEta_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->GetXaxis()->SetTitle("DEta_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching");
  TH1F*  DPhi_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching   = new TH1F ("DPhi_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching",
                                                                         "DPhi_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->GetXaxis()->SetTitle("DPhi_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching");
  TH1F*  DR_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching     = new TH1F ("DR_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching",
                                                                         "DR_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching",Nbin,DRMin,DRMatchingCut);
  DR_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->GetXaxis()->SetTitle("DR_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching");
  TH1F*  DMass_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching  = new TH1F ("DMass_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching",
                                                                         "DMass_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->GetXaxis()->SetTitle("DMass_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching");
  TH1F*  DPt_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching    = new TH1F ("DPt_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching",
                                                                         "DPt_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching",Nbin,PtEMassMin,PtEMassMax);

  DPt_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->GetXaxis()->SetTitle("DPt_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching");
  TH1F*  DEta_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching   = new TH1F ("DEta_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching",
                                                                         "DEta_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->GetXaxis()->SetTitle("DEta_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching");
  TH1F*  DPhi_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching   = new TH1F ("DPhi_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching",
                                                                         "DPhi_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->GetXaxis()->SetTitle("DPhi_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching");
  TH1F*  DR_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching     = new TH1F ("DR_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching",
                                                                         "DR_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching",Nbin,DRMin,DRMatchingCut);
  DR_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->GetXaxis()->SetTitle("DR_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching");
  TH1F*  DMass_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching  = new TH1F ("DMass_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching",
                                                                         "DMass_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->GetXaxis()->SetTitle("DMass_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching");
  TH1F*  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching    = new TH1F ("DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching",
                                                                                "DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching",Nbin,PtEMassMin,PtEMassMax);

  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetXaxis()->SetTitle("DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching");
  TH1F*  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching   = new TH1F ("DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching",
                                                                                "DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetXaxis()->SetTitle("DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching");
  TH1F*  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching   = new TH1F ("DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching",
                                                                                "DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetXaxis()->SetTitle("DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching");
  TH1F*  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching     = new TH1F ("DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching",
                                                                                "DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching",Nbin,DRMin,DRMatchingCut);
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetXaxis()->SetTitle("DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching");
  TH1F*  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching  = new TH1F ("DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching",
                                                                                "DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetXaxis()->SetTitle("DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching");

  TH1F*  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching    = new TH1F ("DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching",
                                                                                "DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetXaxis()->SetTitle("DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching");
  TH1F*  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching   = new TH1F ("DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching",
                                                                                "DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetXaxis()->SetTitle("DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching");
  TH1F*  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching   = new TH1F ("DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching",
                                                                                "DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetXaxis()->SetTitle("DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching");
  TH1F*  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching     = new TH1F ("DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching",
                                                                                "DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching",Nbin,DRMin,DRMatchingCut);
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetXaxis()->SetTitle("DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching");
  TH1F*  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching  = new TH1F ("DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching",
                                                                                "DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetXaxis()->SetTitle("DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching");


  //   Histogramm after matching selections : Double Matching in event with three jets

  TH1F*  DPt_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching    = new TH1F ("DPt_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching",
                                                                         "DPt_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->GetXaxis()->SetTitle("DPt_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching");
  TH1F*  DEta_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching   = new TH1F ("DEta_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching",
                                                                         "DEta_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->GetXaxis()->SetTitle("DEta_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching");
  TH1F*  DPhi_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching   = new TH1F ("DPhi_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching",
                                                                         "DPhi_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->GetXaxis()->SetTitle("DPhi_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching");
  TH1F*  DR_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching     = new TH1F ("DR_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching",
                                                                         "DR_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching",Nbin,DRMin,DRMatchingCut);
  DR_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->GetXaxis()->SetTitle("DR_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching");
  TH1F*  DMass_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching  = new TH1F ("DMass_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching",
                                                                         "DMass_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->GetXaxis()->SetTitle("DMass_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching");
  TH1F*  DPt_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching    = new TH1F ("DPt_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching",
                                                                         "DPt_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching",Nbin,PtEMassMin,PtEMassMax);

  DPt_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->GetXaxis()->SetTitle("DPt_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching");
  TH1F*  DEta_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching   = new TH1F ("DEta_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching",
                                                                         "DEta_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->GetXaxis()->SetTitle("DEta_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching");
  TH1F*  DPhi_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching   = new TH1F ("DPhi_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching",
                                                                         "DPhi_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->GetXaxis()->SetTitle("DPhi_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching");
  TH1F*  DR_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching     = new TH1F ("DR_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching",
                                                                         "DR_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching",Nbin,DRMin,DRMatchingCut);
  DR_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->GetXaxis()->SetTitle("DR_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching");
  TH1F*  DMass_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching  = new TH1F ("DMass_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching",
                                                                         "DMass_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->GetXaxis()->SetTitle("DMass_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching");
  TH1F*  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching    = new TH1F ("DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching",
                                                                                "DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching",Nbin,PtEMassMin,PtEMassMax);

  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetXaxis()->SetTitle("DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching");
  TH1F*  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching   = new TH1F ("DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching",
                                                                                "DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetXaxis()->SetTitle("DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching");
  TH1F*  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching   = new TH1F ("DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching",
                                                                                "DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetXaxis()->SetTitle("DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching");
  TH1F*  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching     = new TH1F ("DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching",
                                                                                "DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching",Nbin,DRMin,DRMatchingCut);
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetXaxis()->SetTitle("DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching");
  TH1F*  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching  = new TH1F ("DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching",
                                                                                "DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetXaxis()->SetTitle("DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching");

  TH1F*  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching    = new TH1F ("DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching",
                                                                                "DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetXaxis()->SetTitle("DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching");
  TH1F*  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching   = new TH1F ("DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching",
                                                                                "DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetXaxis()->SetTitle("DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching");
  TH1F*  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching   = new TH1F ("DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching",
                                                                                "DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetXaxis()->SetTitle("DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching");
  TH1F*  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching     = new TH1F ("DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching",
                                                                                "DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching",Nbin,DRMin,DRMatchingCut);
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetXaxis()->SetTitle("DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching");
  TH1F*  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching  = new TH1F ("DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching",
                                                                                "DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetXaxis()->SetTitle("DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching");

  //   Histogramm after matching selections : No  Matching in event with three jets

  TH1F*  DPt_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching    = new TH1F ("DPt_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching",
                                                                         "DPt_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->GetXaxis()->SetTitle("DPt_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching");
  TH1F*  DEta_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching   = new TH1F ("DEta_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching",
                                                                         "DEta_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->GetXaxis()->SetTitle("DEta_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching");
  TH1F*  DPhi_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching   = new TH1F ("DPhi_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching",
                                                                         "DPhi_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->GetXaxis()->SetTitle("DPhi_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching");
  TH1F*  DR_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching     = new TH1F ("DR_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching",
                                                                         "DR_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching",Nbin,DRMin,DRMax);
  DR_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->GetXaxis()->SetTitle("DR_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching");
  TH1F*  DMass_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching  = new TH1F ("DMass_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching",
                                                                         "DMass_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->GetXaxis()->SetTitle("DMass_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching");
  TH1F*  DPt_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching    = new TH1F ("DPt_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching",
                                                                         "DPt_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching",Nbin,PtEMassMin,PtEMassMax);

  DPt_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->GetXaxis()->SetTitle("DPt_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching");
  TH1F*  DEta_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching   = new TH1F ("DEta_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching",
                                                                         "DEta_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->GetXaxis()->SetTitle("DEta_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching");
  TH1F*  DPhi_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching   = new TH1F ("DPhi_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching",
                                                                         "DPhi_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->GetXaxis()->SetTitle("DPhi_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching");
  TH1F*  DR_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching     = new TH1F ("DR_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching",
                                                                         "DR_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching",Nbin,DRMin,DRMax);
  DR_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->GetXaxis()->SetTitle("DR_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching");
  TH1F*  DMass_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching  = new TH1F ("DMass_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching",
                                                                         "DMass_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->GetXaxis()->SetTitle("DMass_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching");
  TH1F*  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching    = new TH1F ("DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching",
                                                                                "DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching",Nbin,PtEMassMin,PtEMassMax);

  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetXaxis()->SetTitle("DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching");
  TH1F*  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching   = new TH1F ("DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching",
                                                                                "DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetXaxis()->SetTitle("DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching");
  TH1F*  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching   = new TH1F ("DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching",
                                                                                "DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetXaxis()->SetTitle("DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching");
  TH1F*  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching     = new TH1F ("DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching",
                                                                                "DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching",Nbin,DRMin,DRMax);
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetXaxis()->SetTitle("DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching");
  TH1F*  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching  = new TH1F ("DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching",
                                                                                "DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetXaxis()->SetTitle("DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching");

  TH1F*  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching    = new TH1F ("DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching",
                                                                                "DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching",Nbin,PtEMassMin,PtEMassMax);
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetXaxis()->SetTitle("DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching");
  TH1F*  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching   = new TH1F ("DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching",
                                                                                "DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching",Nbin,-EtaMinMax,EtaMinMax);
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetXaxis()->SetTitle("DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching");
  TH1F*  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching   = new TH1F ("DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching",
                                                                                "DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching",Nbin,-PhiMinMax,PhiMinMax);
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetXaxis()->SetTitle("DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching");
  TH1F*  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching     = new TH1F ("DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching",
                                                                                "DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching",Nbin,DRMin,DRMax);
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetXaxis()->SetTitle("DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching");
  TH1F*  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching  = new TH1F ("DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching",
                                                                                "DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching",Nbin,PtEMassMin,PtEMassMax);
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetXaxis()->SetTitle("DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching");


  // Loop on the selected events

  int nEntry = fTree->GetEntries();

  int nExclusive_1Jet_CA8 = 0 ;
  int nExclusive_2Jet_CA8 = 0 ;
  int nExclusive_3Jet_CA8 = 0 ;

  for(int iEntry =0; iEntry<nEntry ; iEntry++){
 
    if( iEntry%1000 == 0 ) std::cout << "reading saved entry " << iEntry << "\r" << std::flush;
 
    fTree->GetEntry(iEntry);

    bool isExclusive_1Jet_CA8 = false ;
    bool isExclusive_2Jet_CA8 = false ;
    bool isExclusive_3Jet_CA8 = false ;
         

    // Basic Selection
    if(!(ggdboostedWevt==1 && W_pt > 200 && GroomedJet_CA8_pt[0] > 200 && event_met_pfmet > 50)) continue ;

    TLorentzVector Parton1    = TLorentzVector(W_Parton_px[0],W_Parton_py[0],W_Parton_pz[0],W_Parton_e[0]);
    TLorentzVector Parton2    = TLorentzVector(W_Parton_px[1],W_Parton_py[1],W_Parton_pz[1],W_Parton_e[1]);

    TLorentzVector GenLepton  = TLorentzVector(W_Lepton_px,W_Lepton_py,W_Lepton_pz,W_Lepton_e);
    TLorentzVector GenNeutrino  = TLorentzVector(W_Met_px,W_Met_py,W_Met_pz,W_Met_e);

    if(Parton1.Pt()<Parton2.Pt()) std::swap(Parton1,Parton2);

    TLorentzVector GenWhad   = Parton1 + Parton2 ;
    TLorentzVector GenWLep   = GenLepton + GenNeutrino ;

    // Distribution of GenWhad and GenWLep
    if(GenWhad.Pt()!=0 && GenWhad.E()!=0 && GenWhad.M()!=0){

      DPt_WLep_GenWhad->Fill(GenWLep.Pt()/GenWhad.Pt());
      DEta_WLep_GenWhad->Fill(GenWLep.Eta()-GenWhad.Eta());
      DE_WLep_GenWhad->Fill(GenWLep.E()/GenWhad.E());
      DMass_WLep_GenWhad->Fill(GenWLep.M()/GenWhad.M());

      DPhi_WLep_GenWhad->Fill(deltaPhi(GenWLep.Phi(),GenWhad.Phi())); 
      DR_WLep_GenWhad->Fill(deltaR(GenWLep.Phi(),GenWhad.Phi(),GenWLep.Eta(),GenWhad.Eta())); 
    }

    // Distribution of 1 Reco Jet and Leptonic Reco W + Matching with hadronic W + Matching Gen hard Jet and hadronic W

    TLorentzVector WLep = TLorentzVector(W_px,W_py,W_pz,W_e);
    
    for( int i =0; i<NJet ; i++){

     TLorentzVector GenGroomedJetCA8    = TLorentzVector(0.,0.,0.,0.);
     TLorentzVector GroomedJetCA8       = TLorentzVector(0.,0.,0.,0.);
     TLorentzVector GenGroomedJetCA8_pr = TLorentzVector(0.,0.,0.,0.);
     TLorentzVector GroomedJetCA8_pr    = TLorentzVector(0.,0.,0.,0.);
     TLorentzVector GenGroomedJetCA8_tr = TLorentzVector(0.,0.,0.,0.);
     TLorentzVector GroomedJetCA8_tr    = TLorentzVector(0.,0.,0.,0.);
     TLorentzVector GenGroomedJetCA8_ft = TLorentzVector(0.,0.,0.,0.);
     TLorentzVector GroomedJetCA8_ft    = TLorentzVector(0.,0.,0.,0.);

     GenGroomedJetCA8.SetPtEtaPhiE(GenGroomedJet_CA8_pt[i],GenGroomedJet_CA8_eta[i],GenGroomedJet_CA8_phi[i],GenGroomedJet_CA8_e[i]);
     GroomedJetCA8.SetPtEtaPhiE(GroomedJet_CA8_pt[i],GroomedJet_CA8_eta[i],GroomedJet_CA8_phi[i],GroomedJet_CA8_e[i]);
     GenGroomedJetCA8_tr.SetPtEtaPhiE(GenGroomedJet_CA8_pt_tr[i],GenGroomedJet_CA8_eta_tr[i],GenGroomedJet_CA8_phi_tr[i],GenGroomedJet_CA8_e_tr[i]);
     GroomedJetCA8_tr.SetPtEtaPhiE(GroomedJet_CA8_pt_tr[i],GroomedJet_CA8_eta_tr[i],GroomedJet_CA8_phi_tr[i],GroomedJet_CA8_e_tr[i]);
     GenGroomedJetCA8_pr.SetPtEtaPhiE(GenGroomedJet_CA8_pt_pr[i],GenGroomedJet_CA8_eta_pr[i],GenGroomedJet_CA8_phi_pr[i],GenGroomedJet_CA8_e_pr[i]);
     GroomedJetCA8_pr.SetPtEtaPhiE(GroomedJet_CA8_pt_pr[i],GroomedJet_CA8_eta_pr[i],GroomedJet_CA8_phi_pr[i],GroomedJet_CA8_e_pr[i]);
     GenGroomedJetCA8_ft.SetPtEtaPhiE(GenGroomedJet_CA8_pt_ft[i],GenGroomedJet_CA8_eta_ft[i],GenGroomedJet_CA8_phi_ft[i],GenGroomedJet_CA8_e_ft[i]);
     GroomedJetCA8_ft.SetPtEtaPhiE(GroomedJet_CA8_pt_ft[i],GroomedJet_CA8_eta_ft[i],GroomedJet_CA8_phi_ft[i],GroomedJet_CA8_e_ft[i]);

     TLorentzVector GenGroomedJetAK5 = TLorentzVector(0.,0.,0.,0.);
     TLorentzVector GroomedJetAK5    = TLorentzVector(0.,0.,0.,0.);

     GenGroomedJetAK5.SetPtEtaPhiE(GenGroomedJet_AK5_pt[i],GenGroomedJet_AK5_eta[i],GenGroomedJet_AK5_phi[i],GenGroomedJet_AK5_e[i]);
     GroomedJetAK5.SetPtEtaPhiE(GroomedJet_AK5_pt[i],GroomedJet_AK5_eta[i],GroomedJet_AK5_phi[i],GroomedJet_AK5_e[i]);
 
     if(GroomedJetCA8.Pt()>PtCut && GroomedJetCA8.M()>0 && GroomedJetCA8.E()>0 && GenWhad.Pt()>0 && GenWhad.E()>0 && GenWhad.M()>0){

       // Comparison with W leptonic 

       DPt_GroomedJet_CA8_WLep[i]    -> Fill(WLep.Pt()/GroomedJetCA8.Pt());
       DEta_GroomedJet_CA8_WLep[i]   -> Fill(WLep.Eta()-GroomedJetCA8.Eta());
       DE_GroomedJet_CA8_WLep[i]     -> Fill(WLep.E()/GroomedJetCA8.E());
       DMass_GroomedJet_CA8_WLep[i]  -> Fill(WLep.M()/GroomedJetCA8.M());

       DPhi_GroomedJet_CA8_WLep[i]   -> Fill(deltaPhi(WLep.Phi(),GroomedJetCA8.Phi())); 
       DR_GroomedJet_CA8_WLep[i]     -> Fill(deltaR(WLep.Eta(),WLep.Phi(),GroomedJetCA8.Eta(),GroomedJetCA8.Phi()));

       DPt_DR_GroomedJet_CA8_WLep[i]    -> Fill(WLep.Pt()/GroomedJetCA8.Pt(),deltaR(WLep.Phi(),GroomedJetCA8.Phi(),WLep.Eta(),GroomedJetCA8.Eta()));
       DMass_DR_GroomedJet_CA8_WLep[i]  -> Fill(WLep.M()/GroomedJetCA8.M(),deltaR(WLep.Phi(),GroomedJetCA8.Phi(),WLep.Eta(),GroomedJetCA8.Eta()));
   
       // Gen jet CA8 distribution

       DPt_GenGroomedJet_CA8_GenWhad[i]   -> Fill(GenGroomedJetCA8.Pt()/GenWhad.Pt());
       DEta_GenGroomedJet_CA8_GenWhad[i]  -> Fill(GenGroomedJetCA8.Eta()-GenWhad.Eta());
       DE_GenGroomedJet_CA8_GenWhad[i]    -> Fill(GenGroomedJetCA8.E()/GenWhad.E());
       DMass_GenGroomedJet_CA8_GenWhad[i] -> Fill(GenGroomedJetCA8.M()/GenWhad.M());

       DPt_GenGroomedJet_CA8_GenWhad_pr[i]   -> Fill(GenGroomedJetCA8_pr.Pt()/GenWhad.Pt());
       DEta_GenGroomedJet_CA8_GenWhad_pr[i]  -> Fill(GenGroomedJetCA8_pr.Eta()-GenWhad.Eta());
       DE_GenGroomedJet_CA8_GenWhad_pr[i]    -> Fill(GenGroomedJetCA8_pr.E()/GenWhad.E());
       DMass_GenGroomedJet_CA8_GenWhad_pr[i] -> Fill(GenGroomedJetCA8_pr.M()/GenWhad.M());

       DPt_GenGroomedJet_CA8_GenWhad_tr[i]   -> Fill(GenGroomedJetCA8_tr.Pt()/GenWhad.Pt());
       DEta_GenGroomedJet_CA8_GenWhad_tr[i]  -> Fill(GenGroomedJetCA8_tr.Eta()-GenWhad.Eta());
       DE_GenGroomedJet_CA8_GenWhad_tr[i]    -> Fill(GenGroomedJetCA8_tr.E()/GenWhad.E());
       DMass_GenGroomedJet_CA8_GenWhad_tr[i] -> Fill(GenGroomedJetCA8_tr.M()/GenWhad.M());

       DPt_GenGroomedJet_CA8_GenWhad_ft[i]   -> Fill(GenGroomedJetCA8_ft.Pt()/GenWhad.Pt());
       DEta_GenGroomedJet_CA8_GenWhad_ft[i]  -> Fill(GenGroomedJetCA8_ft.Eta()-GenWhad.Eta());
       DE_GenGroomedJet_CA8_GenWhad_ft[i]    -> Fill(GenGroomedJetCA8_ft.E()/GenWhad.E());
       DMass_GenGroomedJet_CA8_GenWhad_ft[i] -> Fill(GenGroomedJetCA8_ft.M()/GenWhad.M());


      DPhi_GenGroomedJet_CA8_GenWhad[i] -> Fill(deltaPhi(GenWhad.Phi(),GenGroomedJetCA8.Phi())); 
      DR_GenGroomedJet_CA8_GenWhad[i]   -> Fill(deltaR(GenWhad.Phi(),GenGroomedJetCA8.Phi(),GenWhad.Eta(),GenGroomedJetCA8.Eta()));
      DPt_DR_GenGroomedJet_CA8_GenWhad[i]->Fill(GenGroomedJetCA8.Pt()/GenWhad.Pt(),deltaR(GenWhad.Phi(),GenGroomedJetCA8.Phi(),GenWhad.Eta(),GenGroomedJetCA8.Eta()));
      DMass_DR_GenGroomedJet_CA8_GenWhad[i]->Fill(GenGroomedJetCA8.M()/GenWhad.M(),deltaR(GenWhad.Phi(),GenGroomedJetCA8.Phi(),GenWhad.Eta(),GenGroomedJetCA8.Eta()));
    
      DPhi_GenGroomedJet_CA8_GenWhad_pr[i] -> Fill(deltaPhi(GenWhad.Phi(),GenGroomedJetCA8_pr.Phi()));        
      DR_GenGroomedJet_CA8_GenWhad_pr[i]   -> Fill(deltaR(GenWhad.Phi(),GenGroomedJetCA8_pr.Phi(),GenWhad.Eta(),GenGroomedJetCA8_pr.Eta()));
      DPt_DR_GenGroomedJet_CA8_GenWhad_pr[i]->Fill(GenGroomedJetCA8_pr.Pt()/GenWhad.Pt(),deltaR(GenWhad.Phi(),GenGroomedJetCA8_pr.Phi(),GenWhad.Eta(),GenGroomedJetCA8_pr.Eta()));
      DMass_DR_GenGroomedJet_CA8_GenWhad_pr[i]-> Fill(GenGroomedJetCA8_pr.M()/GenWhad.M(),deltaR(GenWhad.Phi(),GenGroomedJetCA8_pr.Phi(),GenWhad.Eta(),GenGroomedJetCA8_pr.Eta()));

      DPhi_GenGroomedJet_CA8_GenWhad_tr[i] -> Fill(deltaPhi(GenWhad.Phi(),GenGroomedJetCA8_tr.Phi()));        
      DR_GenGroomedJet_CA8_GenWhad_tr[i]   -> Fill(deltaR(GenWhad.Phi(),GenGroomedJetCA8_tr.Phi(),GenWhad.Eta(),GenGroomedJetCA8_tr.Eta()));
      DPt_DR_GenGroomedJet_CA8_GenWhad_tr[i]->Fill(GenGroomedJetCA8_tr.Pt()/GenWhad.Pt(),deltaR(GenWhad.Phi(),GenGroomedJetCA8_tr.Phi(),GenWhad.Eta(),GenGroomedJetCA8_tr.Eta()));
      DMass_DR_GenGroomedJet_CA8_GenWhad_tr[i]-> Fill(GenGroomedJetCA8_tr.M()/GenWhad.M(),deltaR(GenWhad.Phi(),GenGroomedJetCA8_tr.Phi(),GenWhad.Eta(),GenGroomedJetCA8_tr.Eta()));

      DPhi_GenGroomedJet_CA8_GenWhad_ft[i] -> Fill(deltaPhi(GenWhad.Phi(),GenGroomedJetCA8_ft.Phi()));        
      DR_GenGroomedJet_CA8_GenWhad_ft[i]   -> Fill(deltaR(GenWhad.Phi(),GenGroomedJetCA8_ft.Phi(),GenWhad.Eta(),GenGroomedJetCA8_ft.Eta()));
      DPt_DR_GenGroomedJet_CA8_GenWhad_ft[i]->Fill(GenGroomedJetCA8_ft.Pt()/GenWhad.Pt(),deltaR(GenWhad.Phi(),GenGroomedJetCA8_ft.Phi(),GenWhad.Eta(),GenGroomedJetCA8_ft.Eta()));
      DMass_DR_GenGroomedJet_CA8_GenWhad_ft[i]-> Fill(GenGroomedJetCA8_ft.M()/GenWhad.M(),deltaR(GenWhad.Phi(),GenGroomedJetCA8_ft.Phi(),GenWhad.Eta(),GenGroomedJetCA8_ft.Eta()));

      TLorentzVector Gensubjet1 = TLorentzVector(GenGroomedJet_CA8_prsubjet1_px[i],GenGroomedJet_CA8_prsubjet1_py[i],
                                                 GenGroomedJet_CA8_prsubjet1_pz[i],GenGroomedJet_CA8_prsubjet1_e[i]);
      TLorentzVector Gensubjet2 = TLorentzVector(GenGroomedJet_CA8_prsubjet2_px[i],GenGroomedJet_CA8_prsubjet2_py[i],
                                                 GenGroomedJet_CA8_prsubjet2_pz[i],GenGroomedJet_CA8_prsubjet2_e[i]);

      if(Gensubjet1.Pt()<Gensubjet2.Pt()) std::swap(Gensubjet1,Gensubjet2);

      if(Gensubjet1.Pt()>0 && Parton1.Pt()>0){

       DPt_GenGroomedJet_CA8_subjet1_WParton[i]   ->Fill(Gensubjet1.Pt()/Parton1.Pt());
       DEta_GenGroomedJet_CA8_subjet1_WParton[i]  ->Fill(Gensubjet1.Eta()-Parton1.Eta());
       DE_GenGroomedJet_CA8_subjet1_WParton[i]    ->Fill(Gensubjet1.E()/Parton1.E());
       DMass_GenGroomedJet_CA8_subjet1_WParton[i] ->Fill(Gensubjet1.M()/Parton1.M());

       DPhi_GenGroomedJet_CA8_subjet1_WParton[i] -> Fill(deltaPhi(Gensubjet1.Phi(),Parton1.Phi())); 
       DR_GenGroomedJet_CA8_subjet1_WParton[i]   -> Fill(deltaR(Gensubjet1.Phi(),Parton1.Phi(),Gensubjet1.Eta(),Parton1.Eta()));
       DPt_DR_GenGroomedJet_CA8_subjet1_WParton[i]->Fill(Gensubjet1.Pt()/Parton1.Pt(),deltaR(Gensubjet1.Phi(),Parton1.Phi(),Gensubjet1.Eta(),Parton1.Eta()));
       DMass_DR_GenGroomedJet_CA8_subjet1_WParton[i]-> Fill(Gensubjet1.M()/Parton1.M(),deltaR(Gensubjet1.Phi(),Parton1.Phi(),Gensubjet1.Eta(),Parton1.Eta()));
    
      }

      if(Gensubjet2.Pt()>0 && Parton2.Pt()>0){
       DPt_GenGroomedJet_CA8_subjet2_WParton[i]   ->Fill(Gensubjet2.Pt()/Parton2.Pt());
       DEta_GenGroomedJet_CA8_subjet2_WParton[i]  ->Fill(Gensubjet2.Eta()-Parton2.Eta());
       DE_GenGroomedJet_CA8_subjet2_WParton[i]    ->Fill(Gensubjet2.E()/Parton2.E());
       DMass_GenGroomedJet_CA8_subjet2_WParton[i] ->Fill(Gensubjet2.M()/Parton2.M());
      
       DPhi_GenGroomedJet_CA8_subjet2_WParton[i]    -> Fill(deltaPhi(Gensubjet1.Phi(),Parton1.Phi())); 
       DR_GenGroomedJet_CA8_subjet2_WParton[i]      -> Fill(deltaR(Gensubjet1.Phi(),Parton1.Phi(),Gensubjet1.Eta(),Parton1.Eta()));
       DPt_DR_GenGroomedJet_CA8_subjet2_WParton[i]  -> Fill(Gensubjet2.Pt()/Parton2.Pt(),deltaR(Gensubjet2.Phi(),Parton2.Phi(),Gensubjet2.Eta(),Parton2.Eta()));
       DMass_DR_GenGroomedJet_CA8_subjet2_WParton[i]-> Fill(Gensubjet2.M()/Parton2.M(),deltaR(Gensubjet2.Phi(),Parton2.Phi(),Gensubjet2.Eta(),Parton2.Eta()));

      }
        
      // Fill Reco CA8 Info

      DPt_GroomedJet_CA8_GenWhad[i]   -> Fill(GroomedJetCA8.Pt()/GenWhad.Pt());
      DEta_GroomedJet_CA8_GenWhad[i]  -> Fill(GroomedJetCA8.Eta()-GenWhad.Eta());
      DE_GroomedJet_CA8_GenWhad[i]    -> Fill(GroomedJetCA8.E()/GenWhad.E());
      DMass_GroomedJet_CA8_GenWhad[i] -> Fill(GroomedJetCA8.M()/GenWhad.M());

      DPt_GroomedJet_CA8_GenWhad_pr[i]   -> Fill(GroomedJetCA8_pr.Pt()/GenWhad.Pt());
      DEta_GroomedJet_CA8_GenWhad_pr[i]  -> Fill(GroomedJetCA8_pr.Eta()-GenWhad.Eta());
      DE_GroomedJet_CA8_GenWhad_pr[i]    -> Fill(GroomedJetCA8_pr.E()/GenWhad.E());
      DMass_GroomedJet_CA8_GenWhad_pr[i] -> Fill(GroomedJetCA8_pr.M()/GenWhad.M());

      DPt_GroomedJet_CA8_GenWhad_ft[i]   -> Fill(GroomedJetCA8_ft.Pt()/GenWhad.Pt());
      DEta_GroomedJet_CA8_GenWhad_ft[i]  -> Fill(GroomedJetCA8_ft.Eta()-GenWhad.Eta());
      DE_GroomedJet_CA8_GenWhad_ft[i]    -> Fill(GroomedJetCA8_ft.E()/GenWhad.E());
      DMass_GroomedJet_CA8_GenWhad_ft[i] -> Fill(GroomedJetCA8_ft.M()/GenWhad.M());

      DPt_GroomedJet_CA8_GenWhad_tr[i]   -> Fill(GroomedJetCA8_tr.Pt()/GenWhad.Pt());
      DEta_GroomedJet_CA8_GenWhad_tr[i]  -> Fill(GroomedJetCA8_tr.Eta()-GenWhad.Eta());
      DE_GroomedJet_CA8_GenWhad_tr[i]    -> Fill(GroomedJetCA8_tr.E()/GenWhad.E());
      DMass_GroomedJet_CA8_GenWhad_tr[i] -> Fill(GroomedJetCA8_tr.M()/GenWhad.M());

      DPhi_GroomedJet_CA8_GenWhad[i] -> Fill(deltaPhi(GenWhad.Phi(),GroomedJetCA8.Phi())); 
      DR_GroomedJet_CA8_GenWhad[i]   -> Fill(deltaR(GenWhad.Phi(),GroomedJetCA8.Phi(),GenWhad.Eta(),GroomedJetCA8.Eta()));   

      DPt_DR_GroomedJet_CA8_GenWhad[i] -> Fill(GroomedJetCA8.Pt()/GenWhad.Pt(),deltaR(GenWhad.Phi(),GroomedJetCA8.Phi(),GenWhad.Eta(),GroomedJetCA8.Eta())); 
      DMass_DR_GroomedJet_CA8_GenWhad[i] -> Fill(GroomedJetCA8.M()/GenWhad.M(),deltaR(GenWhad.Phi(),GroomedJetCA8.Phi(),GenWhad.Eta(),GroomedJetCA8.Eta()));   


      DPhi_GroomedJet_CA8_GenWhad_pr[i] -> Fill(deltaPhi(GenWhad.Phi(),GroomedJetCA8_pr.Phi())); 
      DR_GroomedJet_CA8_GenWhad_pr[i]   -> Fill(deltaR(GenWhad.Phi(),GroomedJetCA8_pr.Phi(),GenWhad.Eta(),GroomedJetCA8_pr.Eta()));      

      DPt_DR_GroomedJet_CA8_GenWhad_pr[i] -> Fill(GroomedJetCA8_pr.Pt()/GenWhad.Pt(),deltaR(GenWhad.Phi(),GroomedJetCA8_pr.Phi(),GenWhad.Eta(),GroomedJetCA8_pr.Eta())); 
      DMass_DR_GroomedJet_CA8_GenWhad_pr[i] -> Fill(GroomedJetCA8_pr.M()/GenWhad.M(),deltaR(GenWhad.Phi(),GroomedJetCA8_pr.Phi(),GenWhad.Eta(),GroomedJetCA8_pr.Eta()));   

      DPhi_GroomedJet_CA8_GenWhad_tr[i] -> Fill(deltaPhi(GenWhad.Phi(),GroomedJetCA8_tr.Phi())); 
      DR_GroomedJet_CA8_GenWhad_tr[i]   -> Fill(deltaR(GenWhad.Phi(),GroomedJetCA8_tr.Phi(),GenWhad.Eta(),GroomedJetCA8_tr.Eta()));

      DPt_DR_GroomedJet_CA8_GenWhad_tr[i] -> Fill(GroomedJetCA8_tr.Pt()/GenWhad.Pt(),deltaR(GenWhad.Phi(),GroomedJetCA8_tr.Phi(),GenWhad.Eta(),GroomedJetCA8_tr.Eta())); 
      DMass_DR_GroomedJet_CA8_GenWhad_tr[i] -> Fill(GroomedJetCA8_tr.M()/GenWhad.M(),deltaR(GenWhad.Phi(),GroomedJetCA8_tr.Phi(),GenWhad.Eta(),GroomedJetCA8_tr.Eta()));   

      DPhi_GroomedJet_CA8_GenWhad_ft[i] -> Fill(deltaPhi(GenWhad.Phi(),GroomedJetCA8_ft.Phi())); 
      DR_GroomedJet_CA8_GenWhad_ft[i]   -> Fill(deltaR(GenWhad.Phi(),GroomedJetCA8_ft.Phi(),GenWhad.Eta(),GroomedJetCA8_ft.Eta()));

      DPt_DR_GroomedJet_CA8_GenWhad_ft[i] -> Fill(GroomedJetCA8_ft.Pt()/GenWhad.Pt(),deltaR(GenWhad.Phi(),GroomedJetCA8_ft.Phi(),GenWhad.Eta(),GroomedJetCA8_ft.Eta())); 
      DMass_DR_GroomedJet_CA8_GenWhad_ft[i] -> Fill(GroomedJetCA8_ft.M()/GenWhad.M(),deltaR(GenWhad.Phi(),GroomedJetCA8_ft.Phi(),GenWhad.Eta(),GroomedJetCA8_ft.Eta()));   


      TLorentzVector subjet1 = TLorentzVector(GroomedJet_CA8_prsubjet1_px[i],GroomedJet_CA8_prsubjet1_py[i],
                                             GroomedJet_CA8_prsubjet1_pz[i],GroomedJet_CA8_prsubjet1_e[i]);
      TLorentzVector subjet2 = TLorentzVector(GroomedJet_CA8_prsubjet2_px[i],GroomedJet_CA8_prsubjet2_py[i],
                                             GroomedJet_CA8_prsubjet2_pz[i],GroomedJet_CA8_prsubjet2_e[i]);


      if(subjet1.Pt()<subjet2.Pt()) std::swap(subjet1,subjet2);

      if(subjet1.Pt()!=0 && Parton1.Pt()!=0){
       DPt_GroomedJet_CA8_subjet1_WParton[i]   ->Fill(subjet1.Pt()/Parton1.Pt());
       DEta_GroomedJet_CA8_subjet1_WParton[i]  ->Fill(subjet1.Eta()-Parton1.Eta());
       DE_GroomedJet_CA8_subjet1_WParton[i]    ->Fill(subjet1.E()/Parton1.E());
       DMass_GroomedJet_CA8_subjet1_WParton[i] ->Fill(subjet1.M()/Parton1.M());

       DPhi_GroomedJet_CA8_subjet1_WParton[i] -> Fill(deltaPhi(subjet1.Phi(),Parton1.Phi())); 
       DR_GroomedJet_CA8_subjet1_WParton[i]   -> Fill(deltaR(subjet1.Phi(),Parton1.Phi(),subjet1.Eta(),Parton1.Eta()));

       DPt_DR_GroomedJet_CA8_subjet1_WParton[i] -> Fill(subjet1.Pt()/Parton1.Pt(),deltaR(subjet1.Phi(),Parton1.Phi(),subjet1.Eta(),Parton1.Eta())); 
       DMass_DR_GroomedJet_CA8_subjet1_WParton[i] -> Fill(subjet1.M()/GroomedJetCA8.M(),deltaR(subjet1.Phi(),Parton1.Phi(),subjet1.Eta(),Parton1.Eta()));   

      }

      if(subjet2.Pt()!=0 && Parton2.Pt()!=0){
       DPt_GroomedJet_CA8_subjet2_WParton[i]   ->Fill(subjet2.Pt()/Parton2.Pt());
       DEta_GroomedJet_CA8_subjet2_WParton[i]  ->Fill(subjet2.Eta()-Parton2.Eta());
       DE_GroomedJet_CA8_subjet2_WParton[i]    ->Fill(subjet2.E()/Parton2.E());
       DMass_GroomedJet_CA8_subjet2_WParton[i] ->Fill(subjet2.M()/Parton2.M());


       DPhi_GroomedJet_CA8_subjet2_WParton[i] -> Fill(deltaPhi(subjet2.Phi(),Parton2.Phi())); 
       DR_GroomedJet_CA8_subjet2_WParton[i]   -> Fill(deltaR(subjet2.Phi(),Parton2.Phi(),subjet2.Eta(),Parton2.Eta()));

       DPt_DR_GroomedJet_CA8_subjet2_WParton[i] -> Fill(subjet2.Pt()/Parton1.Pt(),deltaR(subjet2.Phi(),Parton2.Phi(),subjet2.Eta(),Parton2.Eta())); 
       DMass_DR_GroomedJet_CA8_subjet2_WParton[i] -> Fill(subjet2.M()/GroomedJetCA8.M(),deltaR(subjet2.Phi(),Parton2.Phi(),subjet2.Eta(),Parton2.Eta()));   

      }
   
     }

     if(GroomedJetAK5.Pt()>PtCut && GroomedJetAK5.M()>=0 && GenWhad.Pt()>0 && GenWhad.M()>0 && GenGroomedJetAK5.Pt()>=0 && GenGroomedJetAK5.M()>=0 ){

      DPt_GenGroomedJet_AK5_GenWhad[i]  -> Fill(GenGroomedJetAK5.Pt()/GenWhad.Pt());
      DEta_GenGroomedJet_AK5_GenWhad[i] -> Fill(GenGroomedJetAK5.Eta()-GenWhad.Eta());
      DE_GenGroomedJet_AK5_GenWhad[i]  -> Fill(GenGroomedJetAK5.E()/GenWhad.E());
      DMass_GenGroomedJet_AK5_GenWhad[i] -> Fill(GenGroomedJetAK5.M()/GenWhad.M());

      DPhi_GenGroomedJet_AK5_GenWhad[i] -> Fill(deltaPhi(GenWhad.Phi(),GenGroomedJetAK5.Phi())); 
      DR_GenGroomedJet_AK5_GenWhad[i]   -> Fill(deltaR(GenWhad.Phi(),GenGroomedJetAK5.Phi(),GenWhad.Eta(),GenGroomedJetAK5.Eta()));

      DPt_DR_GenGroomedJet_AK5_GenWhad[i] -> Fill(GenGroomedJetAK5.Pt()/GenWhad.Pt(),deltaR(GenWhad.Phi(),GenGroomedJetAK5.Phi(),GenWhad.Eta(),GenGroomedJetAK5.Eta()));
      DMass_DR_GenGroomedJet_AK5_GenWhad[i] -> Fill(GenGroomedJetAK5.M()/GenWhad.M(),deltaR(GenWhad.Phi(),GenGroomedJetAK5.Phi(),GenWhad.Eta(),GenGroomedJetAK5.Eta()));
   

      DPt_GroomedJet_AK5_GenWhad[i]  -> Fill(GroomedJetAK5.Pt()/GenWhad.Pt());
      DEta_GroomedJet_AK5_GenWhad[i] -> Fill(GroomedJetAK5.Eta()-GenWhad.Eta());
      DE_GroomedJet_AK5_GenWhad[i]  -> Fill(GroomedJetAK5.E()/GenWhad.E());
      DMass_GroomedJet_AK5_GenWhad[i] -> Fill(GroomedJetAK5.M()/GenWhad.M());

      DPhi_GroomedJet_AK5_GenWhad[i] -> Fill(deltaPhi(GenWhad.Phi(),GroomedJetAK5.Phi())); 
      DR_GroomedJet_AK5_GenWhad[i]   -> Fill(deltaR(GenWhad.Phi(),GroomedJetAK5.Phi(),GenWhad.Eta(),GroomedJetAK5.Eta())); 

      DPt_DR_GroomedJet_AK5_GenWhad[i] -> Fill(GroomedJetAK5.Pt()/GenWhad.Pt(),deltaR(GenWhad.Phi(),GroomedJetAK5.Phi(),GenWhad.Eta(),GroomedJetAK5.Eta()));
      DMass_DR_GroomedJet_AK5_GenWhad[i] -> Fill(GroomedJetAK5.M()/GenWhad.M(),deltaR(GenWhad.Phi(),GroomedJetAK5.Phi(),GenWhad.Eta(),GroomedJetAK5.Eta()));
    }

     if(GroomedJetCA8.Pt()>PtCut && GroomedJetCA8.M()>0 && GroomedJetCA8.E()>0 && GenWhad.Pt()>0 && GenWhad.E()>0 && GenWhad.M()>0){

      // One CA8 Jet exclusive Category

      if( i==0 && isExclusive_1Jet_CA8 == false ) isExclusive_1Jet_CA8 = true ;
      if( i!=0 )  isExclusive_1Jet_CA8 = false;

      // Two CA8 Jet exclusive Category

     if( i==1 && isExclusive_2Jet_CA8 == false) isExclusive_2Jet_CA8 = true ;
     if( i>1) isExclusive_2Jet_CA8 = false ;

      // Three CA8 Jet exclusive Category

     if( i==2 && isExclusive_3Jet_CA8 == false) isExclusive_3Jet_CA8 = true ; 
     if( i>2) isExclusive_3Jet_CA8 = false ; 
     
     if(deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta())<DRMatchingCut){

       DPt_GroomedCA8_GenWhad_DRCut[i]->Fill(GroomedJetCA8.Pt()/GenWhad.Pt());
       DEta_GroomedCA8_GenWhad_DRCut[i]->Fill(GroomedJetCA8.Eta()-GenWhad.Eta());
       DPhi_GroomedCA8_GenWhad_DRCut[i]->Fill(deltaPhi(GroomedJetCA8.Phi(),GenWhad.Phi()));
       DR_GroomedCA8_GenWhad_DRCut[i]->Fill(deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta()));
       DMass_GroomedCA8_GenWhad_DRCut[i]->Fill(GroomedJetCA8.M()/GenWhad.M());

       if(fabs(GroomedJetCA8.Pt()/GenWhad.Pt()-1.)<PtMatchingFraction){

       DPt_GroomedCA8_GenWhad_DRCut_DPtCut[i]->Fill(GroomedJetCA8.Pt()/GenWhad.Pt());
       DEta_GroomedCA8_GenWhad_DRCut_DPtCut[i]->Fill(GroomedJetCA8.Eta()-GenWhad.Eta());
       DPhi_GroomedCA8_GenWhad_DRCut_DPtCut[i]->Fill(deltaPhi(GroomedJetCA8.Phi(),GenWhad.Phi()));
       DR_GroomedCA8_GenWhad_DRCut_DPtCut[i]->Fill(deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta()));
       DMass_GroomedCA8_GenWhad_DRCut_DPtCut[i]->Fill(GroomedJetCA8.M()/GenWhad.M());

      }
     }

     if(deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta())<DRMatchingCut){

       DPt_GroomedAK5_GenWhad_DRCut[i]->Fill(GroomedJetAK5.Pt()/GenWhad.Pt());
       DEta_GroomedAK5_GenWhad_DRCut[i]->Fill(GroomedJetAK5.Eta()-GenWhad.Eta());
       DPhi_GroomedAK5_GenWhad_DRCut[i]->Fill(deltaPhi(GroomedJetAK5.Phi(),GenWhad.Phi()));
       DR_GroomedAK5_GenWhad_DRCut[i]->Fill(deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta()));
       DMass_GroomedAK5_GenWhad_DRCut[i]->Fill(GroomedJetAK5.M()/GenWhad.M());

       if(fabs(GroomedJetAK5.Pt()/GenWhad.Pt()-1.)<PtMatchingFraction){

        DPt_GroomedAK5_GenWhad_DRCut_DPtCut[i]->Fill(GroomedJetAK5.Pt()/GenWhad.Pt());
        DEta_GroomedAK5_GenWhad_DRCut_DPtCut[i]->Fill(GroomedJetAK5.Eta()-GenWhad.Eta());
        DPhi_GroomedAK5_GenWhad_DRCut_DPtCut[i]->Fill(deltaPhi(GroomedJetAK5.Phi(),GenWhad.Phi()));
        DR_GroomedAK5_GenWhad_DRCut_DPtCut[i]->Fill(deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta()));
        DMass_GroomedAK5_GenWhad_DRCut_DPtCut[i]->Fill(GroomedJetAK5.M()/GenWhad.M());
       }
     }
    }
   }

   TLorentzVector GroomedJetCA8 = TLorentzVector(0.,0.,0.,0.);
   TLorentzVector GroomedJetAK5 = TLorentzVector(0.,0.,0.,0.);
   GroomedJetCA8.SetPtEtaPhiE(GroomedJet_CA8_pt[0],GroomedJet_CA8_eta[0],GroomedJet_CA8_phi[0],GroomedJet_CA8_e[0]);
   GroomedJetAK5.SetPtEtaPhiE(GroomedJet_AK5_pt[0],GroomedJet_AK5_eta[0],GroomedJet_AK5_phi[0],GroomedJet_AK5_e[0]);
   
   if(isExclusive_1Jet_CA8){

     nExclusive_1Jet_CA8++;
   
     if(deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta())<DRMatchingCut){
     
      DPt_GroomedCA8_GenWhad_DRCut_1Jet  ->Fill(GroomedJetCA8.Pt()/GenWhad.Pt());
      DEta_GroomedCA8_GenWhad_DRCut_1Jet ->Fill(GroomedJetCA8.Eta()-GenWhad.Eta());
      DPhi_GroomedCA8_GenWhad_DRCut_1Jet ->Fill(deltaPhi(GroomedJetCA8.Phi(),GenWhad.Phi()));
      DR_GroomedCA8_GenWhad_DRCut_1Jet   ->Fill(deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta()));
      DMass_GroomedCA8_GenWhad_DRCut_1Jet->Fill(GroomedJetCA8.M()/GenWhad.M());

      if(fabs(GroomedJetCA8.Pt()/GenWhad.Pt()-1.)<PtMatchingFraction){

       DPt_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet ->Fill(GroomedJetCA8.Pt()/GenWhad.Pt());
       DEta_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet ->Fill(GroomedJetCA8.Eta()-GenWhad.Eta());
       DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet ->Fill(deltaPhi(GroomedJetCA8.Phi(),GenWhad.Phi()));
       DR_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet   ->Fill(deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta()));
       DMass_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet->Fill(GroomedJetCA8.M()/GenWhad.M());   
      }
     }

     if(deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta())<DRMatchingCut){

       DPt_GroomedAK5_GenWhad_DRCut_1Jet  ->Fill(GroomedJetAK5.Pt()/GenWhad.Pt());
       DEta_GroomedAK5_GenWhad_DRCut_1Jet ->Fill(GroomedJetAK5.Eta()-GenWhad.Eta());
       DPhi_GroomedAK5_GenWhad_DRCut_1Jet ->Fill(deltaPhi(GroomedJetAK5.Phi(),GenWhad.Phi()));
       DR_GroomedAK5_GenWhad_DRCut_1Jet   ->Fill(deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta()));
       DMass_GroomedAK5_GenWhad_DRCut_1Jet->Fill(GroomedJetAK5.M()/GenWhad.M());

       if(fabs(GroomedJetAK5.Pt()/GenWhad.Pt()-1.)<PtMatchingFraction){

        DPt_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet  ->Fill(GroomedJetAK5.Pt()/GenWhad.Pt());
        DEta_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet ->Fill(GroomedJetAK5.Eta()-GenWhad.Eta());
        DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet ->Fill(deltaPhi(GroomedJetAK5.Phi(),GenWhad.Phi()));
        DR_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet   ->Fill(deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta()));
        DMass_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->Fill(GroomedJetAK5.M()/GenWhad.M());
       }
      }
    }

   if(isExclusive_2Jet_CA8){

     nExclusive_2Jet_CA8++;
     TLorentzVector GroomedJetCA8_1 = TLorentzVector(0.,0.,0.,0.);
     TLorentzVector GroomedJetAK5_1 = TLorentzVector(0.,0.,0.,0.);
     GroomedJetCA8_1.SetPtEtaPhiE(GroomedJet_CA8_pt[1],GroomedJet_CA8_eta[1],GroomedJet_CA8_phi[1],GroomedJet_CA8_e[1]);
     GroomedJetAK5_1.SetPtEtaPhiE(GroomedJet_AK5_pt[1],GroomedJet_AK5_eta[1],GroomedJet_AK5_phi[1],GroomedJet_AK5_e[1]);

     // matching only for the hard (DR)

     if(deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta())<=DRMatchingCut &&
        deltaR(GroomedJetCA8_1.Phi(),GenWhad.Phi(),GroomedJetCA8_1.Eta(),GenWhad.Eta())>DRMatchingCut){
     
     DPt_GroomedCA8_GenWhad_DRCut_2Jet  ->Fill(GroomedJetCA8.Pt()/GenWhad.Pt());
     DEta_GroomedCA8_GenWhad_DRCut_2Jet ->Fill(GroomedJetCA8.Eta()-GenWhad.Eta());
     DPhi_GroomedCA8_GenWhad_DRCut_2Jet ->Fill(deltaPhi(GroomedJetCA8.Phi(),GenWhad.Phi()));
     DR_GroomedCA8_GenWhad_DRCut_2Jet   ->Fill(deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta()));
     DMass_GroomedCA8_GenWhad_DRCut_2Jet->Fill(GroomedJetCA8.M()/GenWhad.M());
   
     }

     // Double matching only for the hard (DR)

     if( deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta())<=DRMatchingCut && 
         deltaR(GroomedJetCA8_1.Phi(),GenWhad.Phi(),GroomedJetCA8_1.Eta(),GenWhad.Eta())<=DRMatchingCut){

     DPt_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching  ->Fill(GroomedJetCA8_1.Pt()/GenWhad.Pt());
     DEta_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching ->Fill(GroomedJetCA8_1.Eta()-GenWhad.Eta());
     DPhi_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching ->Fill(deltaPhi(GroomedJetCA8_1.Phi(),GenWhad.Phi()));
     DR_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching   ->Fill(deltaR(GroomedJetCA8_1.Phi(),GenWhad.Phi(),GroomedJetCA8_1.Eta(),GenWhad.Eta()));
     DMass_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->Fill(GroomedJetCA8_1.M()/GenWhad.M());
     }

    // Inverted Matching (DR)

    if(deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta())>DRMatchingCut && 
       deltaR(GroomedJetCA8_1.Phi(),GenWhad.Phi(),GroomedJetCA8_1.Eta(),GenWhad.Eta())<=DRMatchingCut){

      DPt_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching  ->Fill(GroomedJetCA8_1.Pt()/GenWhad.Pt());
      DEta_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching ->Fill(GroomedJetCA8_1.Eta()-GenWhad.Eta());
      DPhi_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching ->Fill(deltaPhi(GroomedJetCA8_1.Phi(),GenWhad.Phi()));
      DR_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching   ->Fill(deltaR(GroomedJetCA8_1.Phi(),GenWhad.Phi(),GroomedJetCA8_1.Eta(),GenWhad.Eta()));
      DMass_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->Fill(GroomedJetCA8_1.M()/GenWhad.M());
    }

    // No matching (DR)

    if(deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta())>DRMatchingCut && 
       deltaR(GroomedJetCA8_1.Phi(),GenWhad.Phi(),GroomedJetCA8_1.Eta(),GenWhad.Eta())>DRMatchingCut){

      DPt_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching  ->Fill(GroomedJetCA8.Pt()/GenWhad.Pt());
      DEta_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching ->Fill(GroomedJetCA8.Eta()-GenWhad.Eta());
      DPhi_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching ->Fill(deltaPhi(GroomedJetCA8.Phi(),GenWhad.Phi()));
      DR_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching   ->Fill(deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta()));
      DMass_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->Fill(GroomedJetCA8.M()/GenWhad.M());
    }
  
    // Matching hard Pt (DR,DPt)   
 
    if((deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta())<=DRMatchingCut && fabs(GroomedJetCA8.Pt()/GenWhad.Pt()-1.)<=PtMatchingFraction) && 
        (deltaR(GroomedJetCA8_1.Phi(),GenWhad.Phi(),GroomedJetCA8_1.Eta(),GenWhad.Eta())>DRMatchingCut || fabs(GroomedJetCA8_1.Pt()/GenWhad.Pt()-1.)>PtMatchingFraction)){
 
        DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet ->Fill(GroomedJetCA8.Pt()/GenWhad.Pt());
        DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet ->Fill(GroomedJetCA8.Eta()-GenWhad.Eta());
        DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet ->Fill(deltaPhi(GroomedJetCA8.Phi(),GenWhad.Phi()));
        DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet   ->Fill(deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta()));
        DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->Fill(GroomedJetCA8.M()/GenWhad.M());   
     }

    // Doubke Matching (DR,DPt)
   
    if((deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta())<=DRMatchingCut && fabs(GroomedJetCA8.Pt()/GenWhad.Pt()-1.)<=PtMatchingFraction) &&
        (deltaR(GroomedJetCA8_1.Phi(),GenWhad.Phi(),GroomedJetCA8_1.Eta(),GenWhad.Eta())<=DRMatchingCut && fabs(GroomedJetCA8_1.Pt()/GenWhad.Pt()-1.)<=PtMatchingFraction)){

        DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching   ->Fill(GroomedJetCA8_1.Pt()/GenWhad.Pt());
        DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching  ->Fill(GroomedJetCA8_1.Eta()-GenWhad.Eta());
        DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching  ->Fill(deltaPhi(GroomedJetCA8_1.Phi(),GenWhad.Phi()));
        DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching    ->Fill(deltaR(GroomedJetCA8_1.Phi(),GenWhad.Phi(),GroomedJetCA8_1.Eta(),GenWhad.Eta()));
        DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching ->Fill(GroomedJetCA8_1.M()/GenWhad.M());   
    }
    
    // Inverted Matching (DR,DPt)
  
    if( (deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta())>DRMatchingCut || fabs(GroomedJetCA8.Pt()/GenWhad.Pt()-1.)>PtMatchingFraction) && 
        (deltaR(GroomedJetCA8_1.Phi(),GenWhad.Phi(),GroomedJetCA8_1.Eta(),GenWhad.Eta())<=DRMatchingCut && fabs(GroomedJetCA8_1.Pt()/GenWhad.Pt()-1.)<=PtMatchingFraction)){

        DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching   ->Fill(GroomedJetCA8_1.Pt()/GenWhad.Pt());
        DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching  ->Fill(GroomedJetCA8_1.Eta()-GenWhad.Eta());
        DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching  ->Fill(deltaPhi(GroomedJetCA8_1.Phi(),GenWhad.Phi()));
        DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching    ->Fill(deltaR(GroomedJetCA8_1.Phi(),GenWhad.Phi(),GroomedJetCA8_1.Eta(),GenWhad.Eta()));
        DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching ->Fill(GroomedJetCA8_1.M()/GenWhad.M());   
     }
    
    // No Matching (DR,DPt)

    if((deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta())>DRMatchingCut || fabs(GroomedJetCA8.Pt()/GenWhad.Pt()-1.)>PtMatchingFraction  ) &&
       (deltaR(GroomedJetCA8_1.Phi(),GenWhad.Phi(),GroomedJetCA8_1.Eta(),GenWhad.Eta())>DRMatchingCut || fabs(GroomedJetCA8_1.Pt()/GenWhad.Pt()-1.)>PtMatchingFraction )){

        DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching   ->Fill(GroomedJetCA8.Pt()/GenWhad.Pt());
        DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching  ->Fill(GroomedJetCA8.Eta()-GenWhad.Eta());
        DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching  ->Fill(deltaPhi(GroomedJetCA8.Phi(),GenWhad.Phi()));
        DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching    ->Fill(deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta()));
        DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching ->Fill(GroomedJetCA8.M()/GenWhad.M());   
     }


    // Same Category for AK5

     if(deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta())<=DRMatchingCut &&
        deltaR(GroomedJetAK5_1.Phi(),GenWhad.Phi(),GroomedJetAK5_1.Eta(),GenWhad.Eta())>DRMatchingCut){

       DPt_GroomedAK5_GenWhad_DRCut_2Jet  ->Fill(GroomedJetAK5.Pt()/GenWhad.Pt());
       DEta_GroomedAK5_GenWhad_DRCut_2Jet ->Fill(GroomedJetAK5.Eta()-GenWhad.Eta());
       DPhi_GroomedAK5_GenWhad_DRCut_2Jet ->Fill(deltaPhi(GroomedJetAK5.Phi(),GenWhad.Phi()));
       DR_GroomedAK5_GenWhad_DRCut_2Jet   ->Fill(deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta()));
       DMass_GroomedAK5_GenWhad_DRCut_2Jet->Fill(GroomedJetAK5.M()/GenWhad.M());
     }

     if(deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta())<=DRMatchingCut  && 
        deltaR(GroomedJetAK5_1.Phi(),GenWhad.Phi(),GroomedJetAK5_1.Eta(),GenWhad.Eta())<=DRMatchingCut){

     DPt_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching  ->Fill(GroomedJetAK5_1.Pt()/GenWhad.Pt());
     DEta_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching ->Fill(GroomedJetAK5_1.Eta()-GenWhad.Eta());
     DPhi_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching ->Fill(deltaPhi(GroomedJetAK5_1.Phi(),GenWhad.Phi()));
     DR_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching   ->Fill(deltaR(GroomedJetAK5_1.Phi(),GenWhad.Phi(),GroomedJetAK5_1.Eta(),GenWhad.Eta()));
     DMass_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->Fill(GroomedJetAK5_1.M()/GenWhad.M());

     }


    if(deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta())>DRMatchingCut && 
       deltaR(GroomedJetAK5_1.Phi(),GenWhad.Phi(),GroomedJetAK5_1.Eta(),GenWhad.Eta())<=DRMatchingCut){

      DPt_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching  ->Fill(GroomedJetAK5_1.Pt()/GenWhad.Pt());
      DEta_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching ->Fill(GroomedJetAK5_1.Eta()-GenWhad.Eta());
      DPhi_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching ->Fill(deltaPhi(GroomedJetAK5_1.Phi(),GenWhad.Phi()));
      DR_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching   ->Fill(deltaR(GroomedJetAK5_1.Phi(),GenWhad.Phi(),GroomedJetAK5_1.Eta(),GenWhad.Eta()));
      DMass_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->Fill(GroomedJetAK5_1.M()/GenWhad.M());
    }

    if(deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta())>DRMatchingCut &&
       deltaR(GroomedJetAK5_1.Phi(),GenWhad.Phi(),GroomedJetAK5_1.Eta(),GenWhad.Eta())>DRMatchingCut){

      DPt_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching  ->Fill(GroomedJetAK5.Pt()/GenWhad.Pt());
      DEta_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching ->Fill(GroomedJetAK5.Eta()-GenWhad.Eta());
      DPhi_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching ->Fill(deltaPhi(GroomedJetAK5.Phi(),GenWhad.Phi()));
      DR_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching   ->Fill(deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta()));
      DMass_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->Fill(GroomedJetAK5.M()/GenWhad.M());

    }

     if((deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta())<=DRMatchingCut  && fabs(GroomedJetAK5.Pt()/GenWhad.Pt()-1.)<=PtMatchingFraction) && 
        (deltaR(GroomedJetAK5_1.Phi(),GenWhad.Phi(),GroomedJetAK5_1.Eta(),GenWhad.Eta())>DRMatchingCut || fabs(GroomedJetAK5_1.Pt()/GenWhad.Pt()-1.)>PtMatchingFraction)){

        DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet ->Fill(GroomedJetAK5.Pt()/GenWhad.Pt());
        DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet ->Fill(GroomedJetAK5.Eta()-GenWhad.Eta());
        DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet ->Fill(deltaPhi(GroomedJetAK5.Phi(),GenWhad.Phi()));
        DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet   ->Fill(deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta()));
        DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->Fill(GroomedJetAK5.M()/GenWhad.M());   
     }

     if((deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta())<=DRMatchingCut  && fabs(GroomedJetAK5.Pt()/GenWhad.Pt()-1.)<=PtMatchingFraction) &&
	(deltaR(GroomedJetAK5_1.Phi(),GenWhad.Phi(),GroomedJetAK5_1.Eta(),GenWhad.Eta())<=DRMatchingCut && fabs(GroomedJetAK5_1.Pt()/GenWhad.Pt()-1.)<=PtMatchingFraction)){

        DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching   ->Fill(GroomedJetAK5_1.Pt()/GenWhad.Pt());
        DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching  ->Fill(GroomedJetAK5_1.Eta()-GenWhad.Eta());
        DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching  ->Fill(deltaPhi(GroomedJetAK5_1.Phi(),GenWhad.Phi()));
        DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching    ->Fill(deltaR(GroomedJetAK5_1.Phi(),GenWhad.Phi(),GroomedJetAK5_1.Eta(),GenWhad.Eta()));
        DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching ->Fill(GroomedJetAK5_1.M()/GenWhad.M());   

       }
    
    if( (deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta())>DRMatchingCut || fabs(GroomedJetAK5.Pt()/GenWhad.Pt()-1.)>PtMatchingFraction) && 
	(deltaR(GroomedJetAK5_1.Phi(),GenWhad.Phi(),GroomedJetAK5_1.Eta(),GenWhad.Eta())<=DRMatchingCut && fabs(GroomedJetAK5_1.Pt()/GenWhad.Pt()-1.)<=PtMatchingFraction) ){

        DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching   ->Fill(GroomedJetAK5_1.Pt()/GenWhad.Pt());
        DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching  ->Fill(GroomedJetAK5_1.Eta()-GenWhad.Eta());
        DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching  ->Fill(deltaPhi(GroomedJetAK5_1.Phi(),GenWhad.Phi()));
        DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching    ->Fill(deltaR(GroomedJetAK5_1.Phi(),GenWhad.Phi(),GroomedJetAK5_1.Eta(),GenWhad.Eta()));
        DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching ->Fill(GroomedJetAK5_1.M()/GenWhad.M());   
    }



    if((deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta())>DRMatchingCut || fabs(GroomedJetAK5.Pt()/GenWhad.Pt()-1.)>PtMatchingFraction  ) &&
       (deltaR(GroomedJetAK5_1.Phi(),GenWhad.Phi(),GroomedJetAK5_1.Eta(),GenWhad.Eta())>DRMatchingCut || fabs(GroomedJetAK5_1.Pt()/GenWhad.Pt()-1.)>PtMatchingFraction )){
    
        DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching   ->Fill(GroomedJetAK5.Pt()/GenWhad.Pt());
        DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching  ->Fill(GroomedJetAK5.Eta()-GenWhad.Eta());
        DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching  ->Fill(deltaPhi(GroomedJetAK5.Phi(),GenWhad.Phi()));
        DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching    ->Fill(deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta()));
        DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching ->Fill(GroomedJetAK5.M()/GenWhad.M());   
     }

   }

   if(isExclusive_3Jet_CA8){

     nExclusive_3Jet_CA8++;

     TLorentzVector GroomedJetCA8_1 = TLorentzVector(0.,0.,0.,0.);
     TLorentzVector GroomedJetAK5_1 = TLorentzVector(0.,0.,0.,0.);
     GroomedJetCA8_1.SetPtEtaPhiE(GroomedJet_CA8_pt[1],GroomedJet_CA8_eta[1],GroomedJet_CA8_phi[1],GroomedJet_CA8_e[1]);
     GroomedJetAK5_1.SetPtEtaPhiE(GroomedJet_AK5_pt[1],GroomedJet_AK5_eta[1],GroomedJet_AK5_phi[1],GroomedJet_AK5_e[1]);

     TLorentzVector GroomedJetCA8_2 = TLorentzVector(0.,0.,0.,0.);
     TLorentzVector GroomedJetAK5_2 = TLorentzVector(0.,0.,0.,0.);
     GroomedJetCA8_2.SetPtEtaPhiE(GroomedJet_CA8_pt[2],GroomedJet_CA8_eta[2],GroomedJet_CA8_phi[2],GroomedJet_CA8_e[2]);
     GroomedJetAK5_2.SetPtEtaPhiE(GroomedJet_AK5_pt[2],GroomedJet_AK5_eta[2],GroomedJet_AK5_phi[2],GroomedJet_AK5_e[2]);

     // DR Matching CA8

     if(deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta())<=DRMatchingCut    && 
        deltaR(GroomedJetCA8_1.Phi(),GenWhad.Phi(),GroomedJetCA8_1.Eta(),GenWhad.Eta())>DRMatchingCut &&
        deltaR(GroomedJetCA8_2.Phi(),GenWhad.Phi(),GroomedJetCA8_2.Eta(),GenWhad.Eta())>DRMatchingCut    ){
     
       DPt_GroomedCA8_GenWhad_DRCut_3Jet  ->Fill(GroomedJetCA8.Pt()/GenWhad.Pt());
       DEta_GroomedCA8_GenWhad_DRCut_3Jet ->Fill(GroomedJetCA8.Eta()-GenWhad.Eta());
       DPhi_GroomedCA8_GenWhad_DRCut_3Jet ->Fill(deltaPhi(GroomedJetCA8.Phi(),GenWhad.Phi()));
       DR_GroomedCA8_GenWhad_DRCut_3Jet   ->Fill(deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta()));
       DMass_GroomedCA8_GenWhad_DRCut_3Jet->Fill(GroomedJetCA8.M()/GenWhad.M());

     }

     if( deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta())<=DRMatchingCut     &&
	(deltaR(GroomedJetCA8_1.Phi(),GenWhad.Phi(),GroomedJetCA8_1.Eta(),GenWhad.Eta())<=DRMatchingCut || 
         deltaR(GroomedJetCA8_2.Phi(),GenWhad.Phi(),GroomedJetCA8_2.Eta(),GenWhad.Eta())<=DRMatchingCut ) ){
     
       DPt_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching  ->Fill(GroomedJetCA8.Pt()/GenWhad.Pt());
       DEta_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching ->Fill(GroomedJetCA8.Eta()-GenWhad.Eta());
       DPhi_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching ->Fill(deltaPhi(GroomedJetCA8_1.Phi(),GenWhad.Phi()));
       DR_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching   ->Fill(deltaR(GroomedJetCA8_1.Phi(),GenWhad.Phi(),GroomedJetCA8_1.Eta(),GenWhad.Eta()));
       DMass_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->Fill(GroomedJetCA8.M()/GenWhad.M());
     }


     if(deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta())>DRMatchingCut){

       if( deltaR(GroomedJetCA8_1.Phi(),GenWhad.Phi(),GroomedJetCA8_1.Eta(),GenWhad.Eta())<=DRMatchingCut && 
           deltaR(GroomedJetCA8_2.Phi(),GenWhad.Phi(),GroomedJetCA8_2.Eta(),GenWhad.Eta())>DRMatchingCut ){
     
          DPt_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching  ->Fill(GroomedJetCA8_1.Pt()/GenWhad.Pt());
          DEta_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching ->Fill(GroomedJetCA8_1.Eta()-GenWhad.Eta());
          DPhi_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching ->Fill(deltaPhi(GroomedJetCA8_1.Phi(),GenWhad.Phi()));
          DR_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching   ->Fill(deltaR(GroomedJetCA8_1.Phi(),GenWhad.Phi(),GroomedJetCA8_1.Eta(),GenWhad.Eta()));
          DMass_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->Fill(GroomedJetCA8_1.M()/GenWhad.M());
      }

       if( deltaR(GroomedJetCA8_1.Phi(),GenWhad.Phi(),GroomedJetCA8_1.Eta(),GenWhad.Eta())>DRMatchingCut && 
           deltaR(GroomedJetCA8_2.Phi(),GenWhad.Phi(),GroomedJetCA8_2.Eta(),GenWhad.Eta())<=DRMatchingCut ){
     
          DPt_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching  ->Fill(GroomedJetCA8_2.Pt()/GenWhad.Pt());
          DEta_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching ->Fill(GroomedJetCA8_2.Eta()-GenWhad.Eta());
          DPhi_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching ->Fill(deltaPhi(GroomedJetCA8_2.Phi(),GenWhad.Phi()));
          DR_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching   ->Fill(deltaR(GroomedJetCA8_2.Phi(),GenWhad.Phi(),GroomedJetCA8_2.Eta(),GenWhad.Eta()));
          DMass_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->Fill(GroomedJetCA8_2.M()/GenWhad.M());
      }

     }

     if( deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta())>DRMatchingCut     &&
	 deltaR(GroomedJetCA8_1.Phi(),GenWhad.Phi(),GroomedJetCA8_1.Eta(),GenWhad.Eta())>DRMatchingCut &&
         deltaR(GroomedJetCA8_2.Phi(),GenWhad.Phi(),GroomedJetCA8_2.Eta(),GenWhad.Eta())>DRMatchingCut  ){
     
       DPt_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching  ->Fill(GroomedJetCA8.Pt()/GenWhad.Pt());
       DEta_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching ->Fill(GroomedJetCA8.Eta()-GenWhad.Eta());
       DPhi_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching ->Fill(deltaPhi(GroomedJetCA8_1.Phi(),GenWhad.Phi()));
       DR_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching   ->Fill(deltaR(GroomedJetCA8_1.Phi(),GenWhad.Phi(),GroomedJetCA8_1.Eta(),GenWhad.Eta()));
       DMass_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->Fill(GroomedJetCA8.M()/GenWhad.M());
     }

     // DR DPt Matching 

     if( (deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta())<=DRMatchingCut     && fabs(GroomedJetCA8.Pt()/GenWhad.Pt()-1.)<=PtMatchingFraction  ) && 
         ((deltaR(GroomedJetCA8_1.Phi(),GenWhad.Phi(),GroomedJetCA8_1.Eta(),GenWhad.Eta())>DRMatchingCut || fabs(GroomedJetCA8_1.Pt()/GenWhad.Pt()-1.)>PtMatchingFraction ) &&
	  (deltaR(GroomedJetCA8_2.Phi(),GenWhad.Phi(),GroomedJetCA8_2.Eta(),GenWhad.Eta())>DRMatchingCut || fabs(GroomedJetCA8_2.Pt()/GenWhad.Pt()-1.)>PtMatchingFraction ))  ){

        DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet ->Fill(GroomedJetCA8.Pt()/GenWhad.Pt());
        DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet ->Fill(GroomedJetCA8.Eta()-GenWhad.Eta());
        DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet ->Fill(deltaPhi(GroomedJetCA8.Phi(),GenWhad.Phi()));
        DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet   ->Fill(deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta()));
        DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet->Fill(GroomedJetCA8.M()/GenWhad.M());   
     }
     

     if( (deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta())<=DRMatchingCut && fabs(GroomedJetCA8.Pt()/GenWhad.Pt()-1.)<=PtMatchingFraction     ) && 
         ((deltaR(GroomedJetCA8_1.Phi(),GenWhad.Phi(),GroomedJetCA8_1.Eta(),GenWhad.Eta())<=DRMatchingCut && fabs(GroomedJetCA8_1.Pt()/GenWhad.Pt()-1.)<=PtMatchingFraction) ||
	  (deltaR(GroomedJetCA8_2.Phi(),GenWhad.Phi(),GroomedJetCA8_2.Eta(),GenWhad.Eta())<=DRMatchingCut && fabs(GroomedJetCA8_2.Pt()/GenWhad.Pt()-1.)<=PtMatchingFraction))  ){

        DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching  ->Fill(GroomedJetCA8.Pt()/GenWhad.Pt());
        DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching ->Fill(GroomedJetCA8.Eta()-GenWhad.Eta());
        DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching ->Fill(deltaPhi(GroomedJetCA8.Phi(),GenWhad.Phi()));
        DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching   ->Fill(deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta()));
        DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->Fill(GroomedJetCA8.M()/GenWhad.M());   
     }
     

     if(deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta())>DRMatchingCut && fabs(GroomedJetCA8.Pt()/GenWhad.Pt()-1.)>PtMatchingFraction ){

       if( ( deltaR(GroomedJetCA8_1.Phi(),GenWhad.Phi(),GroomedJetCA8_1.Eta(),GenWhad.Eta())<=DRMatchingCut && fabs(GroomedJetCA8_1.Pt()/GenWhad.Pt()-1.)<=PtMatchingFraction ) && 
           ( deltaR(GroomedJetCA8_2.Phi(),GenWhad.Phi(),GroomedJetCA8_2.Eta(),GenWhad.Eta())>DRMatchingCut || fabs(GroomedJetCA8_2.Pt()/GenWhad.Pt()-1.)>PtMatchingFraction ) ){
     
          DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching  ->Fill(GroomedJetCA8_1.Pt()/GenWhad.Pt());
          DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching ->Fill(GroomedJetCA8_1.Eta()-GenWhad.Eta());
          DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching ->Fill(deltaPhi(GroomedJetCA8_1.Phi(),GenWhad.Phi()));
          DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching   ->Fill(deltaR(GroomedJetCA8_1.Phi(),GenWhad.Phi(),GroomedJetCA8_1.Eta(),GenWhad.Eta()));
          DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->Fill(GroomedJetCA8_1.M()/GenWhad.M());
      }

       if( ( deltaR(GroomedJetCA8_1.Phi(),GenWhad.Phi(),GroomedJetCA8_1.Eta(),GenWhad.Eta())>DRMatchingCut || fabs(GroomedJetCA8_1.Pt()/GenWhad.Pt()-1.)>PtMatchingFraction ) && 
           ( deltaR(GroomedJetCA8_2.Phi(),GenWhad.Phi(),GroomedJetCA8_2.Eta(),GenWhad.Eta())<=DRMatchingCut && fabs(GroomedJetCA8_2.Pt()/GenWhad.Pt()-1.)<=PtMatchingFraction ) ){
     
          DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching  ->Fill(GroomedJetCA8_2.Pt()/GenWhad.Pt());
          DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching ->Fill(GroomedJetCA8_2.Eta()-GenWhad.Eta());
          DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching ->Fill(deltaPhi(GroomedJetCA8_2.Phi(),GenWhad.Phi()));
          DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching   ->Fill(deltaR(GroomedJetCA8_2.Phi(),GenWhad.Phi(),GroomedJetCA8_2.Eta(),GenWhad.Eta()));
          DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->Fill(GroomedJetCA8_2.M()/GenWhad.M());
      }

     }
 

     if( (deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta())>DRMatchingCut || fabs(GroomedJetCA8.Pt()/GenWhad.Pt()-1.)>PtMatchingFraction     ) && 
         ((deltaR(GroomedJetCA8_1.Phi(),GenWhad.Phi(),GroomedJetCA8_1.Eta(),GenWhad.Eta())>DRMatchingCut || fabs(GroomedJetCA8_1.Pt()/GenWhad.Pt()-1.)>PtMatchingFraction ) &&
	  (deltaR(GroomedJetCA8_2.Phi(),GenWhad.Phi(),GroomedJetCA8_2.Eta(),GenWhad.Eta())>DRMatchingCut || fabs(GroomedJetCA8_2.Pt()/GenWhad.Pt()-1.)>PtMatchingFraction ))  ){

        DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching ->Fill(GroomedJetCA8.Pt()/GenWhad.Pt());
        DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching ->Fill(GroomedJetCA8.Eta()-GenWhad.Eta());
        DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching ->Fill(deltaPhi(GroomedJetCA8.Phi(),GenWhad.Phi()));
        DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching   ->Fill(deltaR(GroomedJetCA8.Phi(),GenWhad.Phi(),GroomedJetCA8.Eta(),GenWhad.Eta()));
        DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->Fill(GroomedJetCA8.M()/GenWhad.M());   
     }
 
     // AK5 Jets 

     if(  deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta())<=DRMatchingCut && 
          deltaR(GroomedJetAK5_1.Phi(),GenWhad.Phi(),GroomedJetAK5_1.Eta(),GenWhad.Eta())>DRMatchingCut &&
          deltaR(GroomedJetAK5_2.Phi(),GenWhad.Phi(),GroomedJetAK5_2.Eta(),GenWhad.Eta())>DRMatchingCut    ){
     
       DPt_GroomedAK5_GenWhad_DRCut_3Jet  ->Fill(GroomedJetAK5.Pt()/GenWhad.Pt());
       DEta_GroomedAK5_GenWhad_DRCut_3Jet ->Fill(GroomedJetAK5.Eta()-GenWhad.Eta());
       DPhi_GroomedAK5_GenWhad_DRCut_3Jet ->Fill(deltaPhi(GroomedJetAK5.Phi(),GenWhad.Phi()));
       DR_GroomedAK5_GenWhad_DRCut_3Jet   ->Fill(deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta()));
       DMass_GroomedAK5_GenWhad_DRCut_3Jet->Fill(GroomedJetAK5.M()/GenWhad.M());
     }

     if(  deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta())<=DRMatchingCut && 
          ( deltaR(GroomedJetAK5_1.Phi(),GenWhad.Phi(),GroomedJetAK5_1.Eta(),GenWhad.Eta())<=DRMatchingCut ||
            deltaR(GroomedJetAK5_2.Phi(),GenWhad.Phi(),GroomedJetAK5_2.Eta(),GenWhad.Eta())<=DRMatchingCut )   ){
     
       DPt_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching  ->Fill(GroomedJetAK5.Pt()/GenWhad.Pt());
       DEta_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching ->Fill(GroomedJetAK5.Eta()-GenWhad.Eta());
       DPhi_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching ->Fill(deltaPhi(GroomedJetAK5.Phi(),GenWhad.Phi()));
       DR_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching   ->Fill(deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta()));
       DMass_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->Fill(GroomedJetAK5.M()/GenWhad.M());
     }


     if(deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta())>DRMatchingCut){

       if( deltaR(GroomedJetAK5_1.Phi(),GenWhad.Phi(),GroomedJetAK5_1.Eta(),GenWhad.Eta())<=DRMatchingCut && 
           deltaR(GroomedJetAK5_2.Phi(),GenWhad.Phi(),GroomedJetAK5_2.Eta(),GenWhad.Eta())>DRMatchingCut ){
     
          DPt_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching  ->Fill(GroomedJetAK5_1.Pt()/GenWhad.Pt());
          DEta_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching ->Fill(GroomedJetAK5_1.Eta()-GenWhad.Eta());
          DPhi_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching ->Fill(deltaPhi(GroomedJetAK5_1.Phi(),GenWhad.Phi()));
          DR_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching   ->Fill(deltaR(GroomedJetAK5_1.Phi(),GenWhad.Phi(),GroomedJetAK5_1.Eta(),GenWhad.Eta()));
          DMass_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->Fill(GroomedJetAK5_1.M()/GenWhad.M());
      }

       if( deltaR(GroomedJetAK5_1.Phi(),GenWhad.Phi(),GroomedJetAK5_1.Eta(),GenWhad.Eta())>DRMatchingCut && 
           deltaR(GroomedJetAK5_2.Phi(),GenWhad.Phi(),GroomedJetAK5_2.Eta(),GenWhad.Eta())<=DRMatchingCut ){
     
          DPt_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching  ->Fill(GroomedJetAK5_2.Pt()/GenWhad.Pt());
          DEta_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching ->Fill(GroomedJetAK5_2.Eta()-GenWhad.Eta());
          DPhi_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching ->Fill(deltaPhi(GroomedJetAK5_2.Phi(),GenWhad.Phi()));
          DR_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching   ->Fill(deltaR(GroomedJetAK5_2.Phi(),GenWhad.Phi(),GroomedJetAK5_2.Eta(),GenWhad.Eta()));
          DMass_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->Fill(GroomedJetAK5_2.M()/GenWhad.M());
      }

     }

     if( deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta())>DRMatchingCut &&
	 deltaR(GroomedJetAK5_1.Phi(),GenWhad.Phi(),GroomedJetAK5_1.Eta(),GenWhad.Eta())>DRMatchingCut &&
         deltaR(GroomedJetAK5_2.Phi(),GenWhad.Phi(),GroomedJetAK5_2.Eta(),GenWhad.Eta())>DRMatchingCut  ){
     
       DPt_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching  ->Fill(GroomedJetAK5.Pt()/GenWhad.Pt());
       DEta_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching ->Fill(GroomedJetAK5.Eta()-GenWhad.Eta());
       DPhi_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching ->Fill(deltaPhi(GroomedJetAK5_1.Phi(),GenWhad.Phi()));
       DR_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching   ->Fill(deltaR(GroomedJetAK5_1.Phi(),GenWhad.Phi(),GroomedJetAK5_1.Eta(),GenWhad.Eta()));
       DMass_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->Fill(GroomedJetAK5.M()/GenWhad.M());
     }


     if( (deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta())<=DRMatchingCut && fabs(GroomedJetAK5.Pt()/GenWhad.Pt()-1.)<=PtMatchingFraction     ) && 
         ((deltaR(GroomedJetAK5_1.Phi(),GenWhad.Phi(),GroomedJetAK5_1.Eta(),GenWhad.Eta())>DRMatchingCut || fabs(GroomedJetAK5_1.Pt()/GenWhad.Pt()-1.)>PtMatchingFraction ) &&
	  (deltaR(GroomedJetAK5_2.Phi(),GenWhad.Phi(),GroomedJetAK5_2.Eta(),GenWhad.Eta())>DRMatchingCut || fabs(GroomedJetAK5_2.Pt()/GenWhad.Pt()-1.)>PtMatchingFraction ) ) ){

        DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet ->Fill(GroomedJetAK5.Pt()/GenWhad.Pt());
        DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet ->Fill(GroomedJetAK5.Eta()-GenWhad.Eta());
        DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet ->Fill(deltaPhi(GroomedJetAK5.Phi(),GenWhad.Phi()));
        DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet   ->Fill(deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta()));
        DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet->Fill(GroomedJetAK5.M()/GenWhad.M());   
     }
     

    if( (deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta())<=DRMatchingCut && fabs(GroomedJetAK5.Pt()/GenWhad.Pt()-1.)<=PtMatchingFraction     ) && 
        ((deltaR(GroomedJetAK5_1.Phi(),GenWhad.Phi(),GroomedJetAK5_1.Eta(),GenWhad.Eta())<=DRMatchingCut && fabs(GroomedJetAK5_1.Pt()/GenWhad.Pt()-1.)<=PtMatchingFraction ) ||
         (deltaR(GroomedJetAK5_2.Phi(),GenWhad.Phi(),GroomedJetAK5_2.Eta(),GenWhad.Eta())<=DRMatchingCut && fabs(GroomedJetAK5_2.Pt()/GenWhad.Pt()-1.)<=PtMatchingFraction ) ) ){

        DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching ->Fill(GroomedJetAK5.Pt()/GenWhad.Pt());
        DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching ->Fill(GroomedJetAK5.Eta()-GenWhad.Eta());
        DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching ->Fill(deltaPhi(GroomedJetAK5.Phi(),GenWhad.Phi()));
        DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching   ->Fill(deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta()));
        DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->Fill(GroomedJetAK5.M()/GenWhad.M());   
     }

     if(deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta())>DRMatchingCut && fabs(GroomedJetAK5.Pt()/GenWhad.Pt()-1.)>PtMatchingFraction ){

       if( ( deltaR(GroomedJetAK5_1.Phi(),GenWhad.Phi(),GroomedJetAK5_1.Eta(),GenWhad.Eta())<=DRMatchingCut && fabs(GroomedJetAK5_1.Pt()/GenWhad.Pt()-1.)<=PtMatchingFraction ) && 
           ( deltaR(GroomedJetAK5_2.Phi(),GenWhad.Phi(),GroomedJetAK5_2.Eta(),GenWhad.Eta())>DRMatchingCut || fabs(GroomedJetAK5_2.Pt()/GenWhad.Pt()-1.)>PtMatchingFraction ) ){
     
          DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching  ->Fill(GroomedJetAK5_1.Pt()/GenWhad.Pt());
          DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching ->Fill(GroomedJetAK5_1.Eta()-GenWhad.Eta());
          DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching ->Fill(deltaPhi(GroomedJetAK5_1.Phi(),GenWhad.Phi()));
          DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching   ->Fill(deltaR(GroomedJetAK5_1.Phi(),GenWhad.Phi(),GroomedJetAK5_1.Eta(),GenWhad.Eta()));
          DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->Fill(GroomedJetAK5_1.M()/GenWhad.M());
      }

       if( ( deltaR(GroomedJetAK5_1.Phi(),GenWhad.Phi(),GroomedJetAK5_1.Eta(),GenWhad.Eta())>DRMatchingCut || fabs(GroomedJetAK5_1.Pt()/GenWhad.Pt()-1.)>PtMatchingFraction ) && 
           ( deltaR(GroomedJetAK5_2.Phi(),GenWhad.Phi(),GroomedJetAK5_2.Eta(),GenWhad.Eta())<=DRMatchingCut && fabs(GroomedJetAK5_2.Pt()/GenWhad.Pt()-1.)<=PtMatchingFraction ) ){
     
          DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching  ->Fill(GroomedJetAK5_2.Pt()/GenWhad.Pt());
          DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching ->Fill(GroomedJetAK5_2.Eta()-GenWhad.Eta());
          DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching ->Fill(deltaPhi(GroomedJetAK5_2.Phi(),GenWhad.Phi()));
          DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching   ->Fill(deltaR(GroomedJetAK5_2.Phi(),GenWhad.Phi(),GroomedJetAK5_2.Eta(),GenWhad.Eta()));
          DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->Fill(GroomedJetAK5_2.M()/GenWhad.M());
      }

     }

    if( (deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta())>DRMatchingCut || fabs(GroomedJetAK5.Pt()/GenWhad.Pt()-1.)>PtMatchingFraction     ) && 
         ((deltaR(GroomedJetAK5_1.Phi(),GenWhad.Phi(),GroomedJetAK5_1.Eta(),GenWhad.Eta())>DRMatchingCut || fabs(GroomedJetAK5_1.Pt()/GenWhad.Pt()-1.)>PtMatchingFraction ) &&
	  (deltaR(GroomedJetAK5_2.Phi(),GenWhad.Phi(),GroomedJetAK5_2.Eta(),GenWhad.Eta())>DRMatchingCut || fabs(GroomedJetAK5_2.Pt()/GenWhad.Pt()-1.)>PtMatchingFraction ))  ){

        DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching ->Fill(GroomedJetAK5.Pt()/GenWhad.Pt());
        DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching ->Fill(GroomedJetAK5.Eta()-GenWhad.Eta());
        DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching ->Fill(deltaPhi(GroomedJetAK5.Phi(),GenWhad.Phi()));
        DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching   ->Fill(deltaR(GroomedJetAK5.Phi(),GenWhad.Phi(),GroomedJetAK5.Eta(),GenWhad.Eta()));
        DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->Fill(GroomedJetAK5.M()/GenWhad.M());   
     }

     
   }

  }

  outputFile->cd();

  // Gen Jet Canvas for Final Plot

  TCanvas ** DPt_GenGroomedJet_AK5_GenWhad_c  = new TCanvas* [NJet];
  TCanvas ** DEta_GenGroomedJet_AK5_GenWhad_c = new TCanvas* [NJet];
  TCanvas ** DPhi_GenGroomedJet_AK5_GenWhad_c = new TCanvas* [NJet];
  TCanvas ** DR_GenGroomedJet_AK5_GenWhad_c   = new TCanvas* [NJet];
  TCanvas ** DE_GenGroomedJet_AK5_GenWhad_c = new TCanvas* [NJet];
  TCanvas ** DMass_GenGroomedJet_AK5_GenWhad_c   = new TCanvas* [NJet];

  TCanvas ** DPt_DR_GenGroomedJet_AK5_GenWhad_c = new TCanvas* [NJet];
  TCanvas ** DMass_DR_GenGroomedJet_AK5_GenWhad_c   = new TCanvas* [NJet];
 
  TCanvas ** DPt_GenGroomedJet_CA8_GenWhad_c  = new TCanvas* [NJet];
  TCanvas ** DEta_GenGroomedJet_CA8_GenWhad_c = new TCanvas* [NJet];
  TCanvas ** DPhi_GenGroomedJet_CA8_GenWhad_c = new TCanvas* [NJet];
  TCanvas ** DR_GenGroomedJet_CA8_GenWhad_c   = new TCanvas* [NJet];
  TCanvas ** DE_GenGroomedJet_CA8_GenWhad_c = new TCanvas* [NJet];
  TCanvas ** DMass_GenGroomedJet_CA8_GenWhad_c   = new TCanvas* [NJet];

  TCanvas ** DPt_DR_GenGroomedJet_CA8_GenWhad_c = new TCanvas* [NJet];
  TCanvas ** DMass_DR_GenGroomedJet_CA8_GenWhad_c   = new TCanvas* [NJet];

  TCanvas ** DPt_GenGroomedJet_CA8_subjet1_WParton_c  = new TCanvas* [NJet];
  TCanvas ** DEta_GenGroomedJet_CA8_subjet1_WParton_c = new TCanvas* [NJet];
  TCanvas ** DPhi_GenGroomedJet_CA8_subjet1_WParton_c = new TCanvas* [NJet];
  TCanvas ** DR_GenGroomedJet_CA8_subjet1_WParton_c   = new TCanvas* [NJet];
  TCanvas ** DE_GenGroomedJet_CA8_subjet1_WParton_c = new TCanvas* [NJet];
  TCanvas ** DMass_GenGroomedJet_CA8_subjet1_WParton_c   = new TCanvas* [NJet];

  TCanvas ** DPt_DR_GenGroomedJet_CA8_subjet1_WParton_c = new TCanvas* [NJet];
  TCanvas ** DMass_DR_GenGroomedJet_CA8_subjet1_WParton_c   = new TCanvas* [NJet];

  TCanvas ** DPt_GenGroomedJet_CA8_subjet2_WParton_c  = new TCanvas* [NJet];
  TCanvas ** DEta_GenGroomedJet_CA8_subjet2_WParton_c = new TCanvas* [NJet];
  TCanvas ** DPhi_GenGroomedJet_CA8_subjet2_WParton_c = new TCanvas* [NJet];
  TCanvas ** DR_GenGroomedJet_CA8_subjet2_WParton_c   = new TCanvas* [NJet];
  TCanvas ** DE_GenGroomedJet_CA8_subjet2_WParton_c = new TCanvas* [NJet];
  TCanvas ** DMass_GenGroomedJet_CA8_subjet2_WParton_c   = new TCanvas* [NJet];

  TCanvas ** DPt_DR_GenGroomedJet_CA8_subjet2_WParton_c = new TCanvas* [NJet];
  TCanvas ** DMass_DR_GenGroomedJet_CA8_subjet2_WParton_c   = new TCanvas* [NJet];

  TCanvas ** DPt_GenGroomedJet_CA8_GenWhad_pr_c  = new TCanvas* [NJet];
  TCanvas ** DEta_GenGroomedJet_CA8_GenWhad_pr_c = new TCanvas* [NJet];
  TCanvas ** DPhi_GenGroomedJet_CA8_GenWhad_pr_c = new TCanvas* [NJet];
  TCanvas ** DR_GenGroomedJet_CA8_GenWhad_pr_c   = new TCanvas* [NJet];
  TCanvas ** DE_GenGroomedJet_CA8_GenWhad_pr_c = new TCanvas* [NJet];
  TCanvas ** DMass_GenGroomedJet_CA8_GenWhad_pr_c   = new TCanvas* [NJet];

  TCanvas ** DPt_DR_GenGroomedJet_CA8_GenWhad_pr_c = new TCanvas* [NJet];
  TCanvas ** DMass_DR_GenGroomedJet_CA8_GenWhad_pr_c   = new TCanvas* [NJet];

  TCanvas ** DPt_GenGroomedJet_CA8_GenWhad_tr_c  = new TCanvas* [NJet];
  TCanvas ** DEta_GenGroomedJet_CA8_GenWhad_tr_c = new TCanvas* [NJet];
  TCanvas ** DPhi_GenGroomedJet_CA8_GenWhad_tr_c = new TCanvas* [NJet];
  TCanvas ** DR_GenGroomedJet_CA8_GenWhad_tr_c   = new TCanvas* [NJet];
  TCanvas ** DE_GenGroomedJet_CA8_GenWhad_tr_c = new TCanvas* [NJet];
  TCanvas ** DMass_GenGroomedJet_CA8_GenWhad_tr_c   = new TCanvas* [NJet];

  TCanvas ** DPt_DR_GenGroomedJet_CA8_GenWhad_tr_c = new TCanvas* [NJet];
  TCanvas ** DMass_DR_GenGroomedJet_CA8_GenWhad_tr_c   = new TCanvas* [NJet];

  TCanvas ** DPt_GenGroomedJet_CA8_GenWhad_ft_c  = new TCanvas* [NJet];
  TCanvas ** DEta_GenGroomedJet_CA8_GenWhad_ft_c = new TCanvas* [NJet];
  TCanvas ** DPhi_GenGroomedJet_CA8_GenWhad_ft_c = new TCanvas* [NJet];
  TCanvas ** DR_GenGroomedJet_CA8_GenWhad_ft_c   = new TCanvas* [NJet];
  TCanvas ** DE_GenGroomedJet_CA8_GenWhad_ft_c = new TCanvas* [NJet];
  TCanvas ** DMass_GenGroomedJet_CA8_GenWhad_ft_c   = new TCanvas* [NJet];

  TCanvas ** DPt_DR_GenGroomedJet_CA8_GenWhad_ft_c = new TCanvas* [NJet];
  TCanvas ** DMass_DR_GenGroomedJet_CA8_GenWhad_ft_c   = new TCanvas* [NJet];

  // Reco Jet Canvas for Final Plot

  TCanvas ** DPt_GroomedJet_AK5_GenWhad_c  = new TCanvas* [NJet];
  TCanvas ** DEta_GroomedJet_AK5_GenWhad_c = new TCanvas* [NJet];
  TCanvas ** DPhi_GroomedJet_AK5_GenWhad_c = new TCanvas* [NJet];
  TCanvas ** DR_GroomedJet_AK5_GenWhad_c   = new TCanvas* [NJet];
  TCanvas ** DE_GroomedJet_AK5_GenWhad_c = new TCanvas* [NJet];
  TCanvas ** DMass_GroomedJet_AK5_GenWhad_c   = new TCanvas* [NJet];

  TCanvas ** DPt_DR_GroomedJet_AK5_GenWhad_c = new TCanvas* [NJet];
  TCanvas ** DMass_DR_GroomedJet_AK5_GenWhad_c   = new TCanvas* [NJet];

  TCanvas ** DPt_GroomedJet_CA8_GenWhad_c  = new TCanvas* [NJet];
  TCanvas ** DEta_GroomedJet_CA8_GenWhad_c = new TCanvas* [NJet];
  TCanvas ** DPhi_GroomedJet_CA8_GenWhad_c = new TCanvas* [NJet];
  TCanvas ** DR_GroomedJet_CA8_GenWhad_c   = new TCanvas* [NJet];
  TCanvas ** DE_GroomedJet_CA8_GenWhad_c = new TCanvas* [NJet];
  TCanvas ** DMass_GroomedJet_CA8_GenWhad_c   = new TCanvas* [NJet];

  TCanvas ** DPt_DR_GroomedJet_CA8_GenWhad_c = new TCanvas* [NJet];
  TCanvas ** DMass_DR_GroomedJet_CA8_GenWhad_c   = new TCanvas* [NJet];

  TCanvas ** DPt_GroomedJet_CA8_subjet1_WParton_c  = new TCanvas* [NJet];
  TCanvas ** DEta_GroomedJet_CA8_subjet1_WParton_c = new TCanvas* [NJet];
  TCanvas ** DPhi_GroomedJet_CA8_subjet1_WParton_c = new TCanvas* [NJet];
  TCanvas ** DR_GroomedJet_CA8_subjet1_WParton_c   = new TCanvas* [NJet];
  TCanvas ** DE_GroomedJet_CA8_subjet1_WParton_c = new TCanvas* [NJet];
  TCanvas ** DMass_GroomedJet_CA8_subjet1_WParton_c   = new TCanvas* [NJet];

  TCanvas ** DPt_DR_GroomedJet_CA8_subjet1_WParton_c = new TCanvas* [NJet];
  TCanvas ** DMass_DR_GroomedJet_CA8_subjet1_WParton_c   = new TCanvas* [NJet];

  TCanvas ** DPt_GroomedJet_CA8_subjet2_WParton_c  = new TCanvas* [NJet];
  TCanvas ** DEta_GroomedJet_CA8_subjet2_WParton_c = new TCanvas* [NJet];
  TCanvas ** DPhi_GroomedJet_CA8_subjet2_WParton_c = new TCanvas* [NJet];
  TCanvas ** DR_GroomedJet_CA8_subjet2_WParton_c   = new TCanvas* [NJet];
  TCanvas ** DE_GroomedJet_CA8_subjet2_WParton_c = new TCanvas* [NJet];
  TCanvas ** DMass_GroomedJet_CA8_subjet2_WParton_c   = new TCanvas* [NJet];

  TCanvas ** DPt_DR_GroomedJet_CA8_subjet2_WParton_c = new TCanvas* [NJet];
  TCanvas ** DMass_DR_GroomedJet_CA8_subjet2_WParton_c   = new TCanvas* [NJet];

  TCanvas ** DPt_GroomedJet_CA8_GenWhad_pr_c  = new TCanvas* [NJet];
  TCanvas ** DEta_GroomedJet_CA8_GenWhad_pr_c = new TCanvas* [NJet];
  TCanvas ** DPhi_GroomedJet_CA8_GenWhad_pr_c = new TCanvas* [NJet];
  TCanvas ** DR_GroomedJet_CA8_GenWhad_pr_c   = new TCanvas* [NJet];
  TCanvas ** DE_GroomedJet_CA8_GenWhad_pr_c = new TCanvas* [NJet];
  TCanvas ** DMass_GroomedJet_CA8_GenWhad_pr_c   = new TCanvas* [NJet];

  TCanvas ** DPt_DR_GroomedJet_CA8_GenWhad_pr_c = new TCanvas* [NJet];
  TCanvas ** DMass_DR_GroomedJet_CA8_GenWhad_pr_c   = new TCanvas* [NJet];

  TCanvas ** DPt_GroomedJet_CA8_GenWhad_tr_c  = new TCanvas* [NJet];
  TCanvas ** DEta_GroomedJet_CA8_GenWhad_tr_c = new TCanvas* [NJet];
  TCanvas ** DPhi_GroomedJet_CA8_GenWhad_tr_c = new TCanvas* [NJet];
  TCanvas ** DR_GroomedJet_CA8_GenWhad_tr_c   = new TCanvas* [NJet];
  TCanvas ** DE_GroomedJet_CA8_GenWhad_tr_c = new TCanvas* [NJet];
  TCanvas ** DMass_GroomedJet_CA8_GenWhad_tr_c   = new TCanvas* [NJet];

  TCanvas ** DPt_DR_GroomedJet_CA8_GenWhad_tr_c = new TCanvas* [NJet];
  TCanvas ** DMass_DR_GroomedJet_CA8_GenWhad_tr_c   = new TCanvas* [NJet];

  TCanvas ** DPt_GroomedJet_CA8_GenWhad_ft_c  = new TCanvas* [NJet];
  TCanvas ** DEta_GroomedJet_CA8_GenWhad_ft_c = new TCanvas* [NJet];
  TCanvas ** DPhi_GroomedJet_CA8_GenWhad_ft_c = new TCanvas* [NJet];
  TCanvas ** DR_GroomedJet_CA8_GenWhad_ft_c   = new TCanvas* [NJet];
  TCanvas ** DE_GroomedJet_CA8_GenWhad_ft_c = new TCanvas* [NJet];
  TCanvas ** DMass_GroomedJet_CA8_GenWhad_ft_c   = new TCanvas* [NJet];

  TCanvas ** DPt_DR_GroomedJet_CA8_GenWhad_ft_c = new TCanvas* [NJet];
  TCanvas ** DMass_DR_GroomedJet_CA8_GenWhad_ft_c   = new TCanvas* [NJet];

  TCanvas ** DPt_GroomedJet_CA8_WLep_c  = new TCanvas* [NJet];
  TCanvas ** DEta_GroomedJet_CA8_WLep_c = new TCanvas* [NJet];
  TCanvas ** DPhi_GroomedJet_CA8_WLep_c = new TCanvas* [NJet];
  TCanvas ** DR_GroomedJet_CA8_WLep_c   = new TCanvas* [NJet];
  TCanvas ** DE_GroomedJet_CA8_WLep_c = new TCanvas* [NJet];
  TCanvas ** DMass_GroomedJet_CA8_WLep_c   = new TCanvas* [NJet];

  TCanvas ** DPt_DR_GroomedJet_CA8_WLep_c     = new TCanvas* [NJet];
  TCanvas ** DMass_DR_GroomedJet_CA8_WLep_c   = new TCanvas* [NJet];
  
 
  for(int i =0; i<NJet ;i++){

   // Gen AK5

   DPt_GenGroomedJet_AK5_GenWhad_c[i]  = new TCanvas(DPt_GenGroomedJet_AK5_GenWhad[i]->GetName());
   DPt_GenGroomedJet_AK5_GenWhad_c[i]->cd();
   DPt_GenGroomedJet_AK5_GenWhad_c[i]->SetLogy(); 
   DPt_GenGroomedJet_AK5_GenWhad[i]->Draw();
   DPt_GenGroomedJet_AK5_GenWhad[i]->Write();
   DPt_GenGroomedJet_AK5_GenWhad_c[i]->Print((OutputPlotDir+"/"+DPt_GenGroomedJet_AK5_GenWhad[i]->GetName()+".png").c_str(),"png");

   DEta_GenGroomedJet_AK5_GenWhad_c[i] = new TCanvas(DEta_GenGroomedJet_AK5_GenWhad[i]->GetName());
   DEta_GenGroomedJet_AK5_GenWhad_c[i]->cd();
   DEta_GenGroomedJet_AK5_GenWhad[i]->Draw();
   DEta_GenGroomedJet_AK5_GenWhad[i]->Write();
   DEta_GenGroomedJet_AK5_GenWhad_c[i]->Print((OutputPlotDir+"/"+DEta_GenGroomedJet_AK5_GenWhad[i]->GetName()+".png").c_str(),"png");

   DPhi_GenGroomedJet_AK5_GenWhad_c[i] = new TCanvas(DPhi_GenGroomedJet_AK5_GenWhad[i]->GetName());
   DPhi_GenGroomedJet_AK5_GenWhad_c[i]->cd();
   DPhi_GenGroomedJet_AK5_GenWhad[i]->Draw();
   DPhi_GenGroomedJet_AK5_GenWhad[i]->Write();
   DPhi_GenGroomedJet_AK5_GenWhad_c[i]->Print((OutputPlotDir+"/"+DPhi_GenGroomedJet_AK5_GenWhad[i]->GetName()+".png").c_str(),"png");

   DR_GenGroomedJet_AK5_GenWhad_c[i]   = new TCanvas(DR_GenGroomedJet_AK5_GenWhad[i]->GetName());
   DR_GenGroomedJet_AK5_GenWhad_c[i]->cd();
   DR_GenGroomedJet_AK5_GenWhad_c[i]->SetLogy(); 
   DR_GenGroomedJet_AK5_GenWhad[i]->Draw();
   DR_GenGroomedJet_AK5_GenWhad[i]->Write();
   DR_GenGroomedJet_AK5_GenWhad_c[i]->Print((OutputPlotDir+"/"+DR_GenGroomedJet_AK5_GenWhad[i]->GetName()+".png").c_str(),"png");

   DE_GenGroomedJet_AK5_GenWhad_c[i]   = new TCanvas(DE_GenGroomedJet_AK5_GenWhad[i]->GetName());
   DE_GenGroomedJet_AK5_GenWhad_c[i]->cd();
   DE_GenGroomedJet_AK5_GenWhad_c[i]->SetLogy();
   DE_GenGroomedJet_AK5_GenWhad[i]->Draw();
   DE_GenGroomedJet_AK5_GenWhad[i]->Write();
   DE_GenGroomedJet_AK5_GenWhad_c[i]->Print((OutputPlotDir+"/"+DE_GenGroomedJet_AK5_GenWhad[i]->GetName()+".png").c_str(),"png");

   DMass_GenGroomedJet_AK5_GenWhad_c[i]   = new TCanvas(DMass_GenGroomedJet_AK5_GenWhad[i]->GetName());
   DMass_GenGroomedJet_AK5_GenWhad_c[i]->cd();
   DMass_GenGroomedJet_AK5_GenWhad_c[i]->SetLogy();
   DMass_GenGroomedJet_AK5_GenWhad[i]->Draw();
   DMass_GenGroomedJet_AK5_GenWhad[i]->Write();
   DMass_GenGroomedJet_AK5_GenWhad_c[i]->Print((OutputPlotDir+"/"+DMass_GenGroomedJet_AK5_GenWhad[i]->GetName()+".png").c_str(),"png");

   DPt_DR_GenGroomedJet_AK5_GenWhad_c[i]   = new TCanvas(DPt_DR_GenGroomedJet_AK5_GenWhad[i]->GetName());
   DPt_DR_GenGroomedJet_AK5_GenWhad_c[i]->cd();
   DPt_DR_GenGroomedJet_AK5_GenWhad[i]->Draw("colz");
   DPt_DR_GenGroomedJet_AK5_GenWhad[i]->Write();
   DPt_DR_GenGroomedJet_AK5_GenWhad_c[i]->Print((OutputPlotDir+"/"+DPt_DR_GenGroomedJet_AK5_GenWhad[i]->GetName()+".png").c_str(),"png");

   DMass_DR_GenGroomedJet_AK5_GenWhad_c[i]   = new TCanvas(DMass_DR_GenGroomedJet_AK5_GenWhad[i]->GetName());
   DMass_DR_GenGroomedJet_AK5_GenWhad_c[i]->cd();
   DMass_DR_GenGroomedJet_AK5_GenWhad[i]->Draw("colz");
   DMass_DR_GenGroomedJet_AK5_GenWhad[i]->Write();
   DMass_DR_GenGroomedJet_AK5_GenWhad_c[i]->Print((OutputPlotDir+"/"+DMass_DR_GenGroomedJet_AK5_GenWhad[i]->GetName()+".png").c_str(),"png");


   // Gen Groomed CA8

   DPt_GenGroomedJet_CA8_GenWhad_c[i]  = new TCanvas(DPt_GenGroomedJet_CA8_GenWhad[i]->GetName());
   DPt_GenGroomedJet_CA8_GenWhad_c[i]->cd();
   DPt_GenGroomedJet_CA8_GenWhad_c[i]->SetLogy();
   DPt_GenGroomedJet_CA8_GenWhad[i]->Draw();
   DPt_GenGroomedJet_CA8_GenWhad[i]->Write();
   DPt_GenGroomedJet_CA8_GenWhad_c[i]->Print((OutputPlotDir+"/"+DPt_GenGroomedJet_CA8_GenWhad[i]->GetName()+".png").c_str(),"png");

   DEta_GenGroomedJet_CA8_GenWhad_c[i] = new TCanvas(DEta_GenGroomedJet_CA8_GenWhad[i]->GetName());
   DEta_GenGroomedJet_CA8_GenWhad_c[i]->cd();
   DEta_GenGroomedJet_CA8_GenWhad[i]->Draw();
   DEta_GenGroomedJet_CA8_GenWhad[i]->Write();
   DEta_GenGroomedJet_CA8_GenWhad_c[i]->Print((OutputPlotDir+"/"+DEta_GenGroomedJet_CA8_GenWhad[i]->GetName()+".png").c_str(),"png");

   DPhi_GenGroomedJet_CA8_GenWhad_c[i] = new TCanvas(DPhi_GenGroomedJet_CA8_GenWhad[i]->GetName());
   DPhi_GenGroomedJet_CA8_GenWhad_c[i]->cd();
   DPhi_GenGroomedJet_CA8_GenWhad[i]->Draw();
   DPhi_GenGroomedJet_CA8_GenWhad[i]->Write();
   DPhi_GenGroomedJet_CA8_GenWhad_c[i]->Print((OutputPlotDir+"/"+DPhi_GenGroomedJet_CA8_GenWhad[i]->GetName()+".png").c_str(),"png");

   DR_GenGroomedJet_CA8_GenWhad_c[i]   = new TCanvas(DR_GenGroomedJet_CA8_GenWhad[i]->GetName());
   DR_GenGroomedJet_CA8_GenWhad_c[i]->cd();
   DR_GenGroomedJet_CA8_GenWhad_c[i]->SetLogy();
   DR_GenGroomedJet_CA8_GenWhad[i]->Draw();
   DR_GenGroomedJet_CA8_GenWhad[i]->Write();
   DR_GenGroomedJet_CA8_GenWhad_c[i]->Print((OutputPlotDir+"/"+DR_GenGroomedJet_CA8_GenWhad[i]->GetName()+".png").c_str(),"png");

   DE_GenGroomedJet_CA8_GenWhad_c[i]   = new TCanvas(DE_GenGroomedJet_CA8_GenWhad[i]->GetName());
   DE_GenGroomedJet_CA8_GenWhad_c[i]->cd();
   DE_GenGroomedJet_CA8_GenWhad_c[i]->SetLogy();
   DE_GenGroomedJet_CA8_GenWhad[i]->Draw();
   DE_GenGroomedJet_CA8_GenWhad[i]->Write();
   DE_GenGroomedJet_CA8_GenWhad_c[i]->Print((OutputPlotDir+"/"+DE_GenGroomedJet_CA8_GenWhad[i]->GetName()+".png").c_str(),"png");

   DMass_GenGroomedJet_CA8_GenWhad_c[i]   = new TCanvas(DMass_GenGroomedJet_CA8_GenWhad[i]->GetName());
   DMass_GenGroomedJet_CA8_GenWhad_c[i]->cd();
   DMass_GenGroomedJet_CA8_GenWhad_c[i]->SetLogy();
   DMass_GenGroomedJet_CA8_GenWhad[i]->Draw();
   DMass_GenGroomedJet_CA8_GenWhad[i]->Write();
   DMass_GenGroomedJet_CA8_GenWhad_c[i]->Print((OutputPlotDir+"/"+DMass_GenGroomedJet_CA8_GenWhad[i]->GetName()+".png").c_str(),"png");

   DPt_DR_GenGroomedJet_CA8_GenWhad_c[i]   = new TCanvas(DPt_DR_GenGroomedJet_CA8_GenWhad[i]->GetName());
   DPt_DR_GenGroomedJet_CA8_GenWhad_c[i]->cd();
   DPt_DR_GenGroomedJet_CA8_GenWhad[i]->Draw("colz");
   DPt_DR_GenGroomedJet_CA8_GenWhad[i]->Write();
   DPt_DR_GenGroomedJet_CA8_GenWhad_c[i]->Print((OutputPlotDir+"/"+DPt_DR_GenGroomedJet_CA8_GenWhad[i]->GetName()+".png").c_str(),"png");

   DMass_DR_GenGroomedJet_CA8_GenWhad_c[i]   = new TCanvas(DMass_DR_GenGroomedJet_CA8_GenWhad[i]->GetName());
   DMass_DR_GenGroomedJet_CA8_GenWhad_c[i]->cd();
   DMass_DR_GenGroomedJet_CA8_GenWhad[i]->Draw("colz");
   DMass_DR_GenGroomedJet_CA8_GenWhad[i]->Write();
   DMass_DR_GenGroomedJet_CA8_GenWhad_c[i]->Print((OutputPlotDir+"/"+DMass_DR_GenGroomedJet_CA8_GenWhad[i]->GetName()+".png").c_str(),"png");

   DPt_GenGroomedJet_CA8_subjet1_WParton_c[i]  = new TCanvas(DPt_GenGroomedJet_CA8_subjet1_WParton[i]->GetName());
   DPt_GenGroomedJet_CA8_subjet1_WParton_c[i]->cd();
   DPt_GenGroomedJet_CA8_subjet1_WParton_c[i]->SetLogy();
   DPt_GenGroomedJet_CA8_subjet1_WParton[i]->Draw();
   DPt_GenGroomedJet_CA8_subjet1_WParton[i]->Write();
   DPt_GenGroomedJet_CA8_subjet1_WParton_c[i]->Print((OutputPlotDir+"/"+DPt_GenGroomedJet_CA8_subjet1_WParton[i]->GetName()+".png").c_str(),"png");

   DEta_GenGroomedJet_CA8_subjet1_WParton_c[i] = new TCanvas(DEta_GenGroomedJet_CA8_subjet1_WParton[i]->GetName());
   DEta_GenGroomedJet_CA8_subjet1_WParton_c[i]->cd();
   DEta_GenGroomedJet_CA8_subjet1_WParton[i]->Draw();
   DEta_GenGroomedJet_CA8_subjet1_WParton[i]->Write();
   DEta_GenGroomedJet_CA8_subjet1_WParton_c[i]->Print((OutputPlotDir+"/"+DEta_GenGroomedJet_CA8_subjet1_WParton[i]->GetName()+".png").c_str(),"png");

   DPhi_GenGroomedJet_CA8_subjet1_WParton_c[i] = new TCanvas(DPhi_GenGroomedJet_CA8_subjet1_WParton[i]->GetName());
   DPhi_GenGroomedJet_CA8_subjet1_WParton[i]->Draw();
   DPhi_GenGroomedJet_CA8_subjet1_WParton[i]->Write();
   DPhi_GenGroomedJet_CA8_subjet1_WParton_c[i]->Print((OutputPlotDir+"/"+DPhi_GenGroomedJet_CA8_subjet1_WParton[i]->GetName()+".png").c_str(),"png");

   DR_GenGroomedJet_CA8_subjet1_WParton_c[i]   = new TCanvas(DR_GenGroomedJet_CA8_subjet1_WParton[i]->GetName());
   DR_GenGroomedJet_CA8_subjet1_WParton_c[i]->cd();
   DR_GenGroomedJet_CA8_subjet1_WParton_c[i]->SetLogy();
   DR_GenGroomedJet_CA8_subjet1_WParton[i]->Draw();
   DR_GenGroomedJet_CA8_subjet1_WParton[i]->Write();
   DR_GenGroomedJet_CA8_subjet1_WParton_c[i]->Print((OutputPlotDir+"/"+DR_GenGroomedJet_CA8_subjet1_WParton[i]->GetName()+".png").c_str(),"png");

   DE_GenGroomedJet_CA8_subjet1_WParton_c[i]   = new TCanvas(DE_GenGroomedJet_CA8_subjet1_WParton[i]->GetName());
   DE_GenGroomedJet_CA8_subjet1_WParton_c[i]->cd();
   DE_GenGroomedJet_CA8_subjet1_WParton_c[i]->SetLogy();
   DE_GenGroomedJet_CA8_subjet1_WParton[i]->Draw();
   DE_GenGroomedJet_CA8_subjet1_WParton[i]->Write();
   DE_GenGroomedJet_CA8_subjet1_WParton_c[i]->Print((OutputPlotDir+"/"+DE_GenGroomedJet_CA8_subjet1_WParton[i]->GetName()+".png").c_str(),"png");

   DMass_GenGroomedJet_CA8_subjet1_WParton_c[i]   = new TCanvas(DMass_GenGroomedJet_CA8_subjet1_WParton[i]->GetName());
   DMass_GenGroomedJet_CA8_subjet1_WParton_c[i]->cd();
   DMass_GenGroomedJet_CA8_subjet1_WParton_c[i]->SetLogy();
   DMass_GenGroomedJet_CA8_subjet1_WParton[i]->Draw();
   DMass_GenGroomedJet_CA8_subjet1_WParton[i]->Write();
   DMass_GenGroomedJet_CA8_subjet1_WParton_c[i]->Print((OutputPlotDir+"/"+DMass_GenGroomedJet_CA8_subjet1_WParton[i]->GetName()+".png").c_str(),"png");

   DPt_DR_GenGroomedJet_CA8_subjet1_WParton_c[i]   = new TCanvas(DPt_DR_GenGroomedJet_CA8_subjet1_WParton[i]->GetName());
   DPt_DR_GenGroomedJet_CA8_subjet1_WParton_c[i]->cd();
   DPt_DR_GenGroomedJet_CA8_subjet1_WParton[i]->Draw("colz");
   DPt_DR_GenGroomedJet_CA8_subjet1_WParton[i]->Write();
   DPt_DR_GenGroomedJet_CA8_subjet1_WParton_c[i]->Print((OutputPlotDir+"/"+DPt_DR_GenGroomedJet_CA8_subjet1_WParton[i]->GetName()+".png").c_str(),"png");

   DMass_DR_GenGroomedJet_CA8_subjet1_WParton_c[i]   = new TCanvas(DMass_DR_GenGroomedJet_CA8_subjet1_WParton[i]->GetName());
   DMass_DR_GenGroomedJet_CA8_subjet1_WParton_c[i]->cd();
   DMass_DR_GenGroomedJet_CA8_subjet1_WParton[i]->Draw("colz");
   DMass_DR_GenGroomedJet_CA8_subjet1_WParton[i]->Write();
   DMass_DR_GenGroomedJet_CA8_subjet1_WParton_c[i]->Print((OutputPlotDir+"/"+DMass_DR_GenGroomedJet_CA8_subjet1_WParton[i]->GetName()+".png").c_str(),"png");

   DPt_GenGroomedJet_CA8_subjet2_WParton_c[i]  = new TCanvas(DPt_GenGroomedJet_CA8_subjet2_WParton[i]->GetName());
   DPt_GenGroomedJet_CA8_subjet2_WParton_c[i]->cd();
   DPt_GenGroomedJet_CA8_subjet2_WParton_c[i]->SetLogy();
   DPt_GenGroomedJet_CA8_subjet2_WParton[i]->Draw();
   DPt_GenGroomedJet_CA8_subjet2_WParton[i]->Write();
   DPt_GenGroomedJet_CA8_subjet2_WParton_c[i]->Print((OutputPlotDir+"/"+DPt_GenGroomedJet_CA8_subjet2_WParton[i]->GetName()+".png").c_str(),"png");

   DEta_GenGroomedJet_CA8_subjet2_WParton_c[i] = new TCanvas(DEta_GenGroomedJet_CA8_subjet2_WParton[i]->GetName());
   DEta_GenGroomedJet_CA8_subjet2_WParton_c[i]->cd();
   DEta_GenGroomedJet_CA8_subjet2_WParton[i]->Draw();
   DEta_GenGroomedJet_CA8_subjet2_WParton[i]->Write();
   DEta_GenGroomedJet_CA8_subjet2_WParton_c[i]->Print((OutputPlotDir+"/"+DEta_GenGroomedJet_CA8_subjet2_WParton[i]->GetName()+".png").c_str(),"png");

   DPhi_GenGroomedJet_CA8_subjet2_WParton_c[i] = new TCanvas(DPhi_GenGroomedJet_CA8_subjet2_WParton[i]->GetName());
   DPhi_GenGroomedJet_CA8_subjet2_WParton_c[i]->cd();
   DPhi_GenGroomedJet_CA8_subjet2_WParton[i]->Draw();
   DPhi_GenGroomedJet_CA8_subjet2_WParton[i]->Write();
   DPhi_GenGroomedJet_CA8_subjet2_WParton_c[i]->Print((OutputPlotDir+"/"+DPhi_GenGroomedJet_CA8_subjet2_WParton[i]->GetName()+".png").c_str(),"png");

   DR_GenGroomedJet_CA8_subjet2_WParton_c[i]   = new TCanvas(DR_GenGroomedJet_CA8_subjet2_WParton[i]->GetName());
   DR_GenGroomedJet_CA8_subjet2_WParton_c[i]->cd();
   DR_GenGroomedJet_CA8_subjet2_WParton_c[i]->SetLogy();
   DR_GenGroomedJet_CA8_subjet2_WParton[i]->Draw();
   DR_GenGroomedJet_CA8_subjet2_WParton[i]->Write();
   DR_GenGroomedJet_CA8_subjet2_WParton_c[i]->Print((OutputPlotDir+"/"+DR_GenGroomedJet_CA8_subjet2_WParton[i]->GetName()+".png").c_str(),"png");

   DE_GenGroomedJet_CA8_subjet2_WParton_c[i]   = new TCanvas(DE_GenGroomedJet_CA8_subjet2_WParton[i]->GetName());
   DE_GenGroomedJet_CA8_subjet2_WParton_c[i]->cd();
   DE_GenGroomedJet_CA8_subjet2_WParton_c[i]->SetLogy();
   DE_GenGroomedJet_CA8_subjet2_WParton[i]->Draw();
   DE_GenGroomedJet_CA8_subjet2_WParton[i]->Write();
   DE_GenGroomedJet_CA8_subjet2_WParton_c[i]->Print((OutputPlotDir+"/"+DE_GenGroomedJet_CA8_subjet2_WParton[i]->GetName()+".png").c_str(),"png");

   DMass_GenGroomedJet_CA8_subjet2_WParton_c[i]   = new TCanvas(DMass_GenGroomedJet_CA8_subjet2_WParton[i]->GetName());
   DMass_GenGroomedJet_CA8_subjet2_WParton_c[i]->cd();
   DMass_GenGroomedJet_CA8_subjet2_WParton_c[i]->SetLogy();
   DMass_GenGroomedJet_CA8_subjet2_WParton[i]->Draw();
   DMass_GenGroomedJet_CA8_subjet2_WParton[i]->Write();
   DMass_GenGroomedJet_CA8_subjet2_WParton_c[i]->Print((OutputPlotDir+"/"+DMass_GenGroomedJet_CA8_subjet2_WParton[i]->GetName()+".png").c_str(),"png");

   DPt_DR_GenGroomedJet_CA8_subjet2_WParton_c[i]   = new TCanvas(DPt_DR_GenGroomedJet_CA8_subjet2_WParton[i]->GetName());
   DPt_DR_GenGroomedJet_CA8_subjet2_WParton_c[i]->cd();
   DPt_DR_GenGroomedJet_CA8_subjet2_WParton[i]->Draw("colz");
   DPt_DR_GenGroomedJet_CA8_subjet2_WParton[i]->Write();
   DPt_DR_GenGroomedJet_CA8_subjet2_WParton_c[i]->Print((OutputPlotDir+"/"+DPt_DR_GenGroomedJet_CA8_subjet2_WParton[i]->GetName()+".png").c_str(),"png");

   DMass_DR_GenGroomedJet_CA8_subjet2_WParton_c[i]   = new TCanvas(DMass_DR_GenGroomedJet_CA8_subjet2_WParton[i]->GetName());
   DMass_DR_GenGroomedJet_CA8_subjet2_WParton_c[i]->cd();
   DMass_DR_GenGroomedJet_CA8_subjet2_WParton[i]->Draw("colz");
   DMass_DR_GenGroomedJet_CA8_subjet2_WParton[i]->Write();
   DMass_DR_GenGroomedJet_CA8_subjet2_WParton_c[i]->Print((OutputPlotDir+"/"+DMass_DR_GenGroomedJet_CA8_subjet2_WParton[i]->GetName()+".png").c_str(),"png");

   DPt_GenGroomedJet_CA8_GenWhad_pr_c[i]  = new TCanvas(DPt_GenGroomedJet_CA8_GenWhad_pr[i]->GetName());
   DPt_GenGroomedJet_CA8_GenWhad_pr_c[i]->cd();
   DPt_GenGroomedJet_CA8_GenWhad_pr_c[i]->SetLogy();
   DPt_GenGroomedJet_CA8_GenWhad_pr[i]->Draw();
   DPt_GenGroomedJet_CA8_GenWhad_pr[i]->Write();
   DPt_GenGroomedJet_CA8_GenWhad_pr_c[i]->Print((OutputPlotDir+"/"+DPt_GenGroomedJet_CA8_GenWhad_pr[i]->GetName()+".png").c_str(),"png");

   DEta_GenGroomedJet_CA8_GenWhad_pr_c[i] = new TCanvas(DEta_GenGroomedJet_CA8_GenWhad_pr[i]->GetName());
   DEta_GenGroomedJet_CA8_GenWhad_pr_c[i]->cd();
   DEta_GenGroomedJet_CA8_GenWhad_pr[i]->Draw();
   DEta_GenGroomedJet_CA8_GenWhad_pr[i]->Write();
   DEta_GenGroomedJet_CA8_GenWhad_pr_c[i]->Print((OutputPlotDir+"/"+DEta_GenGroomedJet_CA8_GenWhad_pr[i]->GetName()+".png").c_str(),"png");

   DPhi_GenGroomedJet_CA8_GenWhad_pr_c[i] = new TCanvas(DPhi_GenGroomedJet_CA8_GenWhad_pr[i]->GetName());
   DPhi_GenGroomedJet_CA8_GenWhad_pr_c[i]->cd();
   DPhi_GenGroomedJet_CA8_GenWhad_pr[i]->Draw();
   DPhi_GenGroomedJet_CA8_GenWhad_pr[i]->Write();
   DPhi_GenGroomedJet_CA8_GenWhad_pr_c[i]->Print((OutputPlotDir+"/"+DPhi_GenGroomedJet_CA8_GenWhad_pr[i]->GetName()+".png").c_str(),"png");

   DR_GenGroomedJet_CA8_GenWhad_pr_c[i]   = new TCanvas(DR_GenGroomedJet_CA8_GenWhad_pr[i]->GetName());
   DR_GenGroomedJet_CA8_GenWhad_pr_c[i]->cd();
   DR_GenGroomedJet_CA8_GenWhad_pr_c[i]->SetLogy();
   DR_GenGroomedJet_CA8_GenWhad_pr[i]->Draw();
   DR_GenGroomedJet_CA8_GenWhad_pr[i]->Write();
   DR_GenGroomedJet_CA8_GenWhad_pr_c[i]->Print((OutputPlotDir+"/"+DR_GenGroomedJet_CA8_GenWhad_pr[i]->GetName()+".png").c_str(),"png");

   DE_GenGroomedJet_CA8_GenWhad_pr_c[i]   = new TCanvas(DE_GenGroomedJet_CA8_GenWhad_pr[i]->GetName());
   DE_GenGroomedJet_CA8_GenWhad_pr_c[i]->cd();
   DE_GenGroomedJet_CA8_GenWhad_pr_c[i]->SetLogy();
   DE_GenGroomedJet_CA8_GenWhad_pr[i]->Draw();
   DE_GenGroomedJet_CA8_GenWhad_pr[i]->Write();
   DE_GenGroomedJet_CA8_GenWhad_pr_c[i]->Print((OutputPlotDir+"/"+DE_GenGroomedJet_CA8_GenWhad_pr[i]->GetName()+".png").c_str(),"png");

   DMass_GenGroomedJet_CA8_GenWhad_pr_c[i]   = new TCanvas(DMass_GenGroomedJet_CA8_GenWhad_pr[i]->GetName());
   DMass_GenGroomedJet_CA8_GenWhad_pr_c[i]->cd();
   DMass_GenGroomedJet_CA8_GenWhad_pr_c[i]->SetLogy();
   DMass_GenGroomedJet_CA8_GenWhad_pr[i]->Draw();
   DMass_GenGroomedJet_CA8_GenWhad_pr[i]->Write();
   DMass_GenGroomedJet_CA8_GenWhad_pr_c[i]->Print((OutputPlotDir+"/"+DMass_GenGroomedJet_CA8_GenWhad_pr[i]->GetName()+".png").c_str(),"png");

   DPt_DR_GenGroomedJet_CA8_GenWhad_pr_c[i]   = new TCanvas(DPt_DR_GenGroomedJet_CA8_GenWhad_pr[i]->GetName());
   DPt_DR_GenGroomedJet_CA8_GenWhad_pr_c[i]->cd();
   DPt_DR_GenGroomedJet_CA8_GenWhad_pr_c[i]->SetLogy();
   DPt_DR_GenGroomedJet_CA8_GenWhad_pr[i]->Draw("colz");
   DPt_DR_GenGroomedJet_CA8_GenWhad_pr[i]->Write();
   DPt_DR_GenGroomedJet_CA8_GenWhad_pr_c[i]->Print((OutputPlotDir+"/"+DPt_DR_GenGroomedJet_CA8_GenWhad_pr[i]->GetName()+".png").c_str(),"png");

   DMass_DR_GenGroomedJet_CA8_GenWhad_pr_c[i]   = new TCanvas(DMass_DR_GenGroomedJet_CA8_GenWhad_pr[i]->GetName());
   DMass_DR_GenGroomedJet_CA8_GenWhad_pr_c[i]->cd();
   DMass_DR_GenGroomedJet_CA8_GenWhad_pr[i]->Draw("colz");
   DMass_DR_GenGroomedJet_CA8_GenWhad_pr[i]->Write();
   DMass_DR_GenGroomedJet_CA8_GenWhad_pr_c[i]->Print((OutputPlotDir+"/"+DMass_DR_GenGroomedJet_CA8_GenWhad_pr[i]->GetName()+".png").c_str(),"png");

   DPt_GenGroomedJet_CA8_GenWhad_tr_c[i]  = new TCanvas(DPt_GenGroomedJet_CA8_GenWhad_tr[i]->GetName());
   DPt_GenGroomedJet_CA8_GenWhad_tr_c[i]->cd();
   DPt_GenGroomedJet_CA8_GenWhad_tr_c[i]->SetLogy();
   DPt_GenGroomedJet_CA8_GenWhad_tr[i]->Draw();
   DPt_GenGroomedJet_CA8_GenWhad_tr[i]->Write();
   DPt_GenGroomedJet_CA8_GenWhad_tr_c[i]->Print((OutputPlotDir+"/"+DPt_GenGroomedJet_CA8_GenWhad_tr[i]->GetName()+".png").c_str(),"png");

   DEta_GenGroomedJet_CA8_GenWhad_tr_c[i] = new TCanvas(DEta_GenGroomedJet_CA8_GenWhad_tr[i]->GetName());
   DEta_GenGroomedJet_CA8_GenWhad_tr_c[i]->cd();
   DEta_GenGroomedJet_CA8_GenWhad_tr[i]->Draw();
   DEta_GenGroomedJet_CA8_GenWhad_tr[i]->Write();
   DEta_GenGroomedJet_CA8_GenWhad_tr_c[i]->Print((OutputPlotDir+"/"+DEta_GenGroomedJet_CA8_GenWhad_tr[i]->GetName()+".png").c_str(),"png");

   DPhi_GenGroomedJet_CA8_GenWhad_tr_c[i] = new TCanvas(DPhi_GenGroomedJet_CA8_GenWhad_tr[i]->GetName());
   DPhi_GenGroomedJet_CA8_GenWhad_tr_c[i]->cd();
   DPhi_GenGroomedJet_CA8_GenWhad_tr[i]->Draw();
   DPhi_GenGroomedJet_CA8_GenWhad_tr[i]->Write();
   DPhi_GenGroomedJet_CA8_GenWhad_tr_c[i]->Print((OutputPlotDir+"/"+DPhi_GenGroomedJet_CA8_GenWhad_tr[i]->GetName()+".png").c_str(),"png");

   DR_GenGroomedJet_CA8_GenWhad_tr_c[i]   = new TCanvas(DR_GenGroomedJet_CA8_GenWhad_tr[i]->GetName());
   DR_GenGroomedJet_CA8_GenWhad_tr_c[i]->cd();
   DR_GenGroomedJet_CA8_GenWhad_tr_c[i]->SetLogy();
   DR_GenGroomedJet_CA8_GenWhad_tr[i]->Draw();
   DR_GenGroomedJet_CA8_GenWhad_tr[i]->Write();
   DR_GenGroomedJet_CA8_GenWhad_tr_c[i]->Print((OutputPlotDir+"/"+DR_GenGroomedJet_CA8_GenWhad_tr[i]->GetName()+".png").c_str(),"png");

   DE_GenGroomedJet_CA8_GenWhad_tr_c[i]   = new TCanvas(DE_GenGroomedJet_CA8_GenWhad_tr[i]->GetName());
   DE_GenGroomedJet_CA8_GenWhad_tr_c[i]->cd();
   DE_GenGroomedJet_CA8_GenWhad_tr_c[i]->SetLogy();
   DE_GenGroomedJet_CA8_GenWhad_tr[i]->Draw();
   DE_GenGroomedJet_CA8_GenWhad_tr[i]->Write();
   DE_GenGroomedJet_CA8_GenWhad_tr_c[i]->Print((OutputPlotDir+"/"+DE_GenGroomedJet_CA8_GenWhad_tr[i]->GetName()+".png").c_str(),"png");

   DMass_GenGroomedJet_CA8_GenWhad_tr_c[i]   = new TCanvas(DMass_GenGroomedJet_CA8_GenWhad_tr[i]->GetName());
   DMass_GenGroomedJet_CA8_GenWhad_tr_c[i]->cd();
   DMass_GenGroomedJet_CA8_GenWhad_tr_c[i]->SetLogy();
   DMass_GenGroomedJet_CA8_GenWhad_tr[i]->Draw();
   DMass_GenGroomedJet_CA8_GenWhad_tr[i]->Write();
   DMass_GenGroomedJet_CA8_GenWhad_tr_c[i]->Print((OutputPlotDir+"/"+DMass_GenGroomedJet_CA8_GenWhad_tr[i]->GetName()+".png").c_str(),"png");

   DPt_DR_GenGroomedJet_CA8_GenWhad_tr_c[i]   = new TCanvas(DPt_DR_GenGroomedJet_CA8_GenWhad_tr[i]->GetName());
   DPt_DR_GenGroomedJet_CA8_GenWhad_tr_c[i]->cd();
   DPt_DR_GenGroomedJet_CA8_GenWhad_tr[i]->Draw("colz");
   DPt_DR_GenGroomedJet_CA8_GenWhad_tr[i]->Write();
   DPt_DR_GenGroomedJet_CA8_GenWhad_tr_c[i]->Print((OutputPlotDir+"/"+DPt_DR_GenGroomedJet_CA8_GenWhad_tr[i]->GetName()+".png").c_str(),"png");

   DMass_DR_GenGroomedJet_CA8_GenWhad_tr_c[i]   = new TCanvas(DMass_DR_GenGroomedJet_CA8_GenWhad_tr[i]->GetName());
   DMass_DR_GenGroomedJet_CA8_GenWhad_tr_c[i]->cd();
   DMass_DR_GenGroomedJet_CA8_GenWhad_tr[i]->Draw("colz");
   DMass_DR_GenGroomedJet_CA8_GenWhad_tr[i]->Write();
   DMass_DR_GenGroomedJet_CA8_GenWhad_tr_c[i]->Print((OutputPlotDir+"/"+DMass_DR_GenGroomedJet_CA8_GenWhad_tr[i]->GetName()+".png").c_str(),"png");


   DPt_GenGroomedJet_CA8_GenWhad_ft_c[i]  = new TCanvas(DPt_GenGroomedJet_CA8_GenWhad_ft[i]->GetName());
   DPt_GenGroomedJet_CA8_GenWhad_ft_c[i]->cd();
   DPt_GenGroomedJet_CA8_GenWhad_ft_c[i]->SetLogy();
   DPt_GenGroomedJet_CA8_GenWhad_ft[i]->Draw();
   DPt_GenGroomedJet_CA8_GenWhad_ft[i]->Write();
   DPt_GenGroomedJet_CA8_GenWhad_ft_c[i]->Print((OutputPlotDir+"/"+DPt_GenGroomedJet_CA8_GenWhad_ft[i]->GetName()+".png").c_str(),"png");

   DEta_GenGroomedJet_CA8_GenWhad_ft_c[i] = new TCanvas(DEta_GenGroomedJet_CA8_GenWhad_ft[i]->GetName());
   DEta_GenGroomedJet_CA8_GenWhad_ft_c[i]->cd();
   DEta_GenGroomedJet_CA8_GenWhad_ft[i]->Draw();
   DEta_GenGroomedJet_CA8_GenWhad_ft[i]->Write();
   DEta_GenGroomedJet_CA8_GenWhad_ft_c[i]->Print((OutputPlotDir+"/"+DEta_GenGroomedJet_CA8_GenWhad_ft[i]->GetName()+".png").c_str(),"png");

   DPhi_GenGroomedJet_CA8_GenWhad_ft_c[i] = new TCanvas(DPhi_GenGroomedJet_CA8_GenWhad_ft[i]->GetName());
   DPhi_GenGroomedJet_CA8_GenWhad_ft_c[i]->cd();
   DPhi_GenGroomedJet_CA8_GenWhad_ft[i]->Draw();
   DPhi_GenGroomedJet_CA8_GenWhad_ft[i]->Write();
   DPhi_GenGroomedJet_CA8_GenWhad_ft_c[i]->Print((OutputPlotDir+"/"+DPhi_GenGroomedJet_CA8_GenWhad_ft[i]->GetName()+".png").c_str(),"png");

   DR_GenGroomedJet_CA8_GenWhad_ft_c[i]   = new TCanvas(DR_GenGroomedJet_CA8_GenWhad_ft[i]->GetName());
   DR_GenGroomedJet_CA8_GenWhad_ft_c[i]->cd();
   DR_GenGroomedJet_CA8_GenWhad_ft_c[i]->SetLogy();
   DR_GenGroomedJet_CA8_GenWhad_ft[i]->Draw();
   DR_GenGroomedJet_CA8_GenWhad_ft[i]->Write();
   DR_GenGroomedJet_CA8_GenWhad_ft_c[i]->Print((OutputPlotDir+"/"+DR_GenGroomedJet_CA8_GenWhad_ft[i]->GetName()+".png").c_str(),"png");

   DE_GenGroomedJet_CA8_GenWhad_ft_c[i]   = new TCanvas(DE_GenGroomedJet_CA8_GenWhad_ft[i]->GetName());
   DE_GenGroomedJet_CA8_GenWhad_ft_c[i]->cd();
   DE_GenGroomedJet_CA8_GenWhad_ft_c[i]->SetLogy();
   DE_GenGroomedJet_CA8_GenWhad_ft[i]->Draw();
   DE_GenGroomedJet_CA8_GenWhad_ft[i]->Write();
   DE_GenGroomedJet_CA8_GenWhad_ft_c[i]->Print((OutputPlotDir+"/"+DE_GenGroomedJet_CA8_GenWhad_ft[i]->GetName()+".png").c_str(),"png");

   DMass_GenGroomedJet_CA8_GenWhad_ft_c[i]   = new TCanvas(DMass_GenGroomedJet_CA8_GenWhad_ft[i]->GetName());
   DMass_GenGroomedJet_CA8_GenWhad_ft_c[i]->cd();
   DMass_GenGroomedJet_CA8_GenWhad_ft_c[i]->SetLogy();
   DMass_GenGroomedJet_CA8_GenWhad_ft[i]->Draw();
   DMass_GenGroomedJet_CA8_GenWhad_ft[i]->Write();
   DMass_GenGroomedJet_CA8_GenWhad_ft_c[i]->Print((OutputPlotDir+"/"+DMass_GenGroomedJet_CA8_GenWhad_ft[i]->GetName()+".png").c_str(),"png");

   DPt_DR_GenGroomedJet_CA8_GenWhad_ft_c[i]   = new TCanvas(DPt_DR_GenGroomedJet_CA8_GenWhad_ft[i]->GetName());
   DPt_DR_GenGroomedJet_CA8_GenWhad_ft_c[i]->cd();
   DPt_DR_GenGroomedJet_CA8_GenWhad_ft[i]->Draw("colz");
   DPt_DR_GenGroomedJet_CA8_GenWhad_ft[i]->Write();
   DPt_DR_GenGroomedJet_CA8_GenWhad_ft_c[i]->Print((OutputPlotDir+"/"+DPt_DR_GenGroomedJet_CA8_GenWhad_ft[i]->GetName()+".png").c_str(),"png");

   DMass_DR_GenGroomedJet_CA8_GenWhad_ft_c[i]   = new TCanvas(DMass_DR_GenGroomedJet_CA8_GenWhad_ft[i]->GetName());
   DMass_DR_GenGroomedJet_CA8_GenWhad_ft_c[i]->cd();
   DMass_DR_GenGroomedJet_CA8_GenWhad_ft[i]->Draw("colz");
   DMass_DR_GenGroomedJet_CA8_GenWhad_ft[i]->Write();
   DMass_DR_GenGroomedJet_CA8_GenWhad_ft_c[i]->Print((OutputPlotDir+"/"+DMass_DR_GenGroomedJet_CA8_GenWhad_ft[i]->GetName()+".png").c_str(),"png");

   // Reco Groomed Jet AK5 

   DPt_GroomedJet_AK5_GenWhad_c[i]  = new TCanvas(DPt_GroomedJet_AK5_GenWhad[i]->GetName());
   DPt_GroomedJet_AK5_GenWhad_c[i]->cd();
   DPt_GroomedJet_AK5_GenWhad_c[i]->SetLogy();
   DPt_GroomedJet_AK5_GenWhad[i]->Draw();
   DPt_GroomedJet_AK5_GenWhad[i]->Write();
   DPt_GroomedJet_AK5_GenWhad_c[i]->Print((OutputPlotDir+"/"+DPt_GroomedJet_AK5_GenWhad[i]->GetName()+".png").c_str(),"png");

   DEta_GroomedJet_AK5_GenWhad_c[i] = new TCanvas(DEta_GroomedJet_AK5_GenWhad[i]->GetName());
   DEta_GroomedJet_AK5_GenWhad_c[i]->cd();
   DEta_GroomedJet_AK5_GenWhad[i]->Draw();
   DEta_GroomedJet_AK5_GenWhad[i]->Write();
   DEta_GroomedJet_AK5_GenWhad_c[i]->Print((OutputPlotDir+"/"+DEta_GroomedJet_AK5_GenWhad[i]->GetName()+".png").c_str(),"png");

   DPhi_GroomedJet_AK5_GenWhad_c[i] = new TCanvas(DPhi_GroomedJet_AK5_GenWhad[i]->GetName());
   DPhi_GroomedJet_AK5_GenWhad_c[i]->cd();
   DPhi_GroomedJet_AK5_GenWhad[i]->Draw();
   DPhi_GroomedJet_AK5_GenWhad[i]->Write();
   DPhi_GroomedJet_AK5_GenWhad_c[i]->Print((OutputPlotDir+"/"+DPhi_GroomedJet_AK5_GenWhad[i]->GetName()+".png").c_str(),"png");

   DR_GroomedJet_AK5_GenWhad_c[i]   = new TCanvas(DR_GroomedJet_AK5_GenWhad[i]->GetName());
   DR_GroomedJet_AK5_GenWhad_c[i]->cd();
   DR_GroomedJet_AK5_GenWhad_c[i]->SetLogy();
   DR_GroomedJet_AK5_GenWhad[i]->Draw();
   DR_GroomedJet_AK5_GenWhad[i]->Write();
   DR_GroomedJet_AK5_GenWhad_c[i]->Print((OutputPlotDir+"/"+DR_GroomedJet_AK5_GenWhad[i]->GetName()+".png").c_str(),"png");

   DE_GroomedJet_AK5_GenWhad_c[i]   = new TCanvas(DE_GroomedJet_AK5_GenWhad[i]->GetName());
   DE_GroomedJet_AK5_GenWhad_c[i]->cd();
   DE_GroomedJet_AK5_GenWhad_c[i]->SetLogy();
   DE_GroomedJet_AK5_GenWhad[i]->Draw();
   DE_GroomedJet_AK5_GenWhad[i]->Write();
   DE_GroomedJet_AK5_GenWhad_c[i]->Print((OutputPlotDir+"/"+DE_GroomedJet_AK5_GenWhad[i]->GetName()+".png").c_str(),"png");

   DMass_GroomedJet_AK5_GenWhad_c[i]   = new TCanvas(DMass_GroomedJet_AK5_GenWhad[i]->GetName());
   DMass_GroomedJet_AK5_GenWhad_c[i]->cd();
   DMass_GroomedJet_AK5_GenWhad_c[i]->SetLogy();
   DMass_GroomedJet_AK5_GenWhad[i]->Draw();
   DMass_GroomedJet_AK5_GenWhad[i]->Write();
   DMass_GroomedJet_AK5_GenWhad_c[i]->Print((OutputPlotDir+"/"+DMass_GroomedJet_AK5_GenWhad[i]->GetName()+".png").c_str(),"png");

   DPt_DR_GroomedJet_AK5_GenWhad_c[i]   = new TCanvas(DPt_DR_GroomedJet_AK5_GenWhad[i]->GetName());
   DPt_DR_GroomedJet_AK5_GenWhad_c[i]->cd();
   DPt_DR_GroomedJet_AK5_GenWhad[i]->Draw("colz");
   DPt_DR_GroomedJet_AK5_GenWhad[i]->Write();
   DPt_DR_GroomedJet_AK5_GenWhad_c[i]->Print((OutputPlotDir+"/"+DPt_DR_GroomedJet_AK5_GenWhad[i]->GetName()+".png").c_str(),"png");

   DMass_DR_GroomedJet_AK5_GenWhad_c[i]   = new TCanvas(DMass_DR_GroomedJet_AK5_GenWhad[i]->GetName());
   DMass_DR_GroomedJet_AK5_GenWhad_c[i]->cd();
   DMass_DR_GroomedJet_AK5_GenWhad[i]->Draw("colz");
   DMass_DR_GroomedJet_AK5_GenWhad[i]->Write();
   DMass_DR_GroomedJet_AK5_GenWhad_c[i]->Print((OutputPlotDir+"/"+DMass_DR_GroomedJet_AK5_GenWhad[i]->GetName()+".png").c_str(),"png");


   // Reco Groomed Jet CA8 

   DPt_GroomedJet_CA8_GenWhad_c[i]  = new TCanvas(DPt_GroomedJet_CA8_GenWhad[i]->GetName());
   DPt_GroomedJet_CA8_GenWhad_c[i]->cd();
   DPt_GroomedJet_CA8_GenWhad_c[i]->SetLogy();
   DPt_GroomedJet_CA8_GenWhad[i]->Draw();
   DPt_GroomedJet_CA8_GenWhad[i]->Write();
   DPt_GroomedJet_CA8_GenWhad_c[i]->Print((OutputPlotDir+"/"+DPt_GroomedJet_CA8_GenWhad[i]->GetName()+".png").c_str(),"png");

   DEta_GroomedJet_CA8_GenWhad_c[i] = new TCanvas(DEta_GroomedJet_CA8_GenWhad[i]->GetName());
   DEta_GroomedJet_CA8_GenWhad_c[i]->cd();
   DEta_GroomedJet_CA8_GenWhad[i]->Draw();
   DEta_GroomedJet_CA8_GenWhad[i]->Write();
   DEta_GroomedJet_CA8_GenWhad_c[i]->Print((OutputPlotDir+"/"+DEta_GroomedJet_CA8_GenWhad[i]->GetName()+".png").c_str(),"png");

   DPhi_GroomedJet_CA8_GenWhad_c[i] = new TCanvas(DPhi_GroomedJet_CA8_GenWhad[i]->GetName());
   DPhi_GroomedJet_CA8_GenWhad_c[i]->cd();
   DPhi_GroomedJet_CA8_GenWhad[i]->Draw();
   DPhi_GroomedJet_CA8_GenWhad[i]->Write();
   DPhi_GroomedJet_CA8_GenWhad_c[i]->Print((OutputPlotDir+"/"+DPhi_GroomedJet_CA8_GenWhad[i]->GetName()+".png").c_str(),"png");

   DR_GroomedJet_CA8_GenWhad_c[i]   = new TCanvas(DR_GroomedJet_CA8_GenWhad[i]->GetName());
   DR_GroomedJet_CA8_GenWhad_c[i]->cd();
   DR_GroomedJet_CA8_GenWhad_c[i]->SetLogy();
   DR_GroomedJet_CA8_GenWhad[i]->Draw();
   DR_GroomedJet_CA8_GenWhad[i]->Write();
   DR_GroomedJet_CA8_GenWhad_c[i]->Print((OutputPlotDir+"/"+DR_GroomedJet_CA8_GenWhad[i]->GetName()+".png").c_str(),"png");

   DE_GroomedJet_CA8_GenWhad_c[i]   = new TCanvas(DE_GroomedJet_CA8_GenWhad[i]->GetName());
   DE_GroomedJet_CA8_GenWhad_c[i]->cd();
   DE_GroomedJet_CA8_GenWhad_c[i]->SetLogy();
   DE_GroomedJet_CA8_GenWhad[i]->Draw();
   DE_GroomedJet_CA8_GenWhad[i]->Write();
   DE_GroomedJet_CA8_GenWhad_c[i]->Print((OutputPlotDir+"/"+DE_GroomedJet_CA8_GenWhad[i]->GetName()+".png").c_str(),"png");

   DMass_GroomedJet_CA8_GenWhad_c[i]   = new TCanvas(DMass_GroomedJet_CA8_GenWhad[i]->GetName());
   DMass_GroomedJet_CA8_GenWhad_c[i]->cd();
   DMass_GroomedJet_CA8_GenWhad_c[i]->SetLogy();
   DMass_GroomedJet_CA8_GenWhad[i]->Draw();
   DMass_GroomedJet_CA8_GenWhad[i]->Write();
   DMass_GroomedJet_CA8_GenWhad_c[i]->Print((OutputPlotDir+"/"+DMass_GroomedJet_CA8_GenWhad[i]->GetName()+".png").c_str(),"png");

   DPt_DR_GroomedJet_CA8_GenWhad_c[i]   = new TCanvas(DPt_DR_GroomedJet_CA8_GenWhad[i]->GetName());
   DPt_DR_GroomedJet_CA8_GenWhad_c[i]->cd();
   DPt_DR_GroomedJet_CA8_GenWhad[i]->Draw("colz");
   DPt_DR_GroomedJet_CA8_GenWhad[i]->Write();
   DPt_DR_GroomedJet_CA8_GenWhad_c[i]->Print((OutputPlotDir+"/"+DPt_DR_GroomedJet_CA8_GenWhad[i]->GetName()+".png").c_str(),"png");

   DMass_DR_GroomedJet_CA8_GenWhad_c[i]   = new TCanvas(DMass_DR_GroomedJet_CA8_GenWhad[i]->GetName());
   DMass_DR_GroomedJet_CA8_GenWhad_c[i]->cd();
   DMass_DR_GroomedJet_CA8_GenWhad[i]->Draw("colz");
   DMass_DR_GroomedJet_CA8_GenWhad[i]->Write();
   DMass_DR_GroomedJet_CA8_GenWhad_c[i]->Print((OutputPlotDir+"/"+DMass_DR_GroomedJet_CA8_GenWhad[i]->GetName()+".png").c_str(),"png");

   DPt_GroomedJet_CA8_subjet1_WParton_c[i]  = new TCanvas(DPt_GroomedJet_CA8_subjet1_WParton[i]->GetName());
   DPt_GroomedJet_CA8_subjet1_WParton_c[i]->cd();
   DPt_GroomedJet_CA8_subjet1_WParton_c[i]->SetLogy();
   DPt_GroomedJet_CA8_subjet1_WParton[i]->Draw();
   DPt_GroomedJet_CA8_subjet1_WParton[i]->Write();
   DPt_GroomedJet_CA8_subjet1_WParton_c[i]->Print((OutputPlotDir+"/"+DPt_GroomedJet_CA8_subjet1_WParton[i]->GetName()+".png").c_str(),"png");

   DEta_GroomedJet_CA8_subjet1_WParton_c[i] = new TCanvas(DEta_GroomedJet_CA8_subjet1_WParton[i]->GetName());
   DEta_GroomedJet_CA8_subjet1_WParton_c[i]->cd();
   DEta_GroomedJet_CA8_subjet1_WParton[i]->Draw();
   DEta_GroomedJet_CA8_subjet1_WParton[i]->Write();
   DEta_GroomedJet_CA8_subjet1_WParton_c[i]->Print((OutputPlotDir+"/"+DEta_GroomedJet_CA8_subjet1_WParton[i]->GetName()+".png").c_str(),"png");

   DPhi_GroomedJet_CA8_subjet1_WParton_c[i] = new TCanvas(DPhi_GroomedJet_CA8_subjet1_WParton[i]->GetName());
   DPhi_GroomedJet_CA8_subjet1_WParton_c[i]->cd();
   DPhi_GroomedJet_CA8_subjet1_WParton[i]->Draw();
   DPhi_GroomedJet_CA8_subjet1_WParton[i]->Write();
   DPhi_GroomedJet_CA8_subjet1_WParton_c[i]->Print((OutputPlotDir+"/"+DPhi_GroomedJet_CA8_subjet1_WParton[i]->GetName()+".png").c_str(),"png");

   DR_GroomedJet_CA8_subjet1_WParton_c[i]   = new TCanvas(DR_GroomedJet_CA8_subjet1_WParton[i]->GetName());
   DR_GroomedJet_CA8_subjet1_WParton_c[i]->cd();
   DR_GroomedJet_CA8_subjet1_WParton_c[i]->SetLogy();
   DR_GroomedJet_CA8_subjet1_WParton[i]->Draw();
   DR_GroomedJet_CA8_subjet1_WParton[i]->Write();
   DR_GroomedJet_CA8_subjet1_WParton_c[i]->Print((OutputPlotDir+"/"+DR_GroomedJet_CA8_subjet1_WParton[i]->GetName()+".png").c_str(),"png");

   DE_GroomedJet_CA8_subjet1_WParton_c[i]   = new TCanvas(DE_GroomedJet_CA8_subjet1_WParton[i]->GetName());
   DE_GroomedJet_CA8_subjet1_WParton_c[i]->cd();
   DE_GroomedJet_CA8_subjet1_WParton_c[i]->SetLogy();
   DE_GroomedJet_CA8_subjet1_WParton[i]->Draw();
   DE_GroomedJet_CA8_subjet1_WParton[i]->Write();
   DE_GroomedJet_CA8_subjet1_WParton_c[i]->Print((OutputPlotDir+"/"+DE_GroomedJet_CA8_subjet1_WParton[i]->GetName()+".png").c_str(),"png");

   DMass_GroomedJet_CA8_subjet1_WParton_c[i]   = new TCanvas(DMass_GroomedJet_CA8_subjet1_WParton[i]->GetName());
   DMass_GroomedJet_CA8_subjet1_WParton_c[i]->cd();
   DMass_GroomedJet_CA8_subjet1_WParton_c[i]->SetLogy();
   DMass_GroomedJet_CA8_subjet1_WParton[i]->Draw();
   DMass_GroomedJet_CA8_subjet1_WParton[i]->Write();
   DMass_GroomedJet_CA8_subjet1_WParton_c[i]->Print((OutputPlotDir+"/"+DMass_GroomedJet_CA8_subjet1_WParton[i]->GetName()+".png").c_str(),"png");

   DPt_DR_GroomedJet_CA8_subjet1_WParton_c[i]   = new TCanvas(DPt_DR_GroomedJet_CA8_subjet1_WParton[i]->GetName());
   DPt_DR_GroomedJet_CA8_subjet1_WParton_c[i]->cd();
   DPt_DR_GroomedJet_CA8_subjet1_WParton[i]->Draw("colz");
   DPt_DR_GroomedJet_CA8_subjet1_WParton[i]->Write();
   DPt_DR_GroomedJet_CA8_subjet1_WParton_c[i]->Print((OutputPlotDir+"/"+DPt_DR_GroomedJet_CA8_subjet1_WParton[i]->GetName()+".png").c_str(),"png");

   DMass_DR_GroomedJet_CA8_subjet1_WParton_c[i]   = new TCanvas(DMass_DR_GroomedJet_CA8_subjet1_WParton[i]->GetName());
   DMass_DR_GroomedJet_CA8_subjet1_WParton_c[i]->cd();
   DMass_DR_GroomedJet_CA8_subjet1_WParton[i]->Draw("colz");
   DMass_DR_GroomedJet_CA8_subjet1_WParton[i]->Write();
   DMass_DR_GroomedJet_CA8_subjet1_WParton_c[i]->Print((OutputPlotDir+"/"+DMass_DR_GroomedJet_CA8_subjet1_WParton[i]->GetName()+".png").c_str(),"png");

   DPt_GroomedJet_CA8_subjet2_WParton_c[i]  = new TCanvas(DPt_GroomedJet_CA8_subjet2_WParton[i]->GetName());
   DPt_GroomedJet_CA8_subjet2_WParton_c[i]->cd();
   DPt_GroomedJet_CA8_subjet2_WParton_c[i]->SetLogy();
   DPt_GroomedJet_CA8_subjet2_WParton[i]->Draw();
   DPt_GroomedJet_CA8_subjet2_WParton[i]->Write();
   DPt_GroomedJet_CA8_subjet2_WParton_c[i]->Print((OutputPlotDir+"/"+DPt_GroomedJet_CA8_subjet2_WParton[i]->GetName()+".png").c_str(),"png");

   DEta_GroomedJet_CA8_subjet2_WParton_c[i] = new TCanvas(DEta_GroomedJet_CA8_subjet2_WParton[i]->GetName());
   DEta_GroomedJet_CA8_subjet2_WParton_c[i]->cd();
   DEta_GroomedJet_CA8_subjet2_WParton[i]->Draw();
   DEta_GroomedJet_CA8_subjet2_WParton[i]->Write();
   DEta_GroomedJet_CA8_subjet2_WParton_c[i]->Print((OutputPlotDir+"/"+DEta_GroomedJet_CA8_subjet2_WParton[i]->GetName()+".png").c_str(),"png");

   DPhi_GroomedJet_CA8_subjet2_WParton_c[i] = new TCanvas(DPhi_GroomedJet_CA8_subjet2_WParton[i]->GetName());
   DPhi_GroomedJet_CA8_subjet2_WParton_c[i]->cd();
   DPhi_GroomedJet_CA8_subjet2_WParton[i]->Draw();
   DPhi_GroomedJet_CA8_subjet2_WParton[i]->Write();
   DPhi_GroomedJet_CA8_subjet2_WParton_c[i]->Print((OutputPlotDir+"/"+DPhi_GroomedJet_CA8_subjet2_WParton[i]->GetName()+".png").c_str(),"png");

   DR_GroomedJet_CA8_subjet2_WParton_c[i]   = new TCanvas(DR_GroomedJet_CA8_subjet2_WParton[i]->GetName());
   DR_GroomedJet_CA8_subjet2_WParton_c[i]->cd();
   DR_GroomedJet_CA8_subjet2_WParton_c[i]->SetLogy();
   DR_GroomedJet_CA8_subjet2_WParton[i]->Draw();
   DR_GroomedJet_CA8_subjet2_WParton[i]->Write();
   DR_GroomedJet_CA8_subjet2_WParton_c[i]->Print((OutputPlotDir+"/"+DR_GroomedJet_CA8_subjet2_WParton[i]->GetName()+".png").c_str(),"png");

   DE_GroomedJet_CA8_subjet2_WParton_c[i]   = new TCanvas(DE_GroomedJet_CA8_subjet2_WParton[i]->GetName());
   DE_GroomedJet_CA8_subjet2_WParton_c[i]->cd();
   DE_GroomedJet_CA8_subjet2_WParton_c[i]->SetLogy();
   DE_GroomedJet_CA8_subjet2_WParton[i]->Draw();
   DE_GroomedJet_CA8_subjet2_WParton[i]->Write();
   DE_GroomedJet_CA8_subjet2_WParton_c[i]->Print((OutputPlotDir+"/"+DE_GroomedJet_CA8_subjet2_WParton[i]->GetName()+".png").c_str(),"png");

   DMass_GroomedJet_CA8_subjet2_WParton_c[i]   = new TCanvas(DMass_GroomedJet_CA8_subjet2_WParton[i]->GetName());
   DMass_GroomedJet_CA8_subjet2_WParton_c[i]->cd();
   DMass_GroomedJet_CA8_subjet2_WParton_c[i]->SetLogy();
   DMass_GroomedJet_CA8_subjet2_WParton[i]->Draw();
   DMass_GroomedJet_CA8_subjet2_WParton[i]->Write();
   DMass_GroomedJet_CA8_subjet2_WParton_c[i]->Print((OutputPlotDir+"/"+DMass_GroomedJet_CA8_subjet2_WParton[i]->GetName()+".png").c_str(),"png");

   DPt_DR_GroomedJet_CA8_subjet2_WParton_c[i]   = new TCanvas(DPt_DR_GroomedJet_CA8_subjet2_WParton[i]->GetName());
   DPt_DR_GroomedJet_CA8_subjet2_WParton_c[i]->cd();
   DPt_DR_GroomedJet_CA8_subjet2_WParton[i]->Draw("colz");
   DPt_DR_GroomedJet_CA8_subjet2_WParton[i]->Write();
   DPt_DR_GroomedJet_CA8_subjet2_WParton_c[i]->Print((OutputPlotDir+"/"+DPt_DR_GroomedJet_CA8_subjet2_WParton[i]->GetName()+".png").c_str(),"png");

   DMass_DR_GroomedJet_CA8_subjet2_WParton_c[i]   = new TCanvas(DMass_DR_GroomedJet_CA8_subjet2_WParton[i]->GetName());
   DMass_DR_GroomedJet_CA8_subjet2_WParton_c[i]->cd();
   DMass_DR_GroomedJet_CA8_subjet2_WParton[i]->Draw("colz");
   DMass_DR_GroomedJet_CA8_subjet2_WParton[i]->Write();
   DMass_DR_GroomedJet_CA8_subjet2_WParton_c[i]->Print((OutputPlotDir+"/"+DMass_DR_GroomedJet_CA8_subjet2_WParton[i]->GetName()+".png").c_str(),"png");

   DPt_GroomedJet_CA8_GenWhad_pr_c[i]  = new TCanvas(DPt_GroomedJet_CA8_GenWhad_pr[i]->GetName());
   DPt_GroomedJet_CA8_GenWhad_pr_c[i]->cd();
   DPt_GroomedJet_CA8_GenWhad_pr_c[i]->SetLogy();
   DPt_GroomedJet_CA8_GenWhad_pr[i]->Draw();
   DPt_GroomedJet_CA8_GenWhad_pr[i]->Write();
   DPt_GroomedJet_CA8_GenWhad_pr_c[i]->Print((OutputPlotDir+"/"+DPt_GroomedJet_CA8_GenWhad_pr[i]->GetName()+".png").c_str(),"png");

   DEta_GroomedJet_CA8_GenWhad_pr_c[i] = new TCanvas(DEta_GroomedJet_CA8_GenWhad_pr[i]->GetName());
   DEta_GroomedJet_CA8_GenWhad_pr_c[i]->cd();
   DEta_GroomedJet_CA8_GenWhad_pr[i]->Draw();
   DEta_GroomedJet_CA8_GenWhad_pr[i]->Write();
   DEta_GroomedJet_CA8_GenWhad_pr_c[i]->Print((OutputPlotDir+"/"+DEta_GroomedJet_CA8_GenWhad_pr[i]->GetName()+".png").c_str(),"png");

   DPhi_GroomedJet_CA8_GenWhad_pr_c[i] = new TCanvas(DPhi_GroomedJet_CA8_GenWhad_pr[i]->GetName());
   DPhi_GroomedJet_CA8_GenWhad_pr_c[i]->cd();
   DPhi_GroomedJet_CA8_GenWhad_pr[i]->Draw();
   DPhi_GroomedJet_CA8_GenWhad_pr[i]->Write();
   DPhi_GroomedJet_CA8_GenWhad_pr_c[i]->Print((OutputPlotDir+"/"+DPhi_GroomedJet_CA8_GenWhad_pr[i]->GetName()+".png").c_str(),"png");

   DR_GroomedJet_CA8_GenWhad_pr_c[i]   = new TCanvas(DR_GroomedJet_CA8_GenWhad_pr[i]->GetName());
   DR_GroomedJet_CA8_GenWhad_pr_c[i]->cd();
   DR_GroomedJet_CA8_GenWhad_pr_c[i]->SetLogy();
   DR_GroomedJet_CA8_GenWhad_pr[i]->Draw();
   DR_GroomedJet_CA8_GenWhad_pr[i]->Write();
   DR_GroomedJet_CA8_GenWhad_pr_c[i]->Print((OutputPlotDir+"/"+DR_GroomedJet_CA8_GenWhad_pr[i]->GetName()+".png").c_str(),"png");

   DE_GroomedJet_CA8_GenWhad_pr_c[i]   = new TCanvas(DE_GroomedJet_CA8_GenWhad_pr[i]->GetName());
   DE_GroomedJet_CA8_GenWhad_pr_c[i]->cd();
   DE_GroomedJet_CA8_GenWhad_pr_c[i]->SetLogy();
   DE_GroomedJet_CA8_GenWhad_pr[i]->Draw();
   DE_GroomedJet_CA8_GenWhad_pr[i]->Write();
   DE_GroomedJet_CA8_GenWhad_pr_c[i]->Print((OutputPlotDir+"/"+DE_GroomedJet_CA8_GenWhad_pr[i]->GetName()+".png").c_str(),"png");

   DMass_GroomedJet_CA8_GenWhad_pr_c[i]   = new TCanvas(DMass_GroomedJet_CA8_GenWhad_pr[i]->GetName());
   DMass_GroomedJet_CA8_GenWhad_pr_c[i]->cd();
   DMass_GroomedJet_CA8_GenWhad_pr_c[i]->SetLogy();
   DMass_GroomedJet_CA8_GenWhad_pr[i]->Draw();
   DMass_GroomedJet_CA8_GenWhad_pr[i]->Write();
   DMass_GroomedJet_CA8_GenWhad_pr_c[i]->Print((OutputPlotDir+"/"+DMass_GroomedJet_CA8_GenWhad_pr[i]->GetName()+".png").c_str(),"png");

   DPt_DR_GroomedJet_CA8_GenWhad_pr_c[i]   = new TCanvas(DPt_DR_GroomedJet_CA8_GenWhad_pr[i]->GetName());
   DPt_DR_GroomedJet_CA8_GenWhad_pr_c[i]->cd();
   DPt_DR_GroomedJet_CA8_GenWhad_pr[i]->Draw("colz");
   DPt_DR_GroomedJet_CA8_GenWhad_pr[i]->Write();
   DPt_DR_GroomedJet_CA8_GenWhad_pr_c[i]->Print((OutputPlotDir+"/"+DPt_DR_GroomedJet_CA8_GenWhad_pr[i]->GetName()+".png").c_str(),"png");

   DMass_DR_GroomedJet_CA8_GenWhad_pr_c[i]   = new TCanvas(DMass_DR_GroomedJet_CA8_GenWhad_pr[i]->GetName());
   DMass_DR_GroomedJet_CA8_GenWhad_pr_c[i]->cd();
   DMass_DR_GroomedJet_CA8_GenWhad_pr[i]->Draw("colz");
   DMass_DR_GroomedJet_CA8_GenWhad_pr[i]->Write();
   DMass_DR_GroomedJet_CA8_GenWhad_pr_c[i]->Print((OutputPlotDir+"/"+DMass_DR_GroomedJet_CA8_GenWhad_pr[i]->GetName()+".png").c_str(),"png");

   DPt_GroomedJet_CA8_GenWhad_tr_c[i]  = new TCanvas(DPt_GroomedJet_CA8_GenWhad_tr[i]->GetName());
   DPt_GroomedJet_CA8_GenWhad_tr_c[i]->cd();
   DPt_GroomedJet_CA8_GenWhad_tr_c[i]->SetLogy();
   DPt_GroomedJet_CA8_GenWhad_tr[i]->Draw();
   DPt_GroomedJet_CA8_GenWhad_tr[i]->Write();
   DPt_GroomedJet_CA8_GenWhad_tr_c[i]->Print((OutputPlotDir+"/"+DPt_GroomedJet_CA8_GenWhad_tr[i]->GetName()+".png").c_str(),"png");

   DEta_GroomedJet_CA8_GenWhad_tr_c[i] = new TCanvas(DEta_GroomedJet_CA8_GenWhad_tr[i]->GetName());
   DEta_GroomedJet_CA8_GenWhad_tr_c[i]->cd();
   DEta_GroomedJet_CA8_GenWhad_tr[i]->Draw();
   DEta_GroomedJet_CA8_GenWhad_tr[i]->Write();
   DEta_GroomedJet_CA8_GenWhad_tr_c[i]->Print((OutputPlotDir+"/"+DEta_GroomedJet_CA8_GenWhad_tr[i]->GetName()+".png").c_str(),"png");

   DPhi_GroomedJet_CA8_GenWhad_tr_c[i] = new TCanvas(DPhi_GroomedJet_CA8_GenWhad_tr[i]->GetName());
   DPhi_GroomedJet_CA8_GenWhad_tr_c[i]->cd();
   DPhi_GroomedJet_CA8_GenWhad_tr[i]->Draw();
   DPhi_GroomedJet_CA8_GenWhad_tr[i]->Write();
   DPhi_GroomedJet_CA8_GenWhad_tr_c[i]->Print((OutputPlotDir+"/"+DPhi_GroomedJet_CA8_GenWhad_tr[i]->GetName()+".png").c_str(),"png");

   DR_GroomedJet_CA8_GenWhad_tr_c[i]   = new TCanvas(DR_GroomedJet_CA8_GenWhad_tr[i]->GetName());
   DR_GroomedJet_CA8_GenWhad_tr_c[i]->cd();
   DR_GroomedJet_CA8_GenWhad_tr_c[i]->SetLogy();
   DR_GroomedJet_CA8_GenWhad_tr[i]->Draw();
   DR_GroomedJet_CA8_GenWhad_tr[i]->Write();
   DR_GroomedJet_CA8_GenWhad_tr_c[i]->Print((OutputPlotDir+"/"+DR_GroomedJet_CA8_GenWhad_tr[i]->GetName()+".png").c_str(),"png");

   DE_GroomedJet_CA8_GenWhad_tr_c[i]   = new TCanvas(DE_GroomedJet_CA8_GenWhad_tr[i]->GetName());
   DE_GroomedJet_CA8_GenWhad_tr_c[i]->cd();
   DE_GroomedJet_CA8_GenWhad_tr_c[i]->SetLogy();
   DE_GroomedJet_CA8_GenWhad_tr[i]->Draw();
   DE_GroomedJet_CA8_GenWhad_tr[i]->Write();
   DE_GroomedJet_CA8_GenWhad_tr_c[i]->Print((OutputPlotDir+"/"+DE_GroomedJet_CA8_GenWhad_tr[i]->GetName()+".png").c_str(),"png");

   DMass_GroomedJet_CA8_GenWhad_tr_c[i]   = new TCanvas(DMass_GroomedJet_CA8_GenWhad_tr[i]->GetName());
   DMass_GroomedJet_CA8_GenWhad_tr_c[i]->cd();
   DMass_GroomedJet_CA8_GenWhad_tr_c[i]->SetLogy();
   DMass_GroomedJet_CA8_GenWhad_tr[i]->Draw();
   DMass_GroomedJet_CA8_GenWhad_tr[i]->Write();
   DMass_GroomedJet_CA8_GenWhad_tr_c[i]->Print((OutputPlotDir+"/"+DMass_GroomedJet_CA8_GenWhad_tr[i]->GetName()+".png").c_str(),"png");

   DPt_DR_GroomedJet_CA8_GenWhad_tr_c[i]   = new TCanvas(DPt_DR_GroomedJet_CA8_GenWhad_tr[i]->GetName());
   DPt_DR_GroomedJet_CA8_GenWhad_tr_c[i]->cd();
   DPt_DR_GroomedJet_CA8_GenWhad_tr[i]->Draw("colz");
   DPt_DR_GroomedJet_CA8_GenWhad_tr[i]->Write();
   DPt_DR_GroomedJet_CA8_GenWhad_tr_c[i]->Print((OutputPlotDir+"/"+DPt_DR_GroomedJet_CA8_GenWhad_tr[i]->GetName()+".png").c_str(),"png");

   DMass_DR_GroomedJet_CA8_GenWhad_tr_c[i]   = new TCanvas(DMass_DR_GroomedJet_CA8_GenWhad_tr[i]->GetName());
   DMass_DR_GroomedJet_CA8_GenWhad_tr_c[i]->cd();
   DMass_DR_GroomedJet_CA8_GenWhad_tr[i]->Draw("colz");
   DMass_DR_GroomedJet_CA8_GenWhad_tr[i]->Write();
   DMass_DR_GroomedJet_CA8_GenWhad_tr_c[i]->Print((OutputPlotDir+"/"+DMass_DR_GroomedJet_CA8_GenWhad_tr[i]->GetName()+".png").c_str(),"png");

   DPt_GroomedJet_CA8_GenWhad_ft_c[i]  = new TCanvas(DPt_GroomedJet_CA8_GenWhad_ft[i]->GetName());
   DPt_GroomedJet_CA8_GenWhad_ft_c[i]->cd();
   DPt_GroomedJet_CA8_GenWhad_ft_c[i]->SetLogy();
   DPt_GroomedJet_CA8_GenWhad_ft[i]->Draw();
   DPt_GroomedJet_CA8_GenWhad_ft[i]->Write();
   DPt_GroomedJet_CA8_GenWhad_ft_c[i]->Print((OutputPlotDir+"/"+DPt_GroomedJet_CA8_GenWhad_ft[i]->GetName()+".png").c_str(),"png");

   DEta_GroomedJet_CA8_GenWhad_ft_c[i] = new TCanvas(DEta_GroomedJet_CA8_GenWhad_ft[i]->GetName());
   DEta_GroomedJet_CA8_GenWhad_ft_c[i]->cd();
   DEta_GroomedJet_CA8_GenWhad_ft[i]->Draw();
   DEta_GroomedJet_CA8_GenWhad_ft[i]->Write();
   DEta_GroomedJet_CA8_GenWhad_ft_c[i]->Print((OutputPlotDir+"/"+DEta_GroomedJet_CA8_GenWhad_ft[i]->GetName()+".png").c_str(),"png");

   DPhi_GroomedJet_CA8_GenWhad_ft_c[i] = new TCanvas(DPhi_GroomedJet_CA8_GenWhad_ft[i]->GetName());
   DPhi_GroomedJet_CA8_GenWhad_ft_c[i]->cd();
   DPhi_GroomedJet_CA8_GenWhad_ft[i]->Draw();
   DPhi_GroomedJet_CA8_GenWhad_ft[i]->Write();
   DPhi_GroomedJet_CA8_GenWhad_ft_c[i]->Print((OutputPlotDir+"/"+DPhi_GroomedJet_CA8_GenWhad_ft[i]->GetName()+".png").c_str(),"png");

   DR_GroomedJet_CA8_GenWhad_ft_c[i]   = new TCanvas(DR_GroomedJet_CA8_GenWhad_ft[i]->GetName());
   DR_GroomedJet_CA8_GenWhad_ft_c[i]->cd();
   DR_GroomedJet_CA8_GenWhad_ft_c[i]->SetLogy();
   DR_GroomedJet_CA8_GenWhad_ft[i]->Draw();
   DR_GroomedJet_CA8_GenWhad_ft[i]->Write();
   DR_GroomedJet_CA8_GenWhad_ft_c[i]->Print((OutputPlotDir+"/"+DR_GroomedJet_CA8_GenWhad_ft[i]->GetName()+".png").c_str(),"png");

   DE_GroomedJet_CA8_GenWhad_ft_c[i]   = new TCanvas(DE_GroomedJet_CA8_GenWhad_ft[i]->GetName());
   DE_GroomedJet_CA8_GenWhad_ft_c[i]->cd();
   DE_GroomedJet_CA8_GenWhad_ft_c[i]->SetLogy();
   DE_GroomedJet_CA8_GenWhad_ft[i]->Draw();
   DE_GroomedJet_CA8_GenWhad_ft[i]->Write();
   DE_GroomedJet_CA8_GenWhad_ft_c[i]->Print((OutputPlotDir+"/"+DE_GroomedJet_CA8_GenWhad_ft[i]->GetName()+".png").c_str(),"png");

   DMass_GroomedJet_CA8_GenWhad_ft_c[i]   = new TCanvas(DMass_GroomedJet_CA8_GenWhad_ft[i]->GetName());
   DMass_GroomedJet_CA8_GenWhad_ft_c[i]->cd();
   DMass_GroomedJet_CA8_GenWhad_ft_c[i]->SetLogy();
   DMass_GroomedJet_CA8_GenWhad_ft[i]->Draw();
   DMass_GroomedJet_CA8_GenWhad_ft[i]->Write();
   DMass_GroomedJet_CA8_GenWhad_ft_c[i]->Print((OutputPlotDir+"/"+DMass_GroomedJet_CA8_GenWhad_ft[i]->GetName()+".png").c_str(),"png");

   DPt_DR_GroomedJet_CA8_GenWhad_ft_c[i]   = new TCanvas(DPt_DR_GroomedJet_CA8_GenWhad_ft[i]->GetName());
   DPt_DR_GroomedJet_CA8_GenWhad_ft_c[i]->cd();
   DPt_DR_GroomedJet_CA8_GenWhad_ft[i]->Draw("colz");
   DPt_DR_GroomedJet_CA8_GenWhad_ft[i]->Write();
   DPt_DR_GroomedJet_CA8_GenWhad_ft_c[i]->Print((OutputPlotDir+"/"+DPt_DR_GroomedJet_CA8_GenWhad_ft[i]->GetName()+".png").c_str(),"png");

   DMass_DR_GroomedJet_CA8_GenWhad_ft_c[i]   = new TCanvas(DMass_DR_GroomedJet_CA8_GenWhad_ft[i]->GetName());
   DMass_DR_GroomedJet_CA8_GenWhad_ft_c[i]->cd();
   DMass_DR_GroomedJet_CA8_GenWhad_ft[i]->Draw("colz");
   DMass_DR_GroomedJet_CA8_GenWhad_ft[i]->Write();
   DMass_DR_GroomedJet_CA8_GenWhad_ft_c[i]->Print((OutputPlotDir+"/"+DMass_DR_GroomedJet_CA8_GenWhad_ft[i]->GetName()+".png").c_str(),"png");

   ///////////////////////////////////

   DPt_GroomedJet_CA8_WLep_c[i]  = new TCanvas(DPt_GroomedJet_CA8_WLep[i]->GetName());
   DPt_GroomedJet_CA8_WLep_c[i]->cd();
   DPt_GroomedJet_CA8_WLep_c[i]->SetLogy();
   DPt_GroomedJet_CA8_WLep[i]->Draw();
   DPt_GroomedJet_CA8_WLep[i]->Write();
   DPt_GroomedJet_CA8_WLep_c[i]->Print((OutputPlotDir+"/"+DPt_GroomedJet_CA8_WLep[i]->GetName()+".png").c_str(),"png");

   DEta_GroomedJet_CA8_WLep_c[i] = new TCanvas(DEta_GroomedJet_CA8_WLep[i]->GetName());
   DEta_GroomedJet_CA8_WLep_c[i]->cd();
   DEta_GroomedJet_CA8_WLep[i]->Draw();
   DEta_GroomedJet_CA8_WLep[i]->Write();
   DEta_GroomedJet_CA8_WLep_c[i]->Print((OutputPlotDir+"/"+DEta_GroomedJet_CA8_WLep[i]->GetName()+".png").c_str(),"png");

   DPhi_GroomedJet_CA8_WLep_c[i] = new TCanvas(DPhi_GroomedJet_CA8_WLep[i]->GetName());
   DPhi_GroomedJet_CA8_WLep_c[i]->cd();
   DPhi_GroomedJet_CA8_WLep[i]->Draw();
   DPhi_GroomedJet_CA8_WLep[i]->Write();
   DPhi_GroomedJet_CA8_WLep_c[i]->Print((OutputPlotDir+"/"+DPhi_GroomedJet_CA8_WLep[i]->GetName()+".png").c_str(),"png");

   DR_GroomedJet_CA8_WLep_c[i]   = new TCanvas(DR_GroomedJet_CA8_WLep[i]->GetName());
   DR_GroomedJet_CA8_WLep_c[i]->cd();
   DR_GroomedJet_CA8_WLep_c[i]->SetLogy();
   DR_GroomedJet_CA8_WLep[i]->Draw();
   DR_GroomedJet_CA8_WLep[i]->Write();
   DR_GroomedJet_CA8_WLep_c[i]->Print((OutputPlotDir+"/"+DR_GroomedJet_CA8_WLep[i]->GetName()+".png").c_str(),"png");

   DE_GroomedJet_CA8_WLep_c[i]   = new TCanvas(DE_GroomedJet_CA8_WLep[i]->GetName());
   DE_GroomedJet_CA8_WLep_c[i]->cd();
   DE_GroomedJet_CA8_WLep_c[i]->SetLogy();
   DE_GroomedJet_CA8_WLep[i]->Draw();
   DE_GroomedJet_CA8_WLep[i]->Write();
   DE_GroomedJet_CA8_WLep_c[i]->Print((OutputPlotDir+"/"+DE_GroomedJet_CA8_WLep[i]->GetName()+".png").c_str(),"png");

   DMass_GroomedJet_CA8_WLep_c[i]   = new TCanvas(DMass_GroomedJet_CA8_WLep[i]->GetName());
   DMass_GroomedJet_CA8_WLep_c[i]->cd();
   DMass_GroomedJet_CA8_WLep_c[i]->SetLogy();
   DMass_GroomedJet_CA8_WLep[i]->Draw();
   DMass_GroomedJet_CA8_WLep[i]->Write();
   DMass_GroomedJet_CA8_WLep_c[i]->Print((OutputPlotDir+"/"+DMass_GroomedJet_CA8_WLep[i]->GetName()+".png").c_str(),"png");

   DPt_DR_GroomedJet_CA8_WLep_c[i]   = new TCanvas(DPt_DR_GroomedJet_CA8_WLep[i]->GetName());
   DPt_DR_GroomedJet_CA8_WLep_c[i]->cd();
   DPt_DR_GroomedJet_CA8_WLep[i]->Draw("colz");
   DPt_DR_GroomedJet_CA8_WLep[i]->Write();
   DPt_DR_GroomedJet_CA8_WLep_c[i]->Print((OutputPlotDir+"/"+DPt_DR_GroomedJet_CA8_WLep[i]->GetName()+".png").c_str(),"png");

   DMass_DR_GroomedJet_CA8_WLep_c[i]   = new TCanvas(DMass_DR_GroomedJet_CA8_WLep[i]->GetName());
   DMass_DR_GroomedJet_CA8_WLep_c[i]->cd();
   DMass_DR_GroomedJet_CA8_WLep[i]->Draw("colz");
   DMass_DR_GroomedJet_CA8_WLep[i]->Write();
   DMass_DR_GroomedJet_CA8_WLep_c[i]->Print((OutputPlotDir+"/"+DMass_DR_GroomedJet_CA8_WLep[i]->GetName()+".png").c_str(),"png");

  }  

  TCanvas* DPt_WLep_GenWhad_c  = new TCanvas(DPt_WLep_GenWhad->GetName());
  DPt_WLep_GenWhad_c->cd();
  DPt_WLep_GenWhad_c->SetLogy();
  DPt_WLep_GenWhad  ->Draw();
  DPt_WLep_GenWhad  ->Write();
  DPt_WLep_GenWhad_c->Print((OutputPlotDir+"/"+DPt_WLep_GenWhad_c->GetName()+".png").c_str(),"png");

  TCanvas* DEta_WLep_GenWhad_c = new TCanvas(DEta_WLep_GenWhad->GetName());
  DEta_WLep_GenWhad_c->cd();
  DEta_WLep_GenWhad  ->Draw();
  DEta_WLep_GenWhad  ->Write();
  DEta_WLep_GenWhad_c->Print((OutputPlotDir+"/"+DEta_WLep_GenWhad_c->GetName()+".png").c_str(),"png");

  TCanvas* DPhi_WLep_GenWhad_c  = new TCanvas(DPhi_WLep_GenWhad->GetName());
  DPhi_WLep_GenWhad_c->cd();
  DPhi_WLep_GenWhad  ->Draw();
  DPhi_WLep_GenWhad  ->Write();
  DPhi_WLep_GenWhad_c->Print((OutputPlotDir+"/"+DPhi_WLep_GenWhad_c->GetName()+".png").c_str(),"png");

  TCanvas* DR_WLep_GenWhad_c = new TCanvas(DR_WLep_GenWhad->GetName());
  DR_WLep_GenWhad_c->cd();
  DR_WLep_GenWhad_c->SetLogy();
  DR_WLep_GenWhad  ->Draw();
  DR_WLep_GenWhad  ->Write();
  DR_WLep_GenWhad_c->Print((OutputPlotDir+"/"+DR_WLep_GenWhad_c->GetName()+".png").c_str(),"png");

  TCanvas* DE_WLep_GenWhad_c = new TCanvas(DE_WLep_GenWhad->GetName());
  DE_WLep_GenWhad_c->cd();
  DE_WLep_GenWhad_c->SetLogy();
  DE_WLep_GenWhad  ->Draw();
  DE_WLep_GenWhad  ->Write();
  DE_WLep_GenWhad_c->Print((OutputPlotDir+"/"+DE_WLep_GenWhad_c->GetName()+".png").c_str(),"png");

  TCanvas* DMass_WLep_GenWhad_c = new TCanvas(DMass_WLep_GenWhad->GetName());
  DMass_WLep_GenWhad_c->cd();
  DMass_WLep_GenWhad_c->SetLogy();
  DMass_WLep_GenWhad  ->Draw();
  DMass_WLep_GenWhad  ->Write();
  DMass_WLep_GenWhad_c->Print((OutputPlotDir+"/"+DMass_WLep_GenWhad_c->GetName()+".png").c_str(),"png");

  // Canvas for Final Plot after matching selections

  TCanvas**  DPt_GroomedCA8_GenWhad_DRCut_c      = new TCanvas* [NJet];
  TCanvas**  DEta_GroomedCA8_GenWhad_DRCut_c     = new TCanvas* [NJet];
  TCanvas**  DPhi_GroomedCA8_GenWhad_DRCut_c     = new TCanvas* [NJet];
  TCanvas**  DR_GroomedCA8_GenWhad_DRCut_c       = new TCanvas* [NJet];
  TCanvas**  DMass_GroomedCA8_GenWhad_DRCut_c    = new TCanvas* [NJet];

  TCanvas**  DPt_GroomedAK5_GenWhad_DRCut_c      = new TCanvas* [NJet];
  TCanvas**  DEta_GroomedAK5_GenWhad_DRCut_c     = new TCanvas* [NJet];
  TCanvas**  DPhi_GroomedAK5_GenWhad_DRCut_c    = new TCanvas* [NJet];
  TCanvas**  DR_GroomedAK5_GenWhad_DRCut_c       = new TCanvas* [NJet];
  TCanvas**  DMass_GroomedAK5_GenWhad_DRCut_c    = new TCanvas* [NJet];

  TCanvas**  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_c      = new TCanvas* [NJet];
  TCanvas**  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_c     = new TCanvas* [NJet];
  TCanvas**  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_c     = new TCanvas* [NJet];
  TCanvas**  DR_GroomedCA8_GenWhad_DRCut_DPtCut_c       = new TCanvas* [NJet];
  TCanvas**  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_c    = new TCanvas* [NJet];

  TCanvas**  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_c      = new TCanvas* [NJet];
  TCanvas**  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_c     = new TCanvas* [NJet];
  TCanvas**  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_c     = new TCanvas* [NJet];
  TCanvas**  DR_GroomedAK5_GenWhad_DRCut_DPtCut_c       = new TCanvas* [NJet];
  TCanvas**  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_c    = new TCanvas* [NJet];

  for( int i =0; i<NJet ; i++){

   // CA8 After Matching Selections (DR)

   DPt_GroomedCA8_GenWhad_DRCut_c[i]   = new TCanvas(DPt_GroomedCA8_GenWhad_DRCut[i]->GetName());
   DPt_GroomedCA8_GenWhad_DRCut_c[i]->cd();
   DPt_GroomedCA8_GenWhad_DRCut_c[i]->SetLogy();
   DPt_GroomedCA8_GenWhad_DRCut[i]->Draw();
   DPt_GroomedCA8_GenWhad_DRCut[i]->Write();
   DPt_GroomedCA8_GenWhad_DRCut_c[i]->Print((OutputPlotDir+"/"+DPt_GroomedCA8_GenWhad_DRCut[i]->GetName()+".png").c_str(),"png");

   DEta_GroomedCA8_GenWhad_DRCut_c[i]   = new TCanvas(DEta_GroomedCA8_GenWhad_DRCut[i]->GetName());
   DEta_GroomedCA8_GenWhad_DRCut_c[i]->cd();
   DEta_GroomedCA8_GenWhad_DRCut[i]->Draw();
   DEta_GroomedCA8_GenWhad_DRCut[i]->Write();
   DEta_GroomedCA8_GenWhad_DRCut_c[i]->Print((OutputPlotDir+"/"+DEta_GroomedCA8_GenWhad_DRCut[i]->GetName()+".png").c_str(),"png");

   DPhi_GroomedCA8_GenWhad_DRCut_c[i]   = new TCanvas(DPhi_GroomedCA8_GenWhad_DRCut[i]->GetName());
   DPhi_GroomedCA8_GenWhad_DRCut_c[i]->cd();
   DPhi_GroomedCA8_GenWhad_DRCut[i]->Draw();
   DPhi_GroomedCA8_GenWhad_DRCut[i]->Write();
   DPhi_GroomedCA8_GenWhad_DRCut_c[i]->Print((OutputPlotDir+"/"+DPhi_GroomedCA8_GenWhad_DRCut[i]->GetName()+".png").c_str(),"png");

   DR_GroomedCA8_GenWhad_DRCut_c[i]   = new TCanvas(DR_GroomedCA8_GenWhad_DRCut[i]->GetName());
   DR_GroomedCA8_GenWhad_DRCut_c[i]->cd();
   DR_GroomedCA8_GenWhad_DRCut_c[i]->SetLogy();
   DR_GroomedCA8_GenWhad_DRCut[i]->Draw();
   DR_GroomedCA8_GenWhad_DRCut[i]->Write();
   DR_GroomedCA8_GenWhad_DRCut_c[i]->Print((OutputPlotDir+"/"+DR_GroomedCA8_GenWhad_DRCut[i]->GetName()+".png").c_str(),"png");

   DMass_GroomedCA8_GenWhad_DRCut_c[i]   = new TCanvas(DMass_GroomedCA8_GenWhad_DRCut[i]->GetName());
   DMass_GroomedCA8_GenWhad_DRCut_c[i]->cd();
   DMass_GroomedCA8_GenWhad_DRCut_c[i]->SetLogy();
   DMass_GroomedCA8_GenWhad_DRCut[i]->Draw();
   DMass_GroomedCA8_GenWhad_DRCut[i]->Write();
   DMass_GroomedCA8_GenWhad_DRCut_c[i]->Print((OutputPlotDir+"/"+DMass_GroomedCA8_GenWhad_DRCut[i]->GetName()+".png").c_str(),"png");

   // AK5 After Matching Selections  (DR)

   DPt_GroomedAK5_GenWhad_DRCut_c[i]   = new TCanvas(DPt_GroomedAK5_GenWhad_DRCut[i]->GetName());
   DPt_GroomedAK5_GenWhad_DRCut_c[i]->cd();
   DPt_GroomedAK5_GenWhad_DRCut_c[i]->SetLogy();
   DPt_GroomedAK5_GenWhad_DRCut[i]->Draw();
   DPt_GroomedAK5_GenWhad_DRCut[i]->Write();
   DPt_GroomedAK5_GenWhad_DRCut_c[i]->Print((OutputPlotDir+"/"+DPt_GroomedAK5_GenWhad_DRCut[i]->GetName()+".png").c_str(),"png");

   DEta_GroomedAK5_GenWhad_DRCut_c[i]   = new TCanvas(DEta_GroomedAK5_GenWhad_DRCut[i]->GetName());
   DEta_GroomedAK5_GenWhad_DRCut_c[i]->cd();
   DEta_GroomedAK5_GenWhad_DRCut[i]->Draw();
   DEta_GroomedAK5_GenWhad_DRCut[i]->Write();
   DEta_GroomedAK5_GenWhad_DRCut_c[i]->Print((OutputPlotDir+"/"+DEta_GroomedAK5_GenWhad_DRCut[i]->GetName()+".png").c_str(),"png");

   DPhi_GroomedAK5_GenWhad_DRCut_c[i]   = new TCanvas(DPhi_GroomedAK5_GenWhad_DRCut[i]->GetName());
   DPhi_GroomedAK5_GenWhad_DRCut_c[i]->cd();
   DPhi_GroomedAK5_GenWhad_DRCut[i]->Draw();
   DPhi_GroomedAK5_GenWhad_DRCut[i]->Write();
   DPhi_GroomedAK5_GenWhad_DRCut_c[i]->Print((OutputPlotDir+"/"+DPhi_GroomedAK5_GenWhad_DRCut[i]->GetName()+".png").c_str(),"png");

   DR_GroomedAK5_GenWhad_DRCut_c[i]   = new TCanvas(DR_GroomedAK5_GenWhad_DRCut[i]->GetName());
   DR_GroomedAK5_GenWhad_DRCut_c[i]->cd();
   DR_GroomedAK5_GenWhad_DRCut_c[i]->SetLogy();
   DR_GroomedAK5_GenWhad_DRCut[i]->Draw();
   DR_GroomedAK5_GenWhad_DRCut[i]->Write();
   DR_GroomedAK5_GenWhad_DRCut_c[i]->Print((OutputPlotDir+"/"+DR_GroomedAK5_GenWhad_DRCut[i]->GetName()+".png").c_str(),"png");

   DMass_GroomedAK5_GenWhad_DRCut_c[i]   = new TCanvas(DMass_GroomedAK5_GenWhad_DRCut[i]->GetName());
   DMass_GroomedAK5_GenWhad_DRCut_c[i]->cd();
   DMass_GroomedAK5_GenWhad_DRCut_c[i]->SetLogy();
   DMass_GroomedAK5_GenWhad_DRCut[i]->Draw();
   DMass_GroomedAK5_GenWhad_DRCut[i]->Write();
   DMass_GroomedAK5_GenWhad_DRCut_c[i]->Print((OutputPlotDir+"/"+DMass_GroomedAK5_GenWhad_DRCut[i]->GetName()+".png").c_str(),"png");

   // CA8 After Matching Selections (DR,DPt)

   DPt_GroomedCA8_GenWhad_DRCut_DPtCut_c[i]   = new TCanvas(DPt_GroomedCA8_GenWhad_DRCut_DPtCut[i]->GetName());
   DPt_GroomedCA8_GenWhad_DRCut_DPtCut_c[i]->cd();
   DPt_GroomedCA8_GenWhad_DRCut_DPtCut_c[i]->SetLogy();
   DPt_GroomedCA8_GenWhad_DRCut_DPtCut[i]->Draw();
   DPt_GroomedCA8_GenWhad_DRCut_DPtCut[i]->Write();
   DPt_GroomedCA8_GenWhad_DRCut_DPtCut_c[i]->Print((OutputPlotDir+"/"+DPt_GroomedCA8_GenWhad_DRCut_DPtCut[i]->GetName()+".png").c_str(),"png");

   DEta_GroomedCA8_GenWhad_DRCut_DPtCut_c[i]   = new TCanvas(DEta_GroomedCA8_GenWhad_DRCut_DPtCut[i]->GetName());
   DEta_GroomedCA8_GenWhad_DRCut_DPtCut_c[i]->cd();
   DEta_GroomedCA8_GenWhad_DRCut_DPtCut[i]->Draw();
   DEta_GroomedCA8_GenWhad_DRCut_DPtCut[i]->Write();
   DEta_GroomedCA8_GenWhad_DRCut_DPtCut_c[i]->Print((OutputPlotDir+"/"+DEta_GroomedCA8_GenWhad_DRCut_DPtCut[i]->GetName()+".png").c_str(),"png");

   DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_c[i]   = new TCanvas(DPhi_GroomedCA8_GenWhad_DRCut_DPtCut[i]->GetName());
   DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_c[i]->cd();
   DPhi_GroomedCA8_GenWhad_DRCut_DPtCut[i]->Draw();
   DPhi_GroomedCA8_GenWhad_DRCut_DPtCut[i]->Write();
   DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_c[i]->Print((OutputPlotDir+"/"+DPhi_GroomedCA8_GenWhad_DRCut_DPtCut[i]->GetName()+".png").c_str(),"png");

   DR_GroomedCA8_GenWhad_DRCut_DPtCut_c[i]   = new TCanvas(DR_GroomedCA8_GenWhad_DRCut_DPtCut[i]->GetName());
   DR_GroomedCA8_GenWhad_DRCut_DPtCut_c[i]->cd();
   DR_GroomedCA8_GenWhad_DRCut_DPtCut_c[i]->SetLogy();
   DR_GroomedCA8_GenWhad_DRCut_DPtCut[i]->Draw();
   DR_GroomedCA8_GenWhad_DRCut_DPtCut[i]->Write();
   DR_GroomedCA8_GenWhad_DRCut_DPtCut_c[i]->Print((OutputPlotDir+"/"+DR_GroomedCA8_GenWhad_DRCut_DPtCut[i]->GetName()+".png").c_str(),"png");

   DMass_GroomedCA8_GenWhad_DRCut_DPtCut_c[i]   = new TCanvas(DMass_GroomedCA8_GenWhad_DRCut_DPtCut[i]->GetName());
   DMass_GroomedCA8_GenWhad_DRCut_DPtCut_c[i]->cd();
   DMass_GroomedCA8_GenWhad_DRCut_DPtCut_c[i]->SetLogy();
   DMass_GroomedCA8_GenWhad_DRCut_DPtCut[i]->Draw();
   DMass_GroomedCA8_GenWhad_DRCut_DPtCut[i]->Write();
   DMass_GroomedCA8_GenWhad_DRCut_DPtCut_c[i]->Print((OutputPlotDir+"/"+DMass_GroomedCA8_GenWhad_DRCut_DPtCut[i]->GetName()+".png").c_str(),"png");

   // AK5 After Matching Selections (DR,DPt)

   DPt_GroomedAK5_GenWhad_DRCut_DPtCut_c[i]   = new TCanvas(DPt_GroomedAK5_GenWhad_DRCut_DPtCut[i]->GetName());
   DPt_GroomedAK5_GenWhad_DRCut_DPtCut_c[i]->cd();
   DPt_GroomedAK5_GenWhad_DRCut_DPtCut_c[i]->SetLogy();
   DPt_GroomedAK5_GenWhad_DRCut_DPtCut[i]->Draw();
   DPt_GroomedAK5_GenWhad_DRCut_DPtCut[i]->Write();
   DPt_GroomedAK5_GenWhad_DRCut_DPtCut_c[i]->Print((OutputPlotDir+"/"+DPt_GroomedAK5_GenWhad_DRCut_DPtCut[i]->GetName()+".png").c_str(),"png");

   DEta_GroomedAK5_GenWhad_DRCut_DPtCut_c[i]   = new TCanvas(DEta_GroomedAK5_GenWhad_DRCut_DPtCut[i]->GetName());
   DEta_GroomedAK5_GenWhad_DRCut_DPtCut_c[i]->cd();
   DEta_GroomedAK5_GenWhad_DRCut_DPtCut[i]->Draw();
   DEta_GroomedAK5_GenWhad_DRCut_DPtCut[i]->Write();
   DEta_GroomedAK5_GenWhad_DRCut_DPtCut_c[i]->Print((OutputPlotDir+"/"+DEta_GroomedAK5_GenWhad_DRCut_DPtCut[i]->GetName()+".png").c_str(),"png");

   DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_c[i]   = new TCanvas(DPhi_GroomedAK5_GenWhad_DRCut_DPtCut[i]->GetName());
   DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_c[i]->cd();
   DPhi_GroomedAK5_GenWhad_DRCut_DPtCut[i]->Draw();
   DPhi_GroomedAK5_GenWhad_DRCut_DPtCut[i]->Write();
   DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_c[i]->Print((OutputPlotDir+"/"+DPhi_GroomedAK5_GenWhad_DRCut_DPtCut[i]->GetName()+".png").c_str(),"png");

   DR_GroomedAK5_GenWhad_DRCut_DPtCut_c[i]   = new TCanvas(DR_GroomedAK5_GenWhad_DRCut_DPtCut[i]->GetName());
   DR_GroomedAK5_GenWhad_DRCut_DPtCut_c[i]->cd();
   DR_GroomedAK5_GenWhad_DRCut_DPtCut_c[i]->SetLogy();
   DR_GroomedAK5_GenWhad_DRCut_DPtCut[i]->Draw();
   DR_GroomedAK5_GenWhad_DRCut_DPtCut[i]->Write();
   DR_GroomedAK5_GenWhad_DRCut_DPtCut_c[i]->Print((OutputPlotDir+"/"+DR_GroomedAK5_GenWhad_DRCut_DPtCut[i]->GetName()+".png").c_str(),"png");

   DMass_GroomedAK5_GenWhad_DRCut_DPtCut_c[i]   = new TCanvas(DMass_GroomedAK5_GenWhad_DRCut_DPtCut[i]->GetName());
   DMass_GroomedAK5_GenWhad_DRCut_DPtCut_c[i]->cd();
   DMass_GroomedAK5_GenWhad_DRCut_DPtCut_c[i]->SetLogy();
   DMass_GroomedAK5_GenWhad_DRCut_DPtCut[i]->Draw();
   DMass_GroomedAK5_GenWhad_DRCut_DPtCut[i]->Write();
   DMass_GroomedAK5_GenWhad_DRCut_DPtCut_c[i]->Print((OutputPlotDir+"/"+DMass_GroomedAK5_GenWhad_DRCut_DPtCut[i]->GetName()+".png").c_str(),"png");


  }


  // CA8 After exclusive selection of events with only one CA8 Jet (DR)

  TCanvas*  DPt_GroomedCA8_GenWhad_DRCut_1Jet_c    = new TCanvas (DPt_GroomedCA8_GenWhad_DRCut_1Jet->GetName());
  DPt_GroomedCA8_GenWhad_DRCut_1Jet_c->cd();
  DPt_GroomedCA8_GenWhad_DRCut_1Jet_c->SetLogy();
  DPt_GroomedCA8_GenWhad_DRCut_1Jet->Draw();
  DPt_GroomedCA8_GenWhad_DRCut_1Jet->Write();
  DPt_GroomedCA8_GenWhad_DRCut_1Jet_c->Print((OutputPlotDir+"/"+DPt_GroomedCA8_GenWhad_DRCut_1Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedCA8_GenWhad_DRCut_1Jet_c    = new TCanvas (DEta_GroomedCA8_GenWhad_DRCut_1Jet->GetName());
  DEta_GroomedCA8_GenWhad_DRCut_1Jet_c->cd();
  DEta_GroomedCA8_GenWhad_DRCut_1Jet_c->SetLogy();
  DEta_GroomedCA8_GenWhad_DRCut_1Jet->Draw();
  DEta_GroomedCA8_GenWhad_DRCut_1Jet->Write();
  DEta_GroomedCA8_GenWhad_DRCut_1Jet_c->Print((OutputPlotDir+"/"+DEta_GroomedCA8_GenWhad_DRCut_1Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedCA8_GenWhad_DRCut_1Jet_c    = new TCanvas (DPhi_GroomedCA8_GenWhad_DRCut_1Jet->GetName());
  DPhi_GroomedCA8_GenWhad_DRCut_1Jet_c->cd();
  DPhi_GroomedCA8_GenWhad_DRCut_1Jet_c->SetLogy();
  DPhi_GroomedCA8_GenWhad_DRCut_1Jet->Draw();
  DPhi_GroomedCA8_GenWhad_DRCut_1Jet->Write();
  DPhi_GroomedCA8_GenWhad_DRCut_1Jet_c->Print((OutputPlotDir+"/"+DPhi_GroomedCA8_GenWhad_DRCut_1Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedCA8_GenWhad_DRCut_1Jet_c    = new TCanvas (DR_GroomedCA8_GenWhad_DRCut_1Jet->GetName());
  DR_GroomedCA8_GenWhad_DRCut_1Jet_c->cd();
  DR_GroomedCA8_GenWhad_DRCut_1Jet_c->SetLogy();
  DR_GroomedCA8_GenWhad_DRCut_1Jet->Draw();
  DR_GroomedCA8_GenWhad_DRCut_1Jet->Write();
  DR_GroomedCA8_GenWhad_DRCut_1Jet_c->Print((OutputPlotDir+"/"+DR_GroomedCA8_GenWhad_DRCut_1Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedCA8_GenWhad_DRCut_1Jet_c    = new TCanvas (DMass_GroomedCA8_GenWhad_DRCut_1Jet->GetName());
  DMass_GroomedCA8_GenWhad_DRCut_1Jet_c->cd();
  DMass_GroomedCA8_GenWhad_DRCut_1Jet_c->SetLogy();
  DMass_GroomedCA8_GenWhad_DRCut_1Jet->Draw();
  DMass_GroomedCA8_GenWhad_DRCut_1Jet->Write();
  DMass_GroomedCA8_GenWhad_DRCut_1Jet_c->Print((OutputPlotDir+"/"+DMass_GroomedCA8_GenWhad_DRCut_1Jet->GetName()+".png").c_str(),"png");


  // AK5 After exclusive selection of events with only one AK5 Jet (DR)

  TCanvas*  DPt_GroomedAK5_GenWhad_DRCut_1Jet_c    = new TCanvas (DPt_GroomedAK5_GenWhad_DRCut_1Jet->GetName());
  DPt_GroomedAK5_GenWhad_DRCut_1Jet_c->cd();
  DPt_GroomedAK5_GenWhad_DRCut_1Jet_c->SetLogy();
  DPt_GroomedAK5_GenWhad_DRCut_1Jet->Draw();
  DPt_GroomedAK5_GenWhad_DRCut_1Jet->Write();
  DPt_GroomedAK5_GenWhad_DRCut_1Jet_c->Print((OutputPlotDir+"/"+DPt_GroomedAK5_GenWhad_DRCut_1Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedAK5_GenWhad_DRCut_1Jet_c    = new TCanvas (DEta_GroomedAK5_GenWhad_DRCut_1Jet->GetName());
  DEta_GroomedAK5_GenWhad_DRCut_1Jet_c->cd();
  DEta_GroomedAK5_GenWhad_DRCut_1Jet_c->SetLogy();
  DEta_GroomedAK5_GenWhad_DRCut_1Jet->Draw();
  DEta_GroomedAK5_GenWhad_DRCut_1Jet->Write();
  DEta_GroomedAK5_GenWhad_DRCut_1Jet_c->Print((OutputPlotDir+"/"+DEta_GroomedAK5_GenWhad_DRCut_1Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedAK5_GenWhad_DRCut_1Jet_c    = new TCanvas (DPhi_GroomedAK5_GenWhad_DRCut_1Jet->GetName());
  DPhi_GroomedAK5_GenWhad_DRCut_1Jet_c->cd();
  DPhi_GroomedAK5_GenWhad_DRCut_1Jet_c->SetLogy();
  DPhi_GroomedAK5_GenWhad_DRCut_1Jet->Draw();
  DPhi_GroomedAK5_GenWhad_DRCut_1Jet->Write();
  DPhi_GroomedAK5_GenWhad_DRCut_1Jet_c->Print((OutputPlotDir+"/"+DPhi_GroomedAK5_GenWhad_DRCut_1Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedAK5_GenWhad_DRCut_1Jet_c    = new TCanvas (DR_GroomedAK5_GenWhad_DRCut_1Jet->GetName());
  DR_GroomedAK5_GenWhad_DRCut_1Jet_c->cd();
  DR_GroomedAK5_GenWhad_DRCut_1Jet_c->SetLogy();
  DR_GroomedAK5_GenWhad_DRCut_1Jet->Draw();
  DR_GroomedAK5_GenWhad_DRCut_1Jet->Write();
  DR_GroomedAK5_GenWhad_DRCut_1Jet_c->Print((OutputPlotDir+"/"+DR_GroomedAK5_GenWhad_DRCut_1Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedAK5_GenWhad_DRCut_1Jet_c    = new TCanvas (DMass_GroomedAK5_GenWhad_DRCut_1Jet->GetName());
  DMass_GroomedAK5_GenWhad_DRCut_1Jet_c->cd();
  DMass_GroomedAK5_GenWhad_DRCut_1Jet_c->SetLogy();
  DMass_GroomedAK5_GenWhad_DRCut_1Jet->Draw();
  DMass_GroomedAK5_GenWhad_DRCut_1Jet->Write();
  DMass_GroomedAK5_GenWhad_DRCut_1Jet_c->Print((OutputPlotDir+"/"+DMass_GroomedAK5_GenWhad_DRCut_1Jet->GetName()+".png").c_str(),"png");


  // CA8 After exclusive selection of events with only one CA8 Jet (DR)

  TCanvas*  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet_c    = new TCanvas (DPt_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet->GetName());
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet_c->cd();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet_c->SetLogy();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet->Draw();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet->Write();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet_c->Print((OutputPlotDir+"/"+DPt_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet_c    = new TCanvas (DEta_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet->GetName());
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet_c->cd();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet_c->SetLogy();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet->Draw();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet->Write();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet_c->Print((OutputPlotDir+"/"+DEta_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet_c    = new TCanvas (DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet->GetName());
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet_c->cd();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet_c->SetLogy();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet->Draw();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet->Write();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet_c->Print((OutputPlotDir+"/"+DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet_c    = new TCanvas (DR_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet->GetName());
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet_c->cd();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet_c->SetLogy();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet->Draw();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet->Write();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet_c->Print((OutputPlotDir+"/"+DR_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet_c    = new TCanvas (DMass_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet->GetName());
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet_c->cd();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet_c->SetLogy();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet->Draw();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet->Write();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet_c->Print((OutputPlotDir+"/"+DMass_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet->GetName()+".png").c_str(),"png");


  // AK5 After exclusive selection of events with only one AK5 Jet (DR)

  TCanvas*  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet_c    = new TCanvas (DPt_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->GetName());
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet_c->cd();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet_c->SetLogy();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->Draw();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->Write();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet_c->Print((OutputPlotDir+"/"+DPt_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet_c    = new TCanvas (DEta_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->GetName());
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet_c->cd();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet_c->SetLogy();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->Draw();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->Write();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet_c->Print((OutputPlotDir+"/"+DEta_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet_c    = new TCanvas (DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->GetName());
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet_c->cd();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet_c->SetLogy();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->Draw();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->Write();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet_c->Print((OutputPlotDir+"/"+DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet_c    = new TCanvas (DR_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->GetName());
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet_c->cd();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet_c->SetLogy();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->Draw();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->Write();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet_c->Print((OutputPlotDir+"/"+DR_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet_c    = new TCanvas (DMass_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->GetName());
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet_c->cd();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet_c->SetLogy();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->Draw();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->Write();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet_c->Print((OutputPlotDir+"/"+DMass_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->GetName()+".png").c_str(),"png");


  // CA8 After exclusive selection of events with only one CA8 Jet (DR)

  TCanvas*  DPt_GroomedCA8_GenWhad_DRCut_2Jet_c    = new TCanvas (DPt_GroomedCA8_GenWhad_DRCut_2Jet->GetName());
  DPt_GroomedCA8_GenWhad_DRCut_2Jet_c->cd();
  DPt_GroomedCA8_GenWhad_DRCut_2Jet_c->SetLogy();
  DPt_GroomedCA8_GenWhad_DRCut_2Jet->Draw();
  DPt_GroomedCA8_GenWhad_DRCut_2Jet->Write();
  DPt_GroomedCA8_GenWhad_DRCut_2Jet_c->Print((OutputPlotDir+"/"+DPt_GroomedCA8_GenWhad_DRCut_2Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedCA8_GenWhad_DRCut_2Jet_c    = new TCanvas (DEta_GroomedCA8_GenWhad_DRCut_2Jet->GetName());
  DEta_GroomedCA8_GenWhad_DRCut_2Jet_c->cd();
  DEta_GroomedCA8_GenWhad_DRCut_2Jet_c->SetLogy();
  DEta_GroomedCA8_GenWhad_DRCut_2Jet->Draw();
  DEta_GroomedCA8_GenWhad_DRCut_2Jet->Write();
  DEta_GroomedCA8_GenWhad_DRCut_2Jet_c->Print((OutputPlotDir+"/"+DEta_GroomedCA8_GenWhad_DRCut_2Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedCA8_GenWhad_DRCut_2Jet_c    = new TCanvas (DPhi_GroomedCA8_GenWhad_DRCut_2Jet->GetName());
  DPhi_GroomedCA8_GenWhad_DRCut_2Jet_c->cd();
  DPhi_GroomedCA8_GenWhad_DRCut_2Jet_c->SetLogy();
  DPhi_GroomedCA8_GenWhad_DRCut_2Jet->Draw();
  DPhi_GroomedCA8_GenWhad_DRCut_2Jet->Write();
  DPhi_GroomedCA8_GenWhad_DRCut_2Jet_c->Print((OutputPlotDir+"/"+DPhi_GroomedCA8_GenWhad_DRCut_2Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedCA8_GenWhad_DRCut_2Jet_c    = new TCanvas (DR_GroomedCA8_GenWhad_DRCut_2Jet->GetName());
  DR_GroomedCA8_GenWhad_DRCut_2Jet_c->cd();
  DR_GroomedCA8_GenWhad_DRCut_2Jet_c->SetLogy();
  DR_GroomedCA8_GenWhad_DRCut_2Jet->Draw();
  DR_GroomedCA8_GenWhad_DRCut_2Jet->Write();
  DR_GroomedCA8_GenWhad_DRCut_2Jet_c->Print((OutputPlotDir+"/"+DR_GroomedCA8_GenWhad_DRCut_2Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedCA8_GenWhad_DRCut_2Jet_c    = new TCanvas (DMass_GroomedCA8_GenWhad_DRCut_2Jet->GetName());
  DMass_GroomedCA8_GenWhad_DRCut_2Jet_c->cd();
  DMass_GroomedCA8_GenWhad_DRCut_2Jet_c->SetLogy();
  DMass_GroomedCA8_GenWhad_DRCut_2Jet->Draw();
  DMass_GroomedCA8_GenWhad_DRCut_2Jet->Write();
  DMass_GroomedCA8_GenWhad_DRCut_2Jet_c->Print((OutputPlotDir+"/"+DMass_GroomedCA8_GenWhad_DRCut_2Jet->GetName()+".png").c_str(),"png");


  // AK5 After exclusive selection of events with only one AK5 Jet (DR)

  TCanvas*  DPt_GroomedAK5_GenWhad_DRCut_2Jet_c    = new TCanvas (DPt_GroomedAK5_GenWhad_DRCut_2Jet->GetName());
  DPt_GroomedAK5_GenWhad_DRCut_2Jet_c->cd();
  DPt_GroomedAK5_GenWhad_DRCut_2Jet_c->SetLogy();
  DPt_GroomedAK5_GenWhad_DRCut_2Jet->Draw();
  DPt_GroomedAK5_GenWhad_DRCut_2Jet->Write();
  DPt_GroomedAK5_GenWhad_DRCut_2Jet_c->Print((OutputPlotDir+"/"+DPt_GroomedAK5_GenWhad_DRCut_2Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedAK5_GenWhad_DRCut_2Jet_c    = new TCanvas (DEta_GroomedAK5_GenWhad_DRCut_2Jet->GetName());
  DEta_GroomedAK5_GenWhad_DRCut_2Jet_c->cd();
  DEta_GroomedAK5_GenWhad_DRCut_2Jet_c->SetLogy();
  DEta_GroomedAK5_GenWhad_DRCut_2Jet->Draw();
  DEta_GroomedAK5_GenWhad_DRCut_2Jet->Write();
  DEta_GroomedAK5_GenWhad_DRCut_2Jet_c->Print((OutputPlotDir+"/"+DEta_GroomedAK5_GenWhad_DRCut_2Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedAK5_GenWhad_DRCut_2Jet_c    = new TCanvas (DPhi_GroomedAK5_GenWhad_DRCut_2Jet->GetName());
  DPhi_GroomedAK5_GenWhad_DRCut_2Jet_c->cd();
  DPhi_GroomedAK5_GenWhad_DRCut_2Jet_c->SetLogy();
  DPhi_GroomedAK5_GenWhad_DRCut_2Jet->Draw();
  DPhi_GroomedAK5_GenWhad_DRCut_2Jet->Write();
  DPhi_GroomedAK5_GenWhad_DRCut_2Jet_c->Print((OutputPlotDir+"/"+DPhi_GroomedAK5_GenWhad_DRCut_2Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedAK5_GenWhad_DRCut_2Jet_c    = new TCanvas (DR_GroomedAK5_GenWhad_DRCut_2Jet->GetName());
  DR_GroomedAK5_GenWhad_DRCut_2Jet_c->cd();
  DR_GroomedAK5_GenWhad_DRCut_2Jet_c->SetLogy();
  DR_GroomedAK5_GenWhad_DRCut_2Jet->Draw();
  DR_GroomedAK5_GenWhad_DRCut_2Jet->Write();
  DR_GroomedAK5_GenWhad_DRCut_2Jet_c->Print((OutputPlotDir+"/"+DR_GroomedAK5_GenWhad_DRCut_2Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedAK5_GenWhad_DRCut_2Jet_c    = new TCanvas (DMass_GroomedAK5_GenWhad_DRCut_2Jet->GetName());
  DMass_GroomedAK5_GenWhad_DRCut_2Jet_c->cd();
  DMass_GroomedAK5_GenWhad_DRCut_2Jet_c->SetLogy();
  DMass_GroomedAK5_GenWhad_DRCut_2Jet->Draw();
  DMass_GroomedAK5_GenWhad_DRCut_2Jet->Write();
  DMass_GroomedAK5_GenWhad_DRCut_2Jet_c->Print((OutputPlotDir+"/"+DMass_GroomedAK5_GenWhad_DRCut_2Jet->GetName()+".png").c_str(),"png");



  // CA8 After exclusive selection of events with only two CA8 Jet (DR)

  TCanvas*  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_c    = new TCanvas (DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->GetName());
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_c->cd();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_c->SetLogy();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->Draw();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->Write();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_c->Print((OutputPlotDir+"/"+DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_c    = new TCanvas (DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->GetName());
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_c->cd();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_c->SetLogy();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->Draw();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->Write();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_c->Print((OutputPlotDir+"/"+DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_c    = new TCanvas (DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->GetName());
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_c->cd();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_c->SetLogy();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->Draw();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->Write();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_c->Print((OutputPlotDir+"/"+DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_c    = new TCanvas (DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->GetName());
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_c->cd();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_c->SetLogy();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->Draw();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->Write();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_c->Print((OutputPlotDir+"/"+DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_c    = new TCanvas (DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->GetName());
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_c->cd();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_c->SetLogy();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->Draw();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->Write();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_c->Print((OutputPlotDir+"/"+DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->GetName()+".png").c_str(),"png");

  // CA8 After exclusive selection of events with only two CA8 Jet (DR)

  TCanvas*  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_c    = new TCanvas (DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->GetName());
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_c->cd();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_c->SetLogy();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->Draw();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->Write();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_c->Print((OutputPlotDir+"/"+DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_c    = new TCanvas (DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->GetName());
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_c->cd();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_c->SetLogy();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->Draw();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->Write();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_c->Print((OutputPlotDir+"/"+DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_c    = new TCanvas (DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->GetName());
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_c->cd();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_c->SetLogy();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->Draw();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->Write();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_c->Print((OutputPlotDir+"/"+DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_c    = new TCanvas (DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->GetName());
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_c->cd();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_c->SetLogy();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->Draw();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->Write();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_c->Print((OutputPlotDir+"/"+DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_c    = new TCanvas (DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->GetName());
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_c->cd();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_c->SetLogy();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->Draw();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->Write();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_c->Print((OutputPlotDir+"/"+DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->GetName()+".png").c_str(),"png");

  // CA8 After exclusive selection of events with only two CA8 Jet (DR) --> both in the DR --> Double Matching

  TCanvas*  DPt_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching_c    = new TCanvas (DPt_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->GetName());
  DPt_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching_c->cd();
  DPt_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching_c->SetLogy();
  DPt_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->Draw();
  DPt_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->Write();
  DPt_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DPt_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching_c    = new TCanvas (DEta_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->GetName());
  DEta_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching_c->cd();
  DEta_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching_c->SetLogy();
  DEta_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->Draw();
  DEta_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->Write();
  DEta_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DEta_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching_c    = new TCanvas (DPhi_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->GetName());
  DPhi_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching_c->cd();
  DPhi_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching_c->SetLogy();
  DPhi_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->Draw();
  DPhi_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->Write();
  DPhi_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DPhi_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching_c    = new TCanvas (DR_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->GetName());
  DR_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching_c->cd();
  DR_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching_c->SetLogy();
  DR_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->Draw();
  DR_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->Write();
  DR_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DR_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching_c    = new TCanvas (DMass_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->GetName());
  DMass_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching_c->cd();
  DMass_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching_c->SetLogy();
  DMass_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->Draw();
  DMass_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->Write();
  DMass_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DMass_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  // AK5 After exclusive selection of events with only two CA8 Jet (DR) --> both in the DR --> Double Matching

  TCanvas*  DPt_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching_c    = new TCanvas (DPt_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->GetName());
  DPt_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching_c->cd();
  DPt_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching_c->SetLogy();
  DPt_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->Draw();
  DPt_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->Write();
  DPt_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DPt_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching_c    = new TCanvas (DEta_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->GetName());
  DEta_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching_c->cd();
  DEta_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching_c->SetLogy();
  DEta_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->Draw();
  DEta_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->Write();
  DEta_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DEta_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching_c    = new TCanvas (DPhi_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->GetName());
  DPhi_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching_c->cd();
  DPhi_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching_c->SetLogy();
  DPhi_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->Draw();
  DPhi_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->Write();
  DPhi_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DPhi_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching_c    = new TCanvas (DR_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->GetName());
  DR_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching_c->cd();
  DR_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching_c->SetLogy();
  DR_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->Draw();
  DR_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->Write();
  DR_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DR_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching_c    = new TCanvas (DMass_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->GetName());
  DMass_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching_c->cd();
  DMass_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching_c->SetLogy();
  DMass_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->Draw();
  DMass_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->Write();
  DMass_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DMass_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c    = new TCanvas (DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetName());
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->cd();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->SetLogy();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->Draw();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->Write();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c    = new TCanvas (DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetName());
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->cd();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->SetLogy();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->Draw();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->Write();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c    = new TCanvas (DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetName());
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->cd();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->SetLogy();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->Draw();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->Write();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c    = new TCanvas (DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetName());
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->cd();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->SetLogy();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->Draw();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->Write();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c    = new TCanvas (DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetName());
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->cd();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->SetLogy();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->Draw();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->Write();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetName()+".png").c_str(),"png");


  TCanvas*  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c    = new TCanvas (DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetName());
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->cd();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->SetLogy();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->Draw();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->Write();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c    = new TCanvas (DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetName());
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->cd();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->SetLogy();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->Draw();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->Write();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c    = new TCanvas (DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetName());
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->cd();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->SetLogy();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->Draw();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->Write();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c    = new TCanvas (DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetName());
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->cd();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->SetLogy();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->Draw();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->Write();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c    = new TCanvas (DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetName());
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->cd();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->SetLogy();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->Draw();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->Write();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetName()+".png").c_str(),"png");


  // CA8 After exclusive selection of events with only two CA8 Jet (DR) --> both in the DR --> Double Matching

  TCanvas*  DPt_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching_c    = new TCanvas (DPt_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->GetName());
  DPt_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching_c->cd();
  DPt_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching_c->SetLogy();
  DPt_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->Draw();
  DPt_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->Write();
  DPt_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DPt_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching_c    = new TCanvas (DEta_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->GetName());
  DEta_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching_c->cd();
  DEta_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching_c->SetLogy();
  DEta_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->Draw();
  DEta_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->Write();
  DEta_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DEta_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching_c    = new TCanvas (DPhi_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->GetName());
  DPhi_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching_c->cd();
  DPhi_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching_c->SetLogy();
  DPhi_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->Draw();
  DPhi_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->Write();
  DPhi_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DPhi_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching_c    = new TCanvas (DR_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->GetName());
  DR_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching_c->cd();
  DR_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching_c->SetLogy();
  DR_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->Draw();
  DR_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->Write();
  DR_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DR_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching_c    = new TCanvas (DMass_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->GetName());
  DMass_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching_c->cd();
  DMass_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching_c->SetLogy();
  DMass_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->Draw();
  DMass_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->Write();
  DMass_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DMass_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPt_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching_c    = new TCanvas (DPt_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->GetName());
  DPt_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching_c->cd();
  DPt_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching_c->SetLogy();
  DPt_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->Draw();
  DPt_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->Write();
  DPt_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DPt_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching_c    = new TCanvas (DEta_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->GetName());
  DEta_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching_c->cd();
  DEta_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching_c->SetLogy();
  DEta_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->Draw();
  DEta_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->Write();
  DEta_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DEta_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching_c    = new TCanvas (DPhi_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->GetName());
  DPhi_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching_c->cd();
  DPhi_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching_c->SetLogy();
  DPhi_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->Draw();
  DPhi_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->Write();
  DPhi_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DPhi_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching_c    = new TCanvas (DR_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->GetName());
  DR_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching_c->cd();
  DR_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching_c->SetLogy();
  DR_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->Draw();
  DR_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->Write();
  DR_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DR_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching_c    = new TCanvas (DMass_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->GetName());
  DMass_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching_c->cd();
  DMass_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching_c->SetLogy();
  DMass_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->Draw();
  DMass_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->Write();
  DMass_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DMass_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c    = new TCanvas (DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetName());
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->cd();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->SetLogy();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->Draw();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->Write();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c    = new TCanvas (DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetName());
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->cd();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->SetLogy();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->Draw();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->Write();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c    = new TCanvas (DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetName());
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->cd();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->SetLogy();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->Draw();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->Write();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c    = new TCanvas (DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetName());
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->cd();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->SetLogy();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->Draw();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->Write();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c    = new TCanvas (DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetName());
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->cd();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->SetLogy();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->Draw();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->Write();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetName()+".png").c_str(),"png");


  TCanvas*  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c    = new TCanvas (DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetName());
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->cd();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->SetLogy();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->Draw();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->Write();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c    = new TCanvas (DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetName());
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->cd();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->SetLogy();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->Draw();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->Write();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c    = new TCanvas (DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetName());
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->cd();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->SetLogy();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->Draw();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->Write();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c    = new TCanvas (DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetName());
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->cd();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->SetLogy();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->Draw();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->Write();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c    = new TCanvas (DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetName());
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->cd();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->SetLogy();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->Draw();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->Write();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  // CA8 After exclusive selection of events with only two CA8 Jet (DR) --> both out the DR --> No Matching

  TCanvas*  DPt_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching_c    = new TCanvas (DPt_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->GetName());
  DPt_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching_c->cd();
  DPt_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching_c->SetLogy();
  DPt_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->Draw();
  DPt_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->Write();
  DPt_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching_c->Print((OutputPlotDir+"/"+DPt_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching_c    = new TCanvas (DEta_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->GetName());
  DEta_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching_c->cd();
  DEta_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching_c->SetLogy();
  DEta_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->Draw();
  DEta_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->Write();
  DEta_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching_c->Print((OutputPlotDir+"/"+DEta_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching_c    = new TCanvas (DPhi_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->GetName());
  DPhi_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching_c->cd();
  DPhi_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching_c->SetLogy();
  DPhi_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->Draw();
  DPhi_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->Write();
  DPhi_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching_c->Print((OutputPlotDir+"/"+DPhi_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching_c    = new TCanvas (DR_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->GetName());
  DR_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching_c->cd();
  DR_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching_c->SetLogy();
  DR_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->Draw();
  DR_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->Write();
  DR_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching_c->Print((OutputPlotDir+"/"+DR_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching_c    = new TCanvas (DMass_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->GetName());
  DMass_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching_c->cd();
  DMass_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching_c->SetLogy();
  DMass_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->Draw();
  DMass_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->Write();
  DMass_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching_c->Print((OutputPlotDir+"/"+DMass_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPt_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching_c    = new TCanvas (DPt_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->GetName());
  DPt_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching_c->cd();
  DPt_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching_c->SetLogy();
  DPt_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->Draw();
  DPt_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->Write();
  DPt_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching_c->Print((OutputPlotDir+"/"+DPt_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching_c    = new TCanvas (DEta_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->GetName());
  DEta_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching_c->cd();
  DEta_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching_c->SetLogy();
  DEta_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->Draw();
  DEta_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->Write();
  DEta_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching_c->Print((OutputPlotDir+"/"+DEta_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching_c    = new TCanvas (DPhi_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->GetName());
  DPhi_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching_c->cd();
  DPhi_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching_c->SetLogy();
  DPhi_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->Draw();
  DPhi_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->Write();
  DPhi_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching_c->Print((OutputPlotDir+"/"+DPhi_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching_c    = new TCanvas (DR_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->GetName());
  DR_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching_c->cd();
  DR_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching_c->SetLogy();
  DR_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->Draw();
  DR_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->Write();
  DR_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching_c->Print((OutputPlotDir+"/"+DR_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching_c    = new TCanvas (DMass_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->GetName());
  DMass_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching_c->cd();
  DMass_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching_c->SetLogy();
  DMass_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->Draw();
  DMass_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->Write();
  DMass_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching_c->Print((OutputPlotDir+"/"+DMass_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c    = new TCanvas (DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetName());
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->cd();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->SetLogy();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching->Draw();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching->Write();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->Print((OutputPlotDir+"/"+DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c    = new TCanvas (DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetName());
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->cd();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->SetLogy();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching->Draw();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching->Write();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->Print((OutputPlotDir+"/"+DEta_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c    = new TCanvas (DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetName());
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->cd();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->SetLogy();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching->Draw();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching->Write();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->Print((OutputPlotDir+"/"+DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c    = new TCanvas (DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetName());
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->cd();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->SetLogy();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching->Draw();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching->Write();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->Print((OutputPlotDir+"/"+DR_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c    = new TCanvas (DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetName());
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->cd();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->SetLogy();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching->Draw();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching->Write();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->Print((OutputPlotDir+"/"+DMass_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetName()+".png").c_str(),"png");


  TCanvas*  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c    = new TCanvas (DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetName());
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->cd();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->SetLogy();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching->Draw();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching->Write();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->Print((OutputPlotDir+"/"+DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c    = new TCanvas (DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetName());
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->cd();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->SetLogy();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching->Draw();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching->Write();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->Print((OutputPlotDir+"/"+DEta_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c    = new TCanvas (DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetName());
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->cd();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->SetLogy();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching->Draw();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching->Write();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->Print((OutputPlotDir+"/"+DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c    = new TCanvas (DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetName());
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->cd();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->SetLogy();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching->Draw();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching->Write();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->Print((OutputPlotDir+"/"+DR_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c    = new TCanvas (DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetName());
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->cd();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->SetLogy();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching->Draw();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching->Write();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching_c->Print((OutputPlotDir+"/"+DMass_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetName()+".png").c_str(),"png");




  // CA8 After exclusive selection of events with only three CA8 Jet (DR)

  TCanvas*  DPt_GroomedCA8_GenWhad_DRCut_3Jet_c    = new TCanvas (DPt_GroomedCA8_GenWhad_DRCut_3Jet->GetName());
  DPt_GroomedCA8_GenWhad_DRCut_3Jet_c->cd();
  DPt_GroomedCA8_GenWhad_DRCut_3Jet_c->SetLogy();
  DPt_GroomedCA8_GenWhad_DRCut_3Jet->Draw();
  DPt_GroomedCA8_GenWhad_DRCut_3Jet->Write();
  DPt_GroomedCA8_GenWhad_DRCut_3Jet_c->Print((OutputPlotDir+"/"+DPt_GroomedCA8_GenWhad_DRCut_3Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedCA8_GenWhad_DRCut_3Jet_c    = new TCanvas (DEta_GroomedCA8_GenWhad_DRCut_3Jet->GetName());
  DEta_GroomedCA8_GenWhad_DRCut_3Jet_c->cd();
  DEta_GroomedCA8_GenWhad_DRCut_3Jet_c->SetLogy();
  DEta_GroomedCA8_GenWhad_DRCut_3Jet->Draw();
  DEta_GroomedCA8_GenWhad_DRCut_3Jet->Write();
  DEta_GroomedCA8_GenWhad_DRCut_3Jet_c->Print((OutputPlotDir+"/"+DEta_GroomedCA8_GenWhad_DRCut_3Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedCA8_GenWhad_DRCut_3Jet_c    = new TCanvas (DPhi_GroomedCA8_GenWhad_DRCut_3Jet->GetName());
  DPhi_GroomedCA8_GenWhad_DRCut_3Jet_c->cd();
  DPhi_GroomedCA8_GenWhad_DRCut_3Jet_c->SetLogy();
  DPhi_GroomedCA8_GenWhad_DRCut_3Jet->Draw();
  DPhi_GroomedCA8_GenWhad_DRCut_3Jet->Write();
  DPhi_GroomedCA8_GenWhad_DRCut_3Jet_c->Print((OutputPlotDir+"/"+DPhi_GroomedCA8_GenWhad_DRCut_3Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedCA8_GenWhad_DRCut_3Jet_c    = new TCanvas (DR_GroomedCA8_GenWhad_DRCut_3Jet->GetName());
  DR_GroomedCA8_GenWhad_DRCut_3Jet_c->cd();
  DR_GroomedCA8_GenWhad_DRCut_3Jet_c->SetLogy();
  DR_GroomedCA8_GenWhad_DRCut_3Jet->Draw();
  DR_GroomedCA8_GenWhad_DRCut_3Jet->Write();
  DR_GroomedCA8_GenWhad_DRCut_3Jet_c->Print((OutputPlotDir+"/"+DR_GroomedCA8_GenWhad_DRCut_3Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedCA8_GenWhad_DRCut_3Jet_c    = new TCanvas (DMass_GroomedCA8_GenWhad_DRCut_3Jet->GetName());
  DMass_GroomedCA8_GenWhad_DRCut_3Jet_c->cd();
  DMass_GroomedCA8_GenWhad_DRCut_3Jet_c->SetLogy();
  DMass_GroomedCA8_GenWhad_DRCut_3Jet->Draw();
  DMass_GroomedCA8_GenWhad_DRCut_3Jet->Write();
  DMass_GroomedCA8_GenWhad_DRCut_3Jet_c->Print((OutputPlotDir+"/"+DMass_GroomedCA8_GenWhad_DRCut_3Jet->GetName()+".png").c_str(),"png");


  // AK5 After exclusive selection of events with only one AK5 Jet (DR)

  TCanvas*  DPt_GroomedAK5_GenWhad_DRCut_3Jet_c    = new TCanvas (DPt_GroomedAK5_GenWhad_DRCut_3Jet->GetName());
  DPt_GroomedAK5_GenWhad_DRCut_3Jet_c->cd();
  DPt_GroomedAK5_GenWhad_DRCut_3Jet_c->SetLogy();
  DPt_GroomedAK5_GenWhad_DRCut_3Jet->Draw();
  DPt_GroomedAK5_GenWhad_DRCut_3Jet->Write();
  DPt_GroomedAK5_GenWhad_DRCut_3Jet_c->Print((OutputPlotDir+"/"+DPt_GroomedAK5_GenWhad_DRCut_3Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedAK5_GenWhad_DRCut_3Jet_c    = new TCanvas (DEta_GroomedAK5_GenWhad_DRCut_3Jet->GetName());
  DEta_GroomedAK5_GenWhad_DRCut_3Jet_c->cd();
  DEta_GroomedAK5_GenWhad_DRCut_3Jet_c->SetLogy();
  DEta_GroomedAK5_GenWhad_DRCut_3Jet->Draw();
  DEta_GroomedAK5_GenWhad_DRCut_3Jet->Write();
  DEta_GroomedAK5_GenWhad_DRCut_3Jet_c->Print((OutputPlotDir+"/"+DEta_GroomedAK5_GenWhad_DRCut_3Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedAK5_GenWhad_DRCut_3Jet_c    = new TCanvas (DPhi_GroomedAK5_GenWhad_DRCut_3Jet->GetName());
  DPhi_GroomedAK5_GenWhad_DRCut_3Jet_c->cd();
  DPhi_GroomedAK5_GenWhad_DRCut_3Jet_c->SetLogy();
  DPhi_GroomedAK5_GenWhad_DRCut_3Jet->Draw();
  DPhi_GroomedAK5_GenWhad_DRCut_3Jet->Write();
  DPhi_GroomedAK5_GenWhad_DRCut_3Jet_c->Print((OutputPlotDir+"/"+DPhi_GroomedAK5_GenWhad_DRCut_3Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedAK5_GenWhad_DRCut_3Jet_c    = new TCanvas (DR_GroomedAK5_GenWhad_DRCut_3Jet->GetName());
  DR_GroomedAK5_GenWhad_DRCut_3Jet_c->cd();
  DR_GroomedAK5_GenWhad_DRCut_3Jet_c->SetLogy();
  DR_GroomedAK5_GenWhad_DRCut_3Jet->Draw();
  DR_GroomedAK5_GenWhad_DRCut_3Jet->Write();
  DR_GroomedAK5_GenWhad_DRCut_3Jet_c->Print((OutputPlotDir+"/"+DR_GroomedAK5_GenWhad_DRCut_3Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedAK5_GenWhad_DRCut_3Jet_c    = new TCanvas (DMass_GroomedAK5_GenWhad_DRCut_3Jet->GetName());
  DMass_GroomedAK5_GenWhad_DRCut_3Jet_c->cd();
  DMass_GroomedAK5_GenWhad_DRCut_3Jet_c->SetLogy();
  DMass_GroomedAK5_GenWhad_DRCut_3Jet->Draw();
  DMass_GroomedAK5_GenWhad_DRCut_3Jet->Write();
  DMass_GroomedAK5_GenWhad_DRCut_3Jet_c->Print((OutputPlotDir+"/"+DMass_GroomedAK5_GenWhad_DRCut_3Jet->GetName()+".png").c_str(),"png");


  // CA8 After exclusive selection of events with only two CA8 Jet (DR)

  TCanvas*  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_c    = new TCanvas (DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet->GetName());
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_c->cd();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_c->SetLogy();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet->Draw();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet->Write();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_c->Print((OutputPlotDir+"/"+DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_c    = new TCanvas (DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet->GetName());
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_c->cd();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_c->SetLogy();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet->Draw();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet->Write();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_c->Print((OutputPlotDir+"/"+DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_c    = new TCanvas (DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet->GetName());
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_c->cd();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_c->SetLogy();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet->Draw();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet->Write();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_c->Print((OutputPlotDir+"/"+DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_c    = new TCanvas (DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet->GetName());
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_c->cd();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_c->SetLogy();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet->Draw();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet->Write();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_c->Print((OutputPlotDir+"/"+DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_c    = new TCanvas (DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet->GetName());
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_c->cd();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_c->SetLogy();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet->Draw();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet->Write();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_c->Print((OutputPlotDir+"/"+DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_c    = new TCanvas (DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet->GetName());
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_c->cd();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_c->SetLogy();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet->Draw();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet->Write();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_c->Print((OutputPlotDir+"/"+DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_c    = new TCanvas (DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet->GetName());
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_c->cd();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_c->SetLogy();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet->Draw();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet->Write();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_c->Print((OutputPlotDir+"/"+DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_c    = new TCanvas (DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet->GetName());
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_c->cd();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_c->SetLogy();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet->Draw();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet->Write();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_c->Print((OutputPlotDir+"/"+DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_c    = new TCanvas (DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet->GetName());
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_c->cd();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_c->SetLogy();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet->Draw();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet->Write();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_c->Print((OutputPlotDir+"/"+DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_c    = new TCanvas (DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet->GetName());
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_c->cd();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_c->SetLogy();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet->Draw();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet->Write();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_c->Print((OutputPlotDir+"/"+DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet->GetName()+".png").c_str(),"png");



  // CA8 After exclusive selection of events with only three CA8 Jet (DR)

  TCanvas*  DPt_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching_c    = new TCanvas (DPt_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->GetName());
  DPt_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching_c->cd();
  DPt_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching_c->SetLogy();
  DPt_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->Draw();
  DPt_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->Write();
  DPt_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DPt_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching_c    = new TCanvas (DEta_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->GetName());
  DEta_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching_c->cd();
  DEta_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching_c->SetLogy();
  DEta_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->Draw();
  DEta_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->Write();
  DEta_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DEta_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching_c    = new TCanvas (DPhi_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->GetName());
  DPhi_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching_c->cd();
  DPhi_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching_c->SetLogy();
  DPhi_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->Draw();
  DPhi_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->Write();
  DPhi_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DPhi_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching_c    = new TCanvas (DR_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->GetName());
  DR_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching_c->cd();
  DR_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching_c->SetLogy();
  DR_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->Draw();
  DR_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->Write();
  DR_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DR_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching_c    = new TCanvas (DMass_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->GetName());
  DMass_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching_c->cd();
  DMass_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching_c->SetLogy();
  DMass_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->Draw();
  DMass_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->Write();
  DMass_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DMass_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->GetName()+".png").c_str(),"png");


  // AK5 After exclusive selection of events with only one AK5 Jet (DR)

  TCanvas*  DPt_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching_c    = new TCanvas (DPt_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->GetName());
  DPt_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching_c->cd();
  DPt_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching_c->SetLogy();
  DPt_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->Draw();
  DPt_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->Write();
  DPt_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DPt_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching_c    = new TCanvas (DEta_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->GetName());
  DEta_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching_c->cd();
  DEta_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching_c->SetLogy();
  DEta_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->Draw();
  DEta_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->Write();
  DEta_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DEta_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching_c    = new TCanvas (DPhi_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->GetName());
  DPhi_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching_c->cd();
  DPhi_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching_c->SetLogy();
  DPhi_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->Draw();
  DPhi_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->Write();
  DPhi_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DPhi_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching_c    = new TCanvas (DR_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->GetName());
  DR_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching_c->cd();
  DR_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching_c->SetLogy();
  DR_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->Draw();
  DR_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->Write();
  DR_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DR_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching_c    = new TCanvas (DMass_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->GetName());
  DMass_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching_c->cd();
  DMass_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching_c->SetLogy();
  DMass_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->Draw();
  DMass_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->Write();
  DMass_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DMass_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->GetName()+".png").c_str(),"png");


  // CA8 After exclusive selection of events with only two CA8 Jet (DR)

  TCanvas*  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c    = new TCanvas (DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetName());
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->cd();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->SetLogy();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->Draw();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->Write();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c    = new TCanvas (DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetName());
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->cd();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->SetLogy();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->Draw();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->Write();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c    = new TCanvas (DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetName());
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->cd();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->SetLogy();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->Draw();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->Write();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c    = new TCanvas (DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetName());
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->cd();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->SetLogy();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->Draw();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->Write();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c    = new TCanvas (DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetName());
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->cd();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->SetLogy();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->Draw();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->Write();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c    = new TCanvas (DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetName());
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->cd();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->SetLogy();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->Draw();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->Write();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c    = new TCanvas (DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetName());
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->cd();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->SetLogy();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->Draw();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->Write();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c    = new TCanvas (DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetName());
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->cd();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->SetLogy();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->Draw();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->Write();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c    = new TCanvas (DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetName());
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->cd();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->SetLogy();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->Draw();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->Write();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c    = new TCanvas (DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetName());
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->cd();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->SetLogy();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->Draw();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->Write();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching_c->Print((OutputPlotDir+"/"+DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetName()+".png").c_str(),"png");

  // CA8 After exclusive selection of events with only three CA8 Jet (DR)

  TCanvas*  DPt_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching_c    = new TCanvas (DPt_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->GetName());
  DPt_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching_c->cd();
  DPt_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching_c->SetLogy();
  DPt_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->Draw();
  DPt_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->Write();
  DPt_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DPt_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching_c    = new TCanvas (DEta_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->GetName());
  DEta_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching_c->cd();
  DEta_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching_c->SetLogy();
  DEta_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->Draw();
  DEta_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->Write();
  DEta_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DEta_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching_c    = new TCanvas (DPhi_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->GetName());
  DPhi_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching_c->cd();
  DPhi_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching_c->SetLogy();
  DPhi_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->Draw();
  DPhi_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->Write();
  DPhi_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DPhi_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching_c    = new TCanvas (DR_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->GetName());
  DR_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching_c->cd();
  DR_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching_c->SetLogy();
  DR_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->Draw();
  DR_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->Write();
  DR_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DR_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching_c    = new TCanvas (DMass_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->GetName());
  DMass_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching_c->cd();
  DMass_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching_c->SetLogy();
  DMass_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->Draw();
  DMass_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->Write();
  DMass_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DMass_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->GetName()+".png").c_str(),"png");


  // AK5 After exclusive selection of events with only one AK5 Jet (DR)

  TCanvas*  DPt_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching_c    = new TCanvas (DPt_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->GetName());
  DPt_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching_c->cd();
  DPt_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching_c->SetLogy();
  DPt_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->Draw();
  DPt_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->Write();
  DPt_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DPt_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching_c    = new TCanvas (DEta_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->GetName());
  DEta_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching_c->cd();
  DEta_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching_c->SetLogy();
  DEta_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->Draw();
  DEta_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->Write();
  DEta_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DEta_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching_c    = new TCanvas (DPhi_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->GetName());
  DPhi_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching_c->cd();
  DPhi_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching_c->SetLogy();
  DPhi_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->Draw();
  DPhi_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->Write();
  DPhi_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DPhi_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching_c    = new TCanvas (DR_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->GetName());
  DR_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching_c->cd();
  DR_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching_c->SetLogy();
  DR_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->Draw();
  DR_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->Write();
  DR_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DR_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching_c    = new TCanvas (DMass_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->GetName());
  DMass_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching_c->cd();
  DMass_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching_c->SetLogy();
  DMass_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->Draw();
  DMass_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->Write();
  DMass_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DMass_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->GetName()+".png").c_str(),"png");


  // CA8 After exclusive selection of events with only two CA8 Jet (DR)

  TCanvas*  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c    = new TCanvas (DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetName());
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->cd();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->SetLogy();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->Draw();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->Write();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c    = new TCanvas (DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetName());
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->cd();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->SetLogy();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->Draw();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->Write();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c    = new TCanvas (DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetName());
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->cd();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->SetLogy();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->Draw();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->Write();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c    = new TCanvas (DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetName());
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->cd();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->SetLogy();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->Draw();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->Write();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c    = new TCanvas (DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetName());
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->cd();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->SetLogy();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->Draw();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->Write();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c    = new TCanvas (DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetName());
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->cd();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->SetLogy();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->Draw();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->Write();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c    = new TCanvas (DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetName());
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->cd();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->SetLogy();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->Draw();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->Write();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c    = new TCanvas (DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetName());
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->cd();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->SetLogy();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->Draw();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->Write();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c    = new TCanvas (DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetName());
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->cd();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->SetLogy();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->Draw();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->Write();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c    = new TCanvas (DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetName());
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->cd();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->SetLogy();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->Draw();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->Write();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching_c->Print((OutputPlotDir+"/"+DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetName()+".png").c_str(),"png");

  // CA8 After exclusive selection of events with only three CA8 Jet (DR)

  TCanvas*  DPt_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching_c    = new TCanvas (DPt_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->GetName());
  DPt_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching_c->cd();
  DPt_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching_c->SetLogy();
  DPt_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->Draw();
  DPt_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->Write();
  DPt_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching_c->Print((OutputPlotDir+"/"+DPt_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching_c    = new TCanvas (DEta_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->GetName());
  DEta_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching_c->cd();
  DEta_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching_c->SetLogy();
  DEta_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->Draw();
  DEta_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->Write();
  DEta_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching_c->Print((OutputPlotDir+"/"+DEta_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching_c    = new TCanvas (DPhi_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->GetName());
  DPhi_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching_c->cd();
  DPhi_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching_c->SetLogy();
  DPhi_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->Draw();
  DPhi_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->Write();
  DPhi_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching_c->Print((OutputPlotDir+"/"+DPhi_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching_c    = new TCanvas (DR_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->GetName());
  DR_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching_c->cd();
  DR_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching_c->SetLogy();
  DR_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->Draw();
  DR_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->Write();
  DR_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching_c->Print((OutputPlotDir+"/"+DR_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching_c    = new TCanvas (DMass_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->GetName());
  DMass_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching_c->cd();
  DMass_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching_c->SetLogy();
  DMass_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->Draw();
  DMass_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->Write();
  DMass_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching_c->Print((OutputPlotDir+"/"+DMass_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->GetName()+".png").c_str(),"png");


  // AK5 After exclusive selection of events with only one AK5 Jet (DR)

  TCanvas*  DPt_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching_c    = new TCanvas (DPt_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->GetName());
  DPt_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching_c->cd();
  DPt_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching_c->SetLogy();
  DPt_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->Draw();
  DPt_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->Write();
  DPt_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching_c->Print((OutputPlotDir+"/"+DPt_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching_c    = new TCanvas (DEta_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->GetName());
  DEta_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching_c->cd();
  DEta_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching_c->SetLogy();
  DEta_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->Draw();
  DEta_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->Write();
  DEta_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching_c->Print((OutputPlotDir+"/"+DEta_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching_c    = new TCanvas (DPhi_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->GetName());
  DPhi_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching_c->cd();
  DPhi_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching_c->SetLogy();
  DPhi_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->Draw();
  DPhi_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->Write();
  DPhi_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching_c->Print((OutputPlotDir+"/"+DPhi_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching_c    = new TCanvas (DR_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->GetName());
  DR_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching_c->cd();
  DR_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching_c->SetLogy();
  DR_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->Draw();
  DR_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->Write();
  DR_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching_c->Print((OutputPlotDir+"/"+DR_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching_c    = new TCanvas (DMass_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->GetName());
  DMass_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching_c->cd();
  DMass_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching_c->SetLogy();
  DMass_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->Draw();
  DMass_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->Write();
  DMass_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching_c->Print((OutputPlotDir+"/"+DMass_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->GetName()+".png").c_str(),"png");


  // CA8 After exclusive selection of events with only two CA8 Jet (DR)

  TCanvas*  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c    = new TCanvas (DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetName());
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->cd();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->SetLogy();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->Draw();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->Write();
  DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->Print((OutputPlotDir+"/"+DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c    = new TCanvas (DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetName());
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->cd();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->SetLogy();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->Draw();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->Write();
  DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->Print((OutputPlotDir+"/"+DEta_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c    = new TCanvas (DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetName());
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->cd();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->SetLogy();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->Draw();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->Write();
  DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->Print((OutputPlotDir+"/"+DPhi_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c    = new TCanvas (DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetName());
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->cd();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->SetLogy();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->Draw();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->Write();
  DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->Print((OutputPlotDir+"/"+DR_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c    = new TCanvas (DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetName());
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->cd();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->SetLogy();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->Draw();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->Write();
  DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->Print((OutputPlotDir+"/"+DMass_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c    = new TCanvas (DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetName());
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->cd();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->SetLogy();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->Draw();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->Write();
  DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->Print((OutputPlotDir+"/"+DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c    = new TCanvas (DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetName());
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->cd();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->SetLogy();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->Draw();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->Write();
  DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->Print((OutputPlotDir+"/"+DEta_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c    = new TCanvas (DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetName());
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->cd();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->SetLogy();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->Draw();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->Write();
  DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->Print((OutputPlotDir+"/"+DPhi_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c    = new TCanvas (DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetName());
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->cd();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->SetLogy();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->Draw();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->Write();
  DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->Print((OutputPlotDir+"/"+DR_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetName()+".png").c_str(),"png");

  TCanvas*  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c    = new TCanvas (DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetName());
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->cd();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->SetLogy();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->Draw();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->Write();
  DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching_c->Print((OutputPlotDir+"/"+DMass_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetName()+".png").c_str(),"png");




  // Print Out the measured Efficiency

  std::cout<<" ########################################################################### "<<std::endl;
  std::cout<<" ###########  Measured Efficiency Applying DR + DPt  Selection ############# "<<std::endl;
  std::cout<<" ########################################################################### "<<std::endl;

  std::cout<<std::endl;
  std::cout<< " ################################## " <<std::endl;
  std::cout<< " Groomed CA8 Jet vs Gen W inclusive " <<std::endl;
  std::cout<< " ################################## " <<std::endl;
  std::cout<<std::endl;
      
  for(int i =0; i<NJet ; i++){

    std::cout<<" Tree Event : "<<fTree->GetEntries()<<" Jet Index : "<<i<<" Boosted Cut + Pt > 30 : "<<DPt_GroomedJet_CA8_GenWhad[i]->Integral(0,DPt_GroomedJet_CA8_GenWhad[i]->GetNbinsX()+1)<<
    " DR Selection : "<<DPt_GroomedCA8_GenWhad_DRCut[i]->Integral(0,DPt_GroomedCA8_GenWhad_DRCut[i]->GetNbinsX()+1)<<" Matching Efficiency : "<<
    DPt_GroomedCA8_GenWhad_DRCut[i]->Integral(0,DPt_GroomedCA8_GenWhad_DRCut[i]->GetNbinsX())/DPt_GroomedJet_CA8_GenWhad[i]->Integral(0,DPt_GroomedJet_CA8_GenWhad[i]->GetNbinsX()+1)<<
    std::endl;

    std::cout<<" Tree Event : "<<fTree->GetEntries()<<" Jet Index : "<<i<<" Boosted Cut + Pt > 30 : "<<DPt_GroomedJet_CA8_GenWhad[i]->Integral(0,DPt_GroomedJet_CA8_GenWhad[i]->GetNbinsX()+1)<<
    " DR + DPt Selection : "<<DPt_GroomedCA8_GenWhad_DRCut_DPtCut[i]->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_DPtCut[i]->GetNbinsX()+1)<<" Matching Efficiency : "<<
    DPt_GroomedCA8_GenWhad_DRCut_DPtCut[i]->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_DPtCut[i]->GetNbinsX()+1)/DPt_GroomedJet_CA8_GenWhad[i]->
    Integral(0,DPt_GroomedJet_CA8_GenWhad[i]->GetNbinsX()+1)<<std::endl;

  }

  std::cout<<std::endl;
  std::cout<< " ################################## " <<std::endl;
  std::cout<< " Groomed AK5 Jet vs Gen W inclusive " <<std::endl;
  std::cout<< " ################################## " <<std::endl;
  std::cout<<std::endl;
    
  for(int i =0; i<NJet ; i++){

    std::cout<<" Tree Event : "<<fTree->GetEntries()<<" Jet Index : "<<i<<" Boosted Cut + Pt > 30 : "<<DPt_GroomedJet_AK5_GenWhad[i]->Integral(0,DPt_GroomedJet_AK5_GenWhad[i]->GetNbinsX()+1)<<
    " DR Selection : "<<DPt_GroomedAK5_GenWhad_DRCut[i]->Integral(0,DPt_GroomedAK5_GenWhad_DRCut[i]->GetNbinsX()+1)<<" Matching Efficiency : "<<
    DPt_GroomedAK5_GenWhad_DRCut[i]->Integral(0,DPt_GroomedAK5_GenWhad_DRCut[i]->GetNbinsX()+1)/DPt_GroomedJet_AK5_GenWhad[i]->Integral(0,DPt_GroomedJet_AK5_GenWhad[i]->GetNbinsX()+1)<<
    std::endl;

    std::cout<<" Tree Event : "<<fTree->GetEntries()<<" Jet Index : "<<i<<" Boosted Cut + Pt > 30 : "<<DPt_GroomedJet_AK5_GenWhad[i]->Integral(0,DPt_GroomedJet_AK5_GenWhad[i]->GetNbinsX()+1)
    <<" DR + DPt Selection : "<<DPt_GroomedAK5_GenWhad_DRCut_DPtCut[i]->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_DPtCut[i]->GetNbinsX()+1)<<" Matching Efficiency : "<<
    DPt_GroomedAK5_GenWhad_DRCut_DPtCut[i]->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_DPtCut[i]->GetNbinsX()+1)/DPt_GroomedJet_AK5_GenWhad[i]->
    Integral(0,DPt_GroomedJet_AK5_GenWhad[i]->GetNbinsX()+1)<<std::endl;

  }

  std::cout<<std::endl;
  std::cout<< " ######################################################################## " <<std::endl;
  std::cout<< " Groomed CA8 Jet vs Gen W Exclusive for events with only one hard CA8 Jet " <<std::endl;
  std::cout<< " ######################################################################## " <<std::endl;
  std::cout<<std::endl;
      
  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<0<<" Excluisive NumEvent : "<<nExclusive_1Jet_CA8<<
    " DR Selection (Only Hard Jet Matched) : "<<DPt_GroomedCA8_GenWhad_DRCut_1Jet->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_1Jet->GetNbinsX()+1)<<" Matching Efficiency : "<<
    DPt_GroomedCA8_GenWhad_DRCut_1Jet->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_1Jet->GetNbinsX()+1)/nExclusive_1Jet_CA8<<std::endl;

  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<0<<" Excluisive NumEvent : "<<nExclusive_1Jet_CA8<<
  " DR + DPt Selection (Only Hard Jet Matched) : "<<DPt_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet->GetNbinsX()+1)<<" Matching Efficiency : "<<
    DPt_GroomedCA8_GenWhad_DRCut_DPtCut_1Jet->GetEffectiveEntries()/nExclusive_1Jet_CA8<<std::endl;

  std::cout<<std::endl;
  std::cout<< " ######################################################################## " <<std::endl;
  std::cout<< " Groomed AK5 Jet vs Gen W Exclusive for events with only one hard CA8 Jet" <<std::endl;
  std::cout<< " ######################################################################## " <<std::endl;
  std::cout<<std::endl;
      
  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<0<<" Excluisive NumEvent : "<<nExclusive_1Jet_CA8<<
  " DR Selection (Only Hard Jet Matched) : "<<DPt_GroomedAK5_GenWhad_DRCut_1Jet->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_1Jet->GetNbinsX()+1)<<" Matching Efficiency : "<<
    DPt_GroomedAK5_GenWhad_DRCut_1Jet->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_1Jet->GetNbinsX()+1)/nExclusive_1Jet_CA8<<std::endl;

  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<0<<" Excluisive NumEvent : "<<nExclusive_1Jet_CA8<<
  " DR + DPt Selection (Only Hard Jet Matched) : "<<DPt_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->GetNbinsX()+1)<<" Matching Efficiency : "<<
    DPt_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_DPtCut_1Jet->GetNbinsX()+1)/nExclusive_1Jet_CA8<<std::endl;

  std::cout<<std::endl;
  std::cout<< " ######################################################################## " <<std::endl;
  std::cout<< " Groomed CA8 Jet vs Gen W Exclusive for events with only two hard CA8 Jet" <<std::endl;
  std::cout<< " ######################################################################## " <<std::endl; 
  std::cout<<std::endl;
      
  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<0<<" Excluisive NumEvent : "<<nExclusive_2Jet_CA8<<
  " DR Selection (Only Hard Jet Matched) : "<<DPt_GroomedCA8_GenWhad_DRCut_2Jet->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_2Jet->GetNbinsX()+1)<<" Matching Efficiency : "<<
    DPt_GroomedCA8_GenWhad_DRCut_2Jet->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_2Jet->GetNbinsX()+1)/nExclusive_2Jet_CA8<<std::endl;
  
  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<0<<" Excluisive NumEvent : "<<nExclusive_2Jet_CA8<<
  " DR + DPt Selection (Only Hard Jet Matched) : "<<DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->GetNbinsX()+1)<<" Matching Efficiency : "<<
    DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->GetNbinsX()+1)/nExclusive_2Jet_CA8<<std::endl;
  
  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<1<<" Excluisive NumEvent : "<<nExclusive_2Jet_CA8<<
    " DR Selection (Double Matching) : "<<DPt_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->GetNbinsX()+1)<<
    " Matching Efficiency : "<<DPt_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_2Jet_DoubleMatching->GetNbinsX()+1)/nExclusive_2Jet_CA8<<std::endl;

  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<1<<" Excluisive NumEvent : "<<nExclusive_2Jet_CA8<<
    " DR + DPt Selection (Double Matching) : "<<DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetNbinsX()+1)<<
    " Matching Efficiency : "<<DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetNbinsX()+1)/nExclusive_2Jet_CA8<<
  std::endl;

  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<1<<" Excluisive NumEvent : "<<nExclusive_2Jet_CA8<<
    " DR Selection (Inverted Matching) : "<<DPt_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->GetNbinsX()+1)<<
    " Matching Efficiency : "<<DPt_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_2Jet_InvertedMatching->GetNbinsX()+1)/nExclusive_2Jet_CA8<<
  std::endl;

  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<1<<" Excluisive NumEvent : "<<nExclusive_2Jet_CA8<<
  " DR + DPt Selection (Inverted Matching) : "<<DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetNbinsX()+1)
  <<" Matching Efficiency : "<<DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetNbinsX()+1)/
  nExclusive_2Jet_CA8<<std::endl;

  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<1<<" Excluisive NumEvent : "<<nExclusive_2Jet_CA8<<
  " DR Selection (No Matching) : "<<DPt_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->GetNbinsX()+1)<<" Matching Efficiency : "<<
    DPt_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_2Jet_NoMatching->GetNbinsX()+1)/nExclusive_2Jet_CA8<<std::endl;

  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<1<<" Excluisive NumEvent : "<<nExclusive_2Jet_CA8<<
    " DR + DPt Selection (No Matching) : "<<DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetNbinsX()+1)<<
    " Matching Efficiency : "<<DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet_NoMatching->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_DPtCut_2Jet->GetNbinsX()+1)/nExclusive_2Jet_CA8<<std::endl;

  std::cout<<std::endl;
  std::cout<< " ######################################################################## " <<std::endl;
  std::cout<< " Groomed AK5 Jet vs Gen W Exclusive for events with only two hard AK5 Jet" <<std::endl;
  std::cout<< " ######################################################################## " <<std::endl; 
  std::cout<<std::endl;
      
  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<0<<" Excluisive NumEvent : "<<nExclusive_2Jet_CA8<<
  " DR Selection (Only Hard Jet Matched) : "<<DPt_GroomedAK5_GenWhad_DRCut_2Jet->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_2Jet->GetNbinsX()+1)<<" Matching Efficiency : "<<
    DPt_GroomedAK5_GenWhad_DRCut_2Jet->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_2Jet->GetNbinsX()+1)/nExclusive_2Jet_CA8<<std::endl;

  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<0<<" Excluisive NumEvent : "<<nExclusive_2Jet_CA8<<
  " DR + DPt Selection (Only Hard Jet Matched) : "<<DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->GetNbinsX()+1)<<" Matching Efficiency : "<<
    DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->GetNbinsX()+1)/nExclusive_2Jet_CA8<<std::endl;

  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<1<<" Excluisive NumEvent : "<<nExclusive_2Jet_CA8<<
    " DR Selection (Double Matching) : "<<DPt_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->GetNbinsX()+1)<<
    " Matching Efficiency : "<<DPt_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_2Jet_DoubleMatching->GetNbinsX()+1)/nExclusive_2Jet_CA8<<std::endl;

  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<1<<" Excluisive NumEvent : "<<nExclusive_2Jet_CA8<<
    " DR + DPt Selection (Double Matching) : "<<DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetNbinsX()+1)<<
    " Matching Efficiency : "<<DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_DoubleMatching->GetNbinsX()+1)/nExclusive_2Jet_CA8<<
  std::endl;

  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<1<<" Excluisive NumEvent : "<<nExclusive_2Jet_CA8<<
    " DR Selection (Inverted Matching) : "<<DPt_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->GetNbinsX()+1)<<
    " Matching Efficiency : "<<DPt_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_2Jet_InvertedMatching->GetNbinsX()+1)/nExclusive_2Jet_CA8<<
  std::endl;

  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<1<<" Excluisive NumEvent : "<<nExclusive_2Jet_CA8<<
  " DR + DPt Selection (Inverted Matching) : "<<DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetNbinsX()+1)
   <<" Matching Efficiency : "<<DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_InvertedMatching->GetNbinsX()+1)/
   nExclusive_2Jet_CA8<<std::endl;

  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<1<<" Excluisive NumEvent : "<<nExclusive_2Jet_CA8<<
    " DR Selection (No Matching) : "<<DPt_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->GetNbinsX()+1)<<" Matching Efficiency : "<<
    DPt_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_2Jet_NoMatching->GetNbinsX()+1)/nExclusive_2Jet_CA8<<std::endl;

  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<1<<" Excluisive NumEvent : "<<nExclusive_2Jet_CA8<<
    " DR + DPt Selection (No Matching) : "<<DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching->GetNbinsX()+1)<<
    " Matching Efficiency : "<<DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet_NoMatching->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_DPtCut_2Jet->GetNbinsX()+1)/nExclusive_2Jet_CA8<<std::endl;


  std::cout<<std::endl;
  std::cout<< " ######################################################################## " <<std::endl;
  std::cout<< " Groomed CA8 Jet vs Gen W Exclusive for events with only three hard CA8 Jet" <<std::endl;
  std::cout<< " ######################################################################## " <<std::endl; 
  std::cout<<std::endl;
      
  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<0<<" Excluisive NumEvent : "<<nExclusive_3Jet_CA8<<
  " DR Selection (Only Hard Jet Matched) : "<<DPt_GroomedCA8_GenWhad_DRCut_3Jet->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_3Jet->GetNbinsX()+1)<<" Matching Efficiency : "<<
    DPt_GroomedCA8_GenWhad_DRCut_3Jet->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_3Jet->GetNbinsX()+1)/nExclusive_3Jet_CA8<<std::endl;

  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<0<<" Excluisive NumEvent : "<<nExclusive_3Jet_CA8<<
  " DR + DPt Selection (Only Hard Jet Matched) : "<<DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet->GetNbinsX()+1)<<" Matching Efficiency : "<<
    DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet->GetEffectiveEntries()/nExclusive_3Jet_CA8<<std::endl;

  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<1<<" Excluisive NumEvent : "<<nExclusive_3Jet_CA8<<
  " DR Selection (Double Matching) : "<<DPt_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->GetNbinsX()+1)<<
  " Matching Efficiency : "<<DPt_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_3Jet_DoubleMatching->GetNbinsX()+1)/nExclusive_3Jet_CA8<<std::endl;

  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<1<<" Excluisive NumEvent : "<<nExclusive_3Jet_CA8<<
  " DR + DPt (Double Matching) : "<<DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetNbinsX()+1)<<
  " Matching Efficiency : "<<DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetNbinsX()+1)/nExclusive_3Jet_CA8
  <<std::endl;

  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<1<<" Excluisive NumEvent : "<<nExclusive_3Jet_CA8<<
  " DR Selection (Inverted Matching) : "<<DPt_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->GetNbinsX()+1)<<
  " Matching Efficiency : "<<DPt_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_3Jet_InvertedMatching->GetNbinsX()+1)/nExclusive_3Jet_CA8<<std::endl;

  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<1<<" Excluisive NumEvent : "<<nExclusive_3Jet_CA8<<
  " DR + DPt (Inverted Matching) : "<<DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetNbinsX()+1)<<
  " Matching Efficiency : "<<DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetNbinsX()+1)/nExclusive_3Jet_CA8
  <<std::endl;

  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<1<<" Excluisive NumEvent : "<<nExclusive_3Jet_CA8<<
  " DR Selection (No Matching) : "<<DPt_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->GetNbinsX()+1)<<" Matching Efficiency : "<<
    DPt_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_3Jet_NoMatching->GetNbinsX()+1)/nExclusive_3Jet_CA8<<std::endl;

  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<1<<" Excluisive NumEvent : "<<nExclusive_3Jet_CA8<<
  " DR + DPt (No Matching) : "<<DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetNbinsX()+1)<<
  " Matching Efficiency : "<<DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->Integral(0,DPt_GroomedCA8_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetNbinsX()+1)/nExclusive_3Jet_CA8<<std::endl;


  std::cout<<std::endl;
  std::cout<< " ######################################################################## " <<std::endl;
  std::cout<< " Groomed AK5 Jet vs Gen W Exclusive for events with only three hard CA8 Jet" <<std::endl;
  std::cout<< " ######################################################################## " <<std::endl; 
  std::cout<<std::endl;
      
  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<0<<" Excluisive NumEvent : "<<nExclusive_3Jet_CA8<<
  " DR Selection (Only Hard Jet Matched) : "<<DPt_GroomedAK5_GenWhad_DRCut_3Jet->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_3Jet->GetNbinsX()+1)<<" Matching Efficiency : "<<
    DPt_GroomedAK5_GenWhad_DRCut_3Jet->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_3Jet->GetNbinsX()+1)/nExclusive_3Jet_CA8<<std::endl;

  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<0<<" Excluisive NumEvent : "<<nExclusive_3Jet_CA8<<
  " DR + DPt Selection (Only Hard Jet Matched) : "<<DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet->GetNbinsX()+1)<<" Matching Efficiency : "<<
    DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet->GetEffectiveEntries()/nExclusive_3Jet_CA8<<std::endl;

  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<1<<" Excluisive NumEvent : "<<nExclusive_3Jet_CA8<<
  " DR Selection (Double Matching) : "<<DPt_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->GetNbinsX()+1)<<
  " Matching Efficiency : "<<DPt_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_3Jet_DoubleMatching->GetNbinsX()+1)/nExclusive_3Jet_CA8<<std::endl;

  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<1<<" Excluisive NumEvent : "<<nExclusive_3Jet_CA8<<
  " DR + DPt (Double Matching) : "<<DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetNbinsX()+1)<<
  " Matching Efficiency : "<<DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_DoubleMatching->GetNbinsX()+1)/nExclusive_3Jet_CA8
  <<std::endl;

  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<1<<" Excluisive NumEvent : "<<nExclusive_3Jet_CA8<<
  " DR Selection (Inverted Matching) : "<<DPt_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->GetNbinsX()+1)<<
  " Matching Efficiency : "<<DPt_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_3Jet_InvertedMatching->GetNbinsX()+1)/nExclusive_3Jet_CA8<<std::endl;

  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<1<<" Excluisive NumEvent : "<<nExclusive_3Jet_CA8<<
  " DR + DPt (Inverted Matching) : "<<DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetNbinsX()+1)<<
  " Matching Efficiency : "<<DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_InvertedMatching->GetNbinsX()+1)/nExclusive_3Jet_CA8
  <<std::endl;

  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<1<<" Excluisive NumEvent : "<<nExclusive_3Jet_CA8<<
  " DR Selection (No Matching) : "<<DPt_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->GetNbinsX()+1)<<" Matching Efficiency : "<<
    DPt_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_3Jet_NoMatching->GetNbinsX()+1)/nExclusive_3Jet_CA8<<std::endl;

  std::cout<<" Tree Event Number : "<<fTree->GetEntries()<<" Jet Index : "<<1<<" Excluisive NumEvent : "<<nExclusive_3Jet_CA8<<
  " DR + DPt (No Matching) : "<<DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetNbinsX()+1)<<
  " Matching Efficiency : "<<DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->Integral(0,DPt_GroomedAK5_GenWhad_DRCut_DPtCut_3Jet_NoMatching->GetNbinsX()+1)/nExclusive_3Jet_CA8<<std::endl;

  
  std::cout<<std::endl;
  outputFile->Close();

  return 0;


}

