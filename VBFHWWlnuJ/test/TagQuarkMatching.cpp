#include <iostream>
#include <vector>
#include <algorithm>
#include <iostream>

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

  // Load TTree Lybrary                                                                                                                                                                          
  gSystem->Load("libTree.so");

  // Set Root style from global enviroment path                                                                                                                                                  
  std::string ROOTStyle =  getenv ("ROOTStyle");

  gROOT->ProcessLine((".x "+ROOTStyle+"/rootLogon.C").c_str());
  gROOT->ProcessLine((".x "+ROOTStyle+"/rootPalette.C").c_str());
  gROOT->ProcessLine((".x "+ROOTStyle+"/rootColors.C").c_str());
  gROOT->ProcessLine((".x "+ROOTStyle+"/setTDRStyle.C").c_str());

  // parse config file parameter                                                                                                                                                                 
  parseConfigFile(argv[1]);

  std::string InputDirectory     = gConfigParser -> readStringOption("Input::InputDirectory");
  std::string InputFileName      = gConfigParser -> readStringOption("Input::InputFileName");

  std::string OutputDirectory     = gConfigParser -> readStringOption("Output::OutputDirectory");
  std::string OutputFileName      = gConfigParser -> readStringOption("Output::OutputFileName");

  std::cout<<"         "<<std::endl;
  std::cout<<" InputDirectory :     "<<InputDirectory<<std::endl;
  std::cout<<"         "<<std::endl;

  std::cout<<" InputFileName :      "<<InputFileName<<std::endl;
  std::cout<<"         "<<std::endl;



  std::cout<<" OutputDirectory :    "<<OutputDirectory<<std::endl;
  std::cout<<"         "<<std::endl;

  std::cout<<" OutputFileName :     "<<OutputFileName<<std::endl;
  std::cout<<"         "<<std::endl;


  std::string TreeName ;
  try{ TreeName = gConfigParser -> readStringOption("Input::TreeName"); }
  catch( char const* exceptionString){ TreeName = "WJet" ;
                                       std::cerr<<" TreeName set by default to --> WJet" <<std::endl;
  }

  float BoostedPtCut =  gConfigParser -> readFloatOption("Input::BoostedPtCut");
  float MetCut       =  gConfigParser -> readFloatOption("Input::MetCut");

  std::cout<<" TreeName :          "<<TreeName<<std::endl;
  std::cout<<"         "<<std::endl;

  float dR_treshold ;
  try{ dR_treshold = gConfigParser -> readFloatOption("Option::DRCut"); }
  catch(char const* exceptionString){ dR_treshold = 0.3 ;
                                      std::cerr<<" DR Cut Set To --> 0.3"<<std::endl;
  }

  std::cout<<" dR_treshold :    "<<dR_treshold<<std::endl;
  std::cout<<"         "<<std::endl;

  float pt_treshold ;
  try{ pt_treshold = gConfigParser -> readFloatOption("Option::PtCut"); }
  catch(char const* exceptionString){ pt_treshold = 0.3 ;
                                      std::cerr<<" pt Cut Set To --> 0.3"<<std::endl;
  }

  std::cout<<" pt_treshold :    "<<pt_treshold<<std::endl;
  std::cout<<"         "<<std::endl;


  int NBin ;

  try{ NBin = gConfigParser -> readIntOption("Option::NBin"); }
  catch(char const* exceptionString){ NBin = 0.6 ;
                                      std::cerr<<" NBin Set To --> 50"<<std::endl;
  }

  std::cout<<" NBin :    "<<NBin<<std::endl;
  std::cout<<"         "<<std::endl;


  float DRMax, DRMin ;

  try{ DRMax = gConfigParser -> readFloatOption("Option::DRMax"); }
  catch(char const* exceptionString){ DRMax = 0.6 ;
                                      std::cerr<<" DRMax Set To --> 0.6"<<std::endl;
  }

  try{ DRMin = gConfigParser -> readFloatOption("Option::DRMin"); }
  catch(char const* exceptionString){ DRMin = 0. ;
                                      std::cerr<<" DRMin Set To --> 0."<<std::endl;
  }

  std::cout<<" DRMin :    "<<DRMin<<std::endl;
  std::cout<<" DRMax :    "<<DRMax<<std::endl;
  std::cout<<"            "<<std::endl;

  float DEtaLimit;

  try{ DEtaLimit = gConfigParser -> readFloatOption("Option::DEtaLimit"); }
  catch(char const* exceptionString){ DEtaLimit = 0.3 ;
                                      std::cerr<<" Eta Limit Set To --> 5"<<std::endl;
  }

  std::cout<<" DEtaLimit :    "<<DEtaLimit<<std::endl;
  std::cout<<"         "<<std::endl;

  float DPhiLimit;

  try{ DPhiLimit = gConfigParser -> readFloatOption("Option::DPhiLimit"); }
  catch(char const* exceptionString){ DEtaLimit = 3.14 ;
                                      std::cerr<<" Phi Limit Set To --> 3.14"<<std::endl;
  }

  std::cout<<" DPhiLimit :    "<<DPhiLimit<<std::endl;
  std::cout<<"         "<<std::endl;

  float DPtMax, DPtMin ;

  try{ DPtMax = gConfigParser -> readFloatOption("Option::DPtMax"); }
  catch(char const* exceptionString){ DPtMax = 2 ;
                                      std::cerr<<" DPt max Set To --> 2"<<std::endl;
  }

  try{ DPtMin = gConfigParser -> readFloatOption("Option::DPtMin"); }
  catch(char const* exceptionString){ DPtMin = 0 ;
                                      std::cerr<<" DPt max Set To --> 0"<<std::endl;
  }


  std::cout<<" DPtMin :    "<<DPtMin<<std::endl;
  std::cout<<" DPtMax :    "<<DPtMax<<std::endl;
  std::cout<<"         "<<std::endl;


  // Open Input File and Get Input Tree

  TFile * inputFile  = new TFile((InputDirectory+"/"+InputFileName).c_str(),"READ");

  if(inputFile==0){ std::cerr<<" No Input File found --> exit "<<std::endl; return -1 ;}
  
  TTree * fTree      = (TTree*) inputFile -> Get(TreeName.c_str());

  if(fTree==0){ std::cerr<<" No Input Tree found --> exit "<<std::endl;  return -1; }

  // Create Output File

  std::string command = " if [ ! -e "+OutputDirectory+"/ ] ; then mkdir "+OutputDirectory+" ; fi " ;
  system(command.c_str());

  TFile * OutputFile  = new TFile((OutputDirectory+"/"+OutputFileName).c_str(),"RECREATE");
  OutputFile->cd();


  // Histogramm in the output


  TH1F * histo_dR1_maxpt   = new TH1F ("histo_dR1_maxpt", "histo_dR1_maxpt", NBin, DRMin, DRMax);
  histo_dR1_maxpt -> GetXaxis() -> SetTitle ("#DeltaR (j1,q1)");
  histo_dR1_maxpt->Sumw2();
 
  TH1F * histo_dR2_maxpt   = new TH1F ("histo_dR2_maxpt", "histo_dR2_maxpt",  NBin, DRMin, DRMax);
  histo_dR2_maxpt -> GetXaxis() -> SetTitle ("#DeltaR (j2, q2)");
  histo_dR2_maxpt->Sumw2();

  TH1F * histo_deta1_maxpt = new TH1F ("histo_deta1_maxpt", "histo_deta1_maxpt", NBin, -DEtaLimit, DEtaLimit);
  histo_deta1_maxpt -> GetXaxis() -> SetTitle ("#Delta#eta (j1, q1)");
  histo_deta1_maxpt->Sumw2();

  TH1F * histo_deta2_maxpt = new TH1F ("histo_deta2_maxpt", "histo_deta2_maxpt", NBin, -DEtaLimit, DEtaLimit);
  histo_deta2_maxpt -> GetXaxis() -> SetTitle ("#Delta#eta (j2, q2)");
  histo_deta2_maxpt->Sumw2();

  TH1F * histo_dphi1_maxpt = new TH1F ("histo_dphi1_maxpt", "histo_deta1_maxpt", NBin, -DPhiLimit, DPhiLimit);
  histo_dphi1_maxpt -> GetXaxis() -> SetTitle ("#Delta#phi (j1,q1)");
  histo_dphi1_maxpt->Sumw2();
  TH1F * histo_dphi2_maxpt = new TH1F ("histo_dphi2_maxpt", "histo_deta2_maxpt",  NBin, -DPhiLimit, DPhiLimit);
  histo_dphi2_maxpt -> GetXaxis() -> SetTitle ("#Delta#phi (j1,q1)");
  histo_dphi2_maxpt->Sumw2();

  TH1F * histo_dpt1_maxpt  = new TH1F ("histo_dpt1_maxpt", "histo_dpt1_maxpt", NBin, DPtMin, DPtMax );
  histo_dpt1_maxpt -> GetXaxis() -> SetTitle ("P_{t} ratio (j1,q1)");
  histo_dpt1_maxpt->Sumw2();

  TH1F * histo_dpt2_maxpt  = new TH1F ("histo_dpt2_maxpt", "histo_dpt2_maxpt", NBin, DPtMin, DPtMax );
  histo_dpt2_maxpt -> GetXaxis() -> SetTitle ("P_{t} ratio (j2,q2)");
  histo_dpt2_maxpt->Sumw2();

  TH2F* histo_dR1_dpt1_maxpt = new TH2F("histo_dR1_dpt1_maxpt","histo_dR1_dpt1_maxpt",NBin,DRMin,DRMax,NBin,DPtMin,DPtMax);
  histo_dR1_dpt1_maxpt -> GetXaxis() -> SetTitle ("#DeltaR (j1,q1)");
  histo_dR1_dpt1_maxpt -> GetYaxis() -> SetTitle ("P_{t} Ratio (j1,q1)");
  histo_dR1_dpt1_maxpt->Sumw2();

  TH2F* histo_dR2_dpt2_maxpt = new TH2F("histo_dR2_dpt2_maxpt","histo_dR2_dpt2_maxpt",NBin,DRMin,DRMax,NBin,DPtMin,DPtMax);
  histo_dR2_dpt2_maxpt -> GetXaxis() -> SetTitle ("#DeltaR (j2,q2)");
  histo_dR2_dpt2_maxpt -> GetYaxis() -> SetTitle ("P_{t} Ratio (j2,q2)");
  histo_dR2_dpt2_maxpt->Sumw2();


  TH1F * histo_dR1_maxDeta   = new TH1F ("histo_dR1_maxDeta", "histo_dR1_maxDeta",  NBin, DRMin, DRMax);
  histo_dR1_maxDeta -> GetXaxis() -> SetTitle ("#DeltaR (j1,q1)");
  histo_dR1_maxDeta->Sumw2();
 
  TH1F * histo_dR2_maxDeta   = new TH1F ("histo_dR2_maxDeta", "histo_dR2_maxDeta",  NBin, DRMin, DRMax);
  histo_dR2_maxDeta -> GetXaxis() -> SetTitle ("#DeltaR (j2, q2)");
  histo_dR2_maxDeta->Sumw2();

  TH1F * histo_deta1_maxDeta = new TH1F ("histo_deta1_maxDeta", "histo_deta1_maxDeta", NBin, -DEtaLimit, DEtaLimit);
  histo_deta1_maxDeta -> GetXaxis() -> SetTitle ("#Delta#eta (j1, q1)");
  histo_deta1_maxDeta->Sumw2();

  TH1F * histo_deta2_maxDeta = new TH1F ("histo_deta2_maxDeta", "histo_deta2_maxDeta", NBin, -DEtaLimit, DEtaLimit);
  histo_deta2_maxDeta -> GetXaxis() -> SetTitle ("#Delta#eta (j2, q2)");
  histo_deta2_maxDeta->Sumw2();

  TH1F * histo_dphi1_maxDeta = new TH1F ("histo_dphi1_maxDeta", "histo_deta1_maxDeta",  NBin, -DPhiLimit, DPhiLimit);
  histo_dphi1_maxDeta -> GetXaxis() -> SetTitle ("#Delta#phi (j1,q1)");
  histo_dphi1_maxDeta->Sumw2();
  TH1F * histo_dphi2_maxDeta = new TH1F ("histo_dphi2_maxDeta", "histo_deta2_maxDeta",  NBin, -DPhiLimit, DPhiLimit);
  histo_dphi2_maxDeta -> GetXaxis() -> SetTitle ("#Delta#phi (j1,q1)");
  histo_dphi2_maxDeta->Sumw2();

  TH1F * histo_dpt1_maxDeta  = new TH1F ("histo_dpt1_maxDeta", "histo_dpt1_maxDeta", NBin, DPtMin, DPtMax );
  histo_dpt1_maxDeta -> GetXaxis() -> SetTitle ("P_{t} ratio (j1,q1)");
  histo_dpt1_maxDeta->Sumw2();

  TH1F * histo_dpt2_maxDeta  = new TH1F ("histo_dpt2_maxDeta", "histo_dpt2_maxDeta", NBin, DPtMin, DPtMax );
  histo_dpt2_maxDeta -> GetXaxis() -> SetTitle ("P_{t} ratio (j2,q2)");
  histo_dpt2_maxDeta->Sumw2();

  TH2F* histo_dR1_dpt1_maxDeta = new TH2F("histo_dR1_dpt1_maxDeta","histo_dR1_dpt1_maxDeta",NBin,DRMin,DRMax,NBin,DPtMin,DPtMax);
  histo_dR1_dpt1_maxDeta -> GetXaxis() -> SetTitle ("#DeltaR (j1,q1)");
  histo_dR1_dpt1_maxDeta -> GetYaxis() -> SetTitle ("P_{t} Ratio (j1,q1)");
  histo_dR1_dpt1_maxDeta->Sumw2();

  TH2F* histo_dR2_dpt2_maxDeta = new TH2F("histo_dR2_dpt2_maxDeta","histo_dR2_dpt2_maxDeta",NBin,DRMin,DRMax,NBin,DPtMin,DPtMax);
  histo_dR2_dpt2_maxDeta -> GetXaxis() -> SetTitle ("#DeltaR (j2,q2)");
  histo_dR2_dpt2_maxDeta -> GetYaxis() -> SetTitle ("P_{t} Ratio (j2,q2)");
  histo_dR2_dpt2_maxDeta->Sumw2();



  TH1F * histo_dR1_maxMjj   = new TH1F ("histo_dR1_maxMjj", "histo_dR1_maxMjj",  NBin, DRMin, DRMax);
  histo_dR1_maxMjj -> GetXaxis() -> SetTitle ("#DeltaR (j1,q1)");
  histo_dR1_maxMjj->Sumw2();
 
  TH1F * histo_dR2_maxMjj   = new TH1F ("histo_dR2_maxMjj", "histo_dR2_maxMjj",  NBin, DRMin, DRMax);
  histo_dR2_maxMjj -> GetXaxis() -> SetTitle ("#DeltaR (j2, q2)");
  histo_dR2_maxMjj->Sumw2();

  TH1F * histo_deta1_maxMjj = new TH1F ("histo_deta1_maxMjj", "histo_deta1_maxMjj", NBin, -DEtaLimit, DEtaLimit);
  histo_deta1_maxMjj -> GetXaxis() -> SetTitle ("#Delta#eta (j1, q1)");
  histo_deta1_maxMjj->Sumw2();

  TH1F * histo_deta2_maxMjj = new TH1F ("histo_deta2_maxMjj", "histo_deta2_maxMjj", NBin, -DEtaLimit, DEtaLimit);
  histo_deta2_maxMjj -> GetXaxis() -> SetTitle ("#Delta#eta (j2, q2)");
  histo_deta2_maxMjj->Sumw2();

  TH1F * histo_dphi1_maxMjj = new TH1F ("histo_dphi1_maxMjj", "histo_deta1_maxMjj",  NBin, -DPhiLimit, DPhiLimit);
  histo_dphi1_maxMjj -> GetXaxis() -> SetTitle ("#Delta#phi (j1,q1)");
  histo_dphi1_maxMjj->Sumw2();
  TH1F * histo_dphi2_maxMjj = new TH1F ("histo_dphi2_maxMjj", "histo_deta2_maxMjj",  NBin, -DPhiLimit, DPhiLimit);
  histo_dphi2_maxMjj -> GetXaxis() -> SetTitle ("#Delta#phi (j1,q1)");
  histo_dphi2_maxMjj->Sumw2();

  TH1F * histo_dpt1_maxMjj  = new TH1F ("histo_dpt1_maxMjj", "histo_dpt1_maxMjj", NBin, DPtMin, DPtMax );
  histo_dpt1_maxMjj -> GetXaxis() -> SetTitle ("P_{t} ratio (j1,q1)");
  histo_dpt1_maxMjj->Sumw2();

  TH1F * histo_dpt2_maxMjj  = new TH1F ("histo_dpt2_maxMjj", "histo_dpt2_maxMjj", NBin, DPtMin, DPtMax );
  histo_dpt2_maxMjj -> GetXaxis() -> SetTitle ("P_{t} ratio (j2,q2)");
  histo_dpt2_maxMjj->Sumw2();

  TH2F* histo_dR1_dpt1_maxMjj = new TH2F("histo_dR1_dpt1_maxMjj","histo_dR1_dpt1_maxMjj",NBin,DRMin,DRMax,NBin,DPtMin,DPtMax);
  histo_dR1_dpt1_maxMjj -> GetXaxis() -> SetTitle ("#DeltaR (j1,q1)");
  histo_dR1_dpt1_maxMjj -> GetYaxis() -> SetTitle ("P_{t} Ratio (j1,q1)");
  histo_dR1_dpt1_maxMjj->Sumw2();

  TH2F* histo_dR2_dpt2_maxMjj = new TH2F("histo_dR2_dpt2_maxMjj","histo_dR2_dpt2_maxMjj",NBin,DRMin,DRMax,NBin,DPtMin,DPtMax);
  histo_dR2_dpt2_maxMjj -> GetXaxis() -> SetTitle ("#DeltaR (j2,q2)");
  histo_dR2_dpt2_maxMjj -> GetYaxis() -> SetTitle ("P_{t} Ratio (j2,q2)");
  histo_dR2_dpt2_maxMjj->Sumw2();

  //Declaration of variables

  float  vbf_maxpt_j1_e, vbf_maxpt_j1_pt, vbf_maxpt_j1_eta, vbf_maxpt_j1_phi;
  bool   vbf_maxpt_j1_isPileUpLoose, vbf_maxpt_j1_isPileUpMedium, vbf_maxpt_j1_isPileUpTight;
  float  vbf_maxpt_j2_e, vbf_maxpt_j2_pt, vbf_maxpt_j2_eta, vbf_maxpt_j2_phi;
  bool   vbf_maxpt_j2_isPileUpLoose, vbf_maxpt_j2_isPileUpMedium, vbf_maxpt_j2_isPileUpTight;
  
  float vbf_maxDeta_j1_e, vbf_maxDeta_j1_pt, vbf_maxDeta_j1_eta, vbf_maxDeta_j1_phi;
  bool  vbf_maxDeta_j1_isPileUpLoose, vbf_maxDeta_j1_isPileUpMedium, vbf_maxDeta_j1_isPileUpTight;
  float vbf_maxDeta_j2_e, vbf_maxDeta_j2_pt, vbf_maxDeta_j2_eta, vbf_maxDeta_j2_phi;
  bool  vbf_maxDeta_j2_isPileUpLoose, vbf_maxDeta_j2_isPileUpMedium, vbf_maxDeta_j2_isPileUpTight;

  float vbf_maxMjj_j1_e, vbf_maxMjj_j1_pt, vbf_maxMjj_j1_eta, vbf_maxMjj_j1_phi;
  bool  vbf_maxMjj_j1_isPileUpLoose, vbf_maxMjj_j1_isPileUpMedium, vbf_maxMjj_j1_isPileUpTight;
  float vbf_maxMjj_j2_e, vbf_maxMjj_j2_pt, vbf_maxMjj_j2_eta, vbf_maxMjj_j2_phi;
  bool  vbf_maxMjj_j2_isPileUpLoose, vbf_maxMjj_j2_isPileUpMedium, vbf_maxMjj_j2_isPileUpTight; 

  float W_TagQuark_E[2], W_TagQuark_pt[2], W_TagQuark_eta[2], W_TagQuark_phi[2];

  int numPFCorJets, numPFCorVBFTagJets;

  float GroomedJet_CA8_pt[8], W_pt, event_met_pfmet;
  
  //Declaration  of branches

  fTree -> SetBranchAddress ("GroomedJet_CA8_pt", &GroomedJet_CA8_pt);
  fTree -> SetBranchAddress ("W_pt", &W_pt);
  fTree -> SetBranchAddress ("event_met_pfmet", &event_met_pfmet);
  
  fTree -> SetBranchAddress ("vbf_maxpt_j1_e", &vbf_maxpt_j1_e);
  fTree -> SetBranchAddress ("vbf_maxpt_j1_pt", &vbf_maxpt_j1_pt);
  fTree -> SetBranchAddress ("vbf_maxpt_j1_eta", &vbf_maxpt_j1_eta);
  fTree -> SetBranchAddress ("vbf_maxpt_j1_phi", &vbf_maxpt_j1_phi);
  fTree -> SetBranchAddress ("vbf_maxpt_j1_isPileUpLoose", &vbf_maxpt_j1_isPileUpLoose);
  fTree -> SetBranchAddress ("vbf_maxpt_j1_isPileUpMedium", &vbf_maxpt_j1_isPileUpMedium);
  fTree -> SetBranchAddress ("vbf_maxpt_j1_isPileUpTight", &vbf_maxpt_j1_isPileUpTight);

  fTree -> SetBranchAddress ("vbf_maxpt_j2_e", &vbf_maxpt_j2_e);
  fTree -> SetBranchAddress ("vbf_maxpt_j2_pt", &vbf_maxpt_j2_pt);
  fTree -> SetBranchAddress ("vbf_maxpt_j2_eta", &vbf_maxpt_j2_eta);
  fTree -> SetBranchAddress ("vbf_maxpt_j2_phi", &vbf_maxpt_j2_phi);
  fTree -> SetBranchAddress ("vbf_maxpt_j2_isPileUpLoose", &vbf_maxpt_j2_isPileUpLoose);
  fTree -> SetBranchAddress ("vbf_maxpt_j2_isPileUpMedium", &vbf_maxpt_j2_isPileUpMedium);
  fTree -> SetBranchAddress ("vbf_maxpt_j2_isPileUpTight", &vbf_maxpt_j2_isPileUpTight);

  fTree -> SetBranchAddress ("vbf_maxDeta_j1_e", &vbf_maxDeta_j1_e);
  fTree -> SetBranchAddress ("vbf_maxDeta_j1_pt", &vbf_maxDeta_j1_pt);
  fTree -> SetBranchAddress ("vbf_maxDeta_j1_eta", &vbf_maxDeta_j1_eta);
  fTree -> SetBranchAddress ("vbf_maxDeta_j1_phi", &vbf_maxDeta_j1_phi);
  fTree -> SetBranchAddress ("vbf_maxDeta_j1_isPileUpLoose", &vbf_maxDeta_j1_isPileUpLoose);
  fTree -> SetBranchAddress ("vbf_maxDeta_j1_isPileUpMedium", &vbf_maxDeta_j1_isPileUpMedium);
  fTree -> SetBranchAddress ("vbf_maxDeta_j1_isPileUpTight", &vbf_maxDeta_j1_isPileUpTight);

  fTree -> SetBranchAddress ("vbf_maxDeta_j2_e", &vbf_maxDeta_j2_e);
  fTree -> SetBranchAddress ("vbf_maxDeta_j2_pt", &vbf_maxDeta_j2_pt);
  fTree -> SetBranchAddress ("vbf_maxDeta_j2_eta", &vbf_maxDeta_j2_eta);
  fTree -> SetBranchAddress ("vbf_maxDeta_j2_phi", &vbf_maxDeta_j2_phi);
  fTree -> SetBranchAddress ("vbf_maxDeta_j2_isPileUpLoose", &vbf_maxDeta_j2_isPileUpLoose);
  fTree -> SetBranchAddress ("vbf_maxDeta_j2_isPileUpMedium", &vbf_maxDeta_j2_isPileUpMedium);
  fTree -> SetBranchAddress ("vbf_maxDeta_j2_isPileUpTight", &vbf_maxDeta_j2_isPileUpTight);

  fTree -> SetBranchAddress ("vbf_maxMjj_j1_e", &vbf_maxMjj_j1_e);
  fTree -> SetBranchAddress ("vbf_maxMjj_j1_pt", &vbf_maxMjj_j1_pt);
  fTree -> SetBranchAddress ("vbf_maxMjj_j1_eta", &vbf_maxMjj_j1_eta);
  fTree -> SetBranchAddress ("vbf_maxMjj_j1_phi", &vbf_maxMjj_j1_phi);
  fTree -> SetBranchAddress ("vbf_maxMjj_j1_isPileUpLoose", &vbf_maxMjj_j1_isPileUpLoose);
  fTree -> SetBranchAddress ("vbf_maxMjj_j1_isPileUpMedium", &vbf_maxMjj_j1_isPileUpMedium);
  fTree -> SetBranchAddress ("vbf_maxMjj_j1_isPileUpTight", &vbf_maxMjj_j1_isPileUpTight);

  fTree -> SetBranchAddress ("vbf_maxMjj_j2_e", &vbf_maxMjj_j2_e);
  fTree -> SetBranchAddress ("vbf_maxMjj_j2_pt", &vbf_maxMjj_j2_pt);
  fTree -> SetBranchAddress ("vbf_maxMjj_j2_eta", &vbf_maxMjj_j2_eta);
  fTree -> SetBranchAddress ("vbf_maxMjj_j2_phi", &vbf_maxMjj_j2_phi);
  fTree -> SetBranchAddress ("vbf_maxMjj_j2_isPileUpLoose", &vbf_maxMjj_j2_isPileUpLoose);
  fTree -> SetBranchAddress ("vbf_maxMjj_j2_isPileUpMedium", &vbf_maxMjj_j2_isPileUpMedium);
  fTree -> SetBranchAddress ("vbf_maxMjj_j2_isPileUpTight", &vbf_maxMjj_j2_isPileUpTight);

  fTree -> SetBranchAddress ("W_TagQuark_E[2]", W_TagQuark_E);
  fTree -> SetBranchAddress ("W_TagQuark_pt[2]", W_TagQuark_pt);
  fTree -> SetBranchAddress ("W_TagQuark_eta[2]", W_TagQuark_eta);
  fTree -> SetBranchAddress ("W_TagQuark_phi[2]", W_TagQuark_phi);

  fTree -> SetBranchAddress ("numPFCorJets", &numPFCorJets);
  fTree -> SetBranchAddress ("numPFCorVBFTagJets", &numPFCorVBFTagJets);
  
  //declaration of variables

  int nEntry = fTree->GetEntries();
  
  float dR11_maxpt   = 0., dR12_maxpt   = 0.;
  float dR11_maxDeta = 0., dR12_maxDeta = 0.;
  float dR11_maxMjj  = 0., dR12_maxMjj  = 0.; 
  float dR21_maxpt   = 0., dR22_maxpt   = 0.;
  float dR21_maxDeta = 0., dR22_maxDeta = 0.;
  float dR21_maxMjj  = 0., dR22_maxMjj  = 0.;

  int cont_maxpt0 =0, cont_maxDeta0 =0, cont_maxMjj0 =0;
  int cont_maxpt1 =0, cont_maxDeta1 =0, cont_maxMjj1 =0;
  int cont_maxpt2 =0, cont_maxDeta2 =0, cont_maxMjj2 =0;
 
  float dR1_maxpt =0., dR1_maxDeta =0., dR1_maxMjj =0.;
  float dR2_maxpt =0., dR2_maxDeta =0., dR2_maxMjj =0.;
  
  float ptratio1_maxpt=0., ptratio1_maxDeta=0., ptratio1_maxMjj=0.;
  float ptratio2_maxpt=0., ptratio2_maxDeta=0., ptratio2_maxMjj=0.;

  int   cont_dRandpt_maxpt = 0, cont_dRandpt_maxDeta=0, cont_dRandpt_maxMjj=0;
  int   cont_maxpt_maxDeta = 0, cont_maxpt_maxMjj=0, cont_maxMjj_maxDeta =0;

  std::vector<int> cont_nj (12,0) ;

  std::vector<int> cont_nj_match_maxpt (12,0); 
  std::vector<int> cont_nj_match_maxDeta (12,0); 
  std::vector<int> cont_nj_match_maxMjj (12,0);  
  std::vector<int> cont_nj_match_dR_maxpt (12,0);
  std::vector<int> cont_nj_match_dR_maxDeta (12,0);
  std::vector<int> cont_nj_match_dR_maxMjj (12,0) ;

  int cont_tight_maxpt=0, cont_tight_maxDeta=0, cont_tight_maxMjj=0;
  int cont_tight_match_maxpt=0, cont_tight_match_maxDeta=0, cont_tight_match_maxMjj=0;
  int cont_tight_match_dR_maxpt=0, cont_tight_match_dR_maxDeta=0, cont_tight_match_dR_maxMjj=0;
  int cont_medium_maxpt=0, cont_medium_maxDeta=0, cont_medium_maxMjj=0;
  int cont_medium_match_maxpt=0, cont_medium_match_maxDeta=0, cont_medium_match_maxMjj=0;
  int cont_medium_match_dR_maxpt=0, cont_medium_match_dR_maxDeta=0, cont_medium_match_dR_maxMjj=0;
  

  int nPassingEvents =0;
  
  //start reading tree

  std::cout<< " Input Numeber of events in the tree : "<<nEntry<<std::endl;
  std::cout<<"   "<<std::endl;
 
  for(int iEntry =0; iEntry<nEntry ; iEntry++){


    dR11_maxpt   = 0., dR12_maxpt   = 0.; dR11_maxDeta = 0., dR12_maxDeta = 0.; dR11_maxMjj  = 0., dR12_maxMjj  = 0.; 
    dR21_maxpt   = 0., dR22_maxpt   = 0.; dR21_maxDeta = 0., dR22_maxDeta = 0.; dR21_maxMjj  = 0., dR22_maxMjj  = 0.;

    dR1_maxpt =0., dR1_maxDeta =0., dR1_maxMjj =0.;
    dR2_maxpt =0., dR2_maxDeta =0., dR2_maxMjj =0.;
  
    ptratio1_maxpt=0., ptratio1_maxDeta=0., ptratio1_maxMjj=0.; ptratio2_maxpt=0., ptratio2_maxDeta=0., ptratio2_maxMjj=0.;

 
    if( iEntry%1000 == 0 ) std::cout << "reading saved entry " << iEntry << "\r" << std::endl;
 
    fTree->GetEntry(iEntry);
    
    if(GroomedJet_CA8_pt[0]<BoostedPtCut || W_pt < BoostedPtCut || event_met_pfmet < MetCut ) continue;
    if(vbf_maxpt_j1_pt <=0 || vbf_maxpt_j2_pt <=0 || W_TagQuark_pt[0] <=0 || W_TagQuark_pt[1] <=0 ) continue ;

    //calculation of dR for max pt
        
    dR11_maxpt = deltaR(W_TagQuark_phi[0],vbf_maxpt_j1_phi,W_TagQuark_eta[0],vbf_maxpt_j1_eta); 
    dR12_maxpt = deltaR(W_TagQuark_phi[0],vbf_maxpt_j2_phi,W_TagQuark_eta[0],vbf_maxpt_j2_eta); 
    dR21_maxpt = deltaR(W_TagQuark_phi[1],vbf_maxpt_j1_phi,W_TagQuark_eta[1],vbf_maxpt_j1_eta); 
    dR22_maxpt = deltaR(W_TagQuark_phi[1],vbf_maxpt_j2_phi,W_TagQuark_eta[1],vbf_maxpt_j2_eta); 

    if(dR11_maxpt < dR12_maxpt && dR22_maxpt<dR21_maxpt){
   
      dR1_maxpt = dR11_maxpt ;
      dR2_maxpt = dR22_maxpt ;

      ptratio1_maxpt = vbf_maxpt_j1_pt/W_TagQuark_pt[0] ;
      ptratio2_maxpt = vbf_maxpt_j2_pt/W_TagQuark_pt[1] ;

      histo_dR1_maxpt->Fill(dR11_maxpt);
      histo_dR2_maxpt->Fill(dR22_maxpt);
      histo_deta1_maxpt->Fill(W_TagQuark_eta[0]-vbf_maxpt_j1_eta);
      histo_deta2_maxpt->Fill(W_TagQuark_eta[1]-vbf_maxpt_j2_eta);
      histo_dphi1_maxpt->Fill(W_TagQuark_eta[0]-vbf_maxpt_j1_eta);
      histo_dphi2_maxpt->Fill(W_TagQuark_eta[1]-vbf_maxpt_j2_eta);
      histo_dpt1_maxpt->Fill(vbf_maxpt_j1_pt/W_TagQuark_pt[0]);
      histo_dpt2_maxpt->Fill(vbf_maxpt_j2_pt/W_TagQuark_pt[1]);

      histo_dR1_dpt1_maxpt->Fill(dR11_maxpt,vbf_maxpt_j1_pt/W_TagQuark_pt[0]);
      histo_dR2_dpt2_maxpt->Fill(dR22_maxpt,vbf_maxpt_j2_pt/W_TagQuark_pt[1]);
    }
    else if (dR12_maxpt <= dR12_maxpt && dR21_maxpt<=dR22_maxpt) {

      dR1_maxpt = dR12_maxpt ;
      dR2_maxpt = dR21_maxpt ;
 
      ptratio1_maxpt = vbf_maxpt_j2_pt/W_TagQuark_pt[0] ;
      ptratio2_maxpt = vbf_maxpt_j1_pt/W_TagQuark_pt[1] ;
 
      histo_dR1_maxpt->Fill(dR12_maxpt);
      histo_dR2_maxpt->Fill(dR21_maxpt);
      histo_deta1_maxpt->Fill(W_TagQuark_eta[0]-vbf_maxpt_j2_eta);
      histo_deta2_maxpt->Fill(W_TagQuark_eta[1]-vbf_maxpt_j1_eta);
      histo_dphi1_maxpt->Fill(W_TagQuark_eta[0]-vbf_maxpt_j2_eta);
      histo_dphi2_maxpt->Fill(W_TagQuark_eta[1]-vbf_maxpt_j1_eta);
      histo_dpt1_maxpt->Fill(vbf_maxpt_j2_pt/W_TagQuark_pt[0]);
      histo_dpt2_maxpt->Fill(vbf_maxpt_j1_pt/W_TagQuark_pt[1]);

      histo_dR1_dpt1_maxpt->Fill(dR12_maxpt,vbf_maxpt_j2_pt/W_TagQuark_pt[0]);
      histo_dR2_dpt2_maxpt->Fill(dR21_maxpt,vbf_maxpt_j1_pt/W_TagQuark_pt[1]);
 
    }

    else continue ;

    nPassingEvents++;

   //same thing whit pair of jets with maxDeta

    dR11_maxDeta = deltaR(W_TagQuark_phi[0],vbf_maxDeta_j1_phi,W_TagQuark_eta[0],vbf_maxDeta_j1_eta);
    dR12_maxDeta = deltaR(W_TagQuark_phi[0],vbf_maxDeta_j2_phi,W_TagQuark_eta[0],vbf_maxDeta_j2_eta);
    dR21_maxDeta = deltaR(W_TagQuark_phi[1],vbf_maxDeta_j1_phi,W_TagQuark_eta[1],vbf_maxDeta_j1_eta);
    dR22_maxDeta = deltaR(W_TagQuark_phi[1],vbf_maxDeta_j2_phi,W_TagQuark_eta[1],vbf_maxDeta_j2_eta);

    if(dR11_maxDeta < dR12_maxDeta){
 
      dR1_maxDeta = dR11_maxDeta ;
      dR2_maxDeta = dR22_maxDeta ;

      ptratio1_maxDeta = vbf_maxDeta_j1_pt/W_TagQuark_pt[0] ;
      ptratio2_maxDeta = vbf_maxDeta_j2_pt/W_TagQuark_pt[1] ;
 
      histo_dR1_maxDeta->Fill(dR11_maxDeta);
      histo_dR2_maxDeta->Fill(dR22_maxDeta);
      histo_deta1_maxDeta->Fill(W_TagQuark_eta[0]-vbf_maxDeta_j1_eta);
      histo_deta2_maxDeta->Fill(W_TagQuark_eta[1]-vbf_maxDeta_j2_eta);
      histo_dphi1_maxDeta->Fill(W_TagQuark_eta[0]-vbf_maxDeta_j1_eta);
      histo_dphi2_maxDeta->Fill(W_TagQuark_eta[1]-vbf_maxDeta_j2_eta);
      histo_dpt1_maxDeta->Fill(vbf_maxDeta_j1_pt/W_TagQuark_pt[0]);
      histo_dpt2_maxDeta->Fill(vbf_maxDeta_j2_pt/W_TagQuark_pt[1]);

      histo_dR1_dpt1_maxDeta->Fill(dR11_maxDeta,vbf_maxDeta_j1_pt/W_TagQuark_pt[0]);
      histo_dR2_dpt2_maxDeta->Fill(dR22_maxDeta,vbf_maxDeta_j2_pt/W_TagQuark_pt[1]);
    }
    else{
 
      dR1_maxDeta = dR12_maxDeta ;
      dR2_maxDeta = dR21_maxDeta ;
 
      ptratio1_maxDeta = vbf_maxDeta_j2_pt/W_TagQuark_pt[0] ;
      ptratio2_maxDeta = vbf_maxDeta_j1_pt/W_TagQuark_pt[1] ;
 
      histo_dR1_maxDeta->Fill(dR12_maxDeta);
      histo_dR2_maxDeta->Fill(dR21_maxDeta);
      histo_deta1_maxDeta->Fill(W_TagQuark_eta[0]-vbf_maxDeta_j2_eta);
      histo_deta2_maxDeta->Fill(W_TagQuark_eta[1]-vbf_maxDeta_j1_eta);
      histo_dphi1_maxDeta->Fill(W_TagQuark_eta[0]-vbf_maxDeta_j2_eta);
      histo_dphi2_maxDeta->Fill(W_TagQuark_eta[1]-vbf_maxDeta_j1_eta);
      histo_dpt1_maxDeta->Fill(vbf_maxDeta_j2_pt/W_TagQuark_pt[0]);
      histo_dpt2_maxDeta->Fill(vbf_maxDeta_j1_pt/W_TagQuark_pt[1]);

      histo_dR1_dpt1_maxDeta->Fill(dR12_maxDeta,vbf_maxDeta_j2_pt/W_TagQuark_pt[0]);
      histo_dR2_dpt2_maxDeta->Fill(dR21_maxDeta,vbf_maxDeta_j1_pt/W_TagQuark_pt[1]);
 
    }

   //same thing whit pair of jets with maxMjj

    dR11_maxMjj = deltaR(W_TagQuark_phi[0],vbf_maxMjj_j1_phi,W_TagQuark_eta[0],vbf_maxMjj_j1_eta);
    dR12_maxMjj = deltaR(W_TagQuark_phi[0],vbf_maxMjj_j2_phi,W_TagQuark_eta[0],vbf_maxMjj_j2_eta);
    dR21_maxMjj = deltaR(W_TagQuark_phi[1],vbf_maxMjj_j1_phi,W_TagQuark_eta[1],vbf_maxMjj_j1_eta);
    dR22_maxMjj = deltaR(W_TagQuark_phi[1],vbf_maxMjj_j2_phi,W_TagQuark_eta[1],vbf_maxMjj_j2_eta);

    if(dR11_maxMjj < dR12_maxMjj){
 
      dR1_maxMjj = dR11_maxMjj ;
      dR2_maxMjj = dR22_maxMjj ;
 
      ptratio1_maxMjj = vbf_maxMjj_j1_pt/W_TagQuark_pt[0] ;
      ptratio2_maxMjj = vbf_maxMjj_j2_pt/W_TagQuark_pt[1] ;
 
      histo_dR1_maxMjj->Fill(dR11_maxMjj);
      histo_dR2_maxMjj->Fill(dR22_maxMjj);
      histo_deta1_maxMjj->Fill(W_TagQuark_eta[0]-vbf_maxMjj_j1_eta);
      histo_deta2_maxMjj->Fill(W_TagQuark_eta[1]-vbf_maxMjj_j2_eta);
      histo_dphi1_maxMjj->Fill(W_TagQuark_eta[0]-vbf_maxMjj_j1_eta);
      histo_dphi2_maxMjj->Fill(W_TagQuark_eta[1]-vbf_maxMjj_j2_eta);
      histo_dpt1_maxMjj->Fill(vbf_maxMjj_j1_pt/W_TagQuark_pt[0]);
      histo_dpt2_maxMjj->Fill(vbf_maxMjj_j2_pt/W_TagQuark_pt[1]);

      histo_dR1_dpt1_maxMjj->Fill(dR11_maxMjj,vbf_maxMjj_j1_pt/W_TagQuark_pt[0]);
      histo_dR2_dpt2_maxMjj->Fill(dR22_maxMjj,vbf_maxMjj_j2_pt/W_TagQuark_pt[1]);
    }
    else{
 
      dR1_maxMjj = dR12_maxMjj ;
      dR2_maxMjj = dR21_maxMjj ;
  
      ptratio1_maxMjj = vbf_maxMjj_j2_pt/W_TagQuark_pt[0] ;
      ptratio2_maxMjj = vbf_maxMjj_j1_pt/W_TagQuark_pt[1] ;
 
      histo_dR1_maxMjj->Fill(dR12_maxMjj);
      histo_dR2_maxMjj->Fill(dR21_maxMjj);
      histo_deta1_maxMjj->Fill(W_TagQuark_eta[0]-vbf_maxMjj_j2_eta);
      histo_deta2_maxMjj->Fill(W_TagQuark_eta[1]-vbf_maxMjj_j1_eta);
      histo_dphi1_maxMjj->Fill(W_TagQuark_eta[0]-vbf_maxMjj_j2_eta);
      histo_dphi2_maxMjj->Fill(W_TagQuark_eta[1]-vbf_maxMjj_j1_eta);
      histo_dpt1_maxMjj->Fill(vbf_maxMjj_j2_pt/W_TagQuark_pt[0]);
      histo_dpt2_maxMjj->Fill(vbf_maxMjj_j1_pt/W_TagQuark_pt[1]);

      histo_dR1_dpt1_maxMjj->Fill(dR12_maxMjj,vbf_maxMjj_j2_pt/W_TagQuark_pt[0]);
      histo_dR2_dpt2_maxMjj->Fill(dR21_maxMjj,vbf_maxMjj_j1_pt/W_TagQuark_pt[1]);
 
    }

    // Count events which satisfy dR matching condition

    if ( ((dR1_maxpt   <= dR_treshold) && (dR2_maxpt   > dR_treshold)) || ((dR1_maxpt   >  dR_treshold) && (dR2_maxpt   <= dR_treshold)) )  cont_maxpt1++;
    if ( ((dR1_maxDeta <= dR_treshold) && (dR2_maxDeta > dR_treshold)) || ((dR1_maxDeta >  dR_treshold) && (dR2_maxDeta <= dR_treshold)) )  cont_maxDeta1++;
    if ( ((dR1_maxMjj  <= dR_treshold) && (dR2_maxMjj  > dR_treshold)) || ((dR1_maxMjj  >  dR_treshold) && (dR2_maxMjj  <= dR_treshold)) )  cont_maxMjj1++;

    //both jets have matching in dR
    if ( (dR1_maxpt   <= dR_treshold) && (dR2_maxpt   <= dR_treshold) )  cont_maxpt2++;
    if ( (dR1_maxDeta <= dR_treshold) && (dR2_maxDeta <= dR_treshold) )  cont_maxDeta2++;
    if ( (dR1_maxMjj  <= dR_treshold) && (dR2_maxMjj  <= dR_treshold) )  cont_maxMjj2++;
    
    //no jet has matching in dR
    if ( (dR1_maxpt   > dR_treshold) && (dR2_maxpt   > dR_treshold) )   cont_maxpt0++;
    if ( (dR1_maxDeta > dR_treshold) && (dR2_maxDeta > dR_treshold) )   cont_maxDeta0++;
    if ( (dR1_maxMjj  > dR_treshold) && (dR2_maxMjj  > dR_treshold) )   cont_maxMjj0++;

    //both jets have matching in dR and ptratio
    int flag_dRandpt_maxpt=0, flag_dRandpt_maxDeta=0, flag_dRandpt_maxMjj=0;

    if ( (dR1_maxpt <= dR_treshold) && (dR2_maxpt <= dR_treshold) && ( fabs(ptratio1_maxpt-1.) <= pt_treshold) && ( fabs(ptratio2_maxpt-1.) <= pt_treshold) ) {
      cont_dRandpt_maxpt++;
      flag_dRandpt_maxpt++;
    }
    
    if ( (dR1_maxDeta <= dR_treshold) && (dR2_maxDeta <= dR_treshold ) && ( fabs(ptratio1_maxDeta-1.) <= pt_treshold) && ( fabs(ptratio2_maxDeta-1.) <= pt_treshold) ) {
      cont_dRandpt_maxDeta++;
      flag_dRandpt_maxDeta++;
    }

    if ( (dR1_maxMjj <= dR_treshold) && (dR2_maxMjj <= dR_treshold) && ( fabs(ptratio1_maxMjj-1.) <= pt_treshold) && (fabs(ptratio2_maxMjj-1.) <= pt_treshold) ) {
      cont_dRandpt_maxMjj++;
      flag_dRandpt_maxMjj++;
    }

    if ( (flag_dRandpt_maxpt == 1 ) && (flag_dRandpt_maxDeta==1) )  cont_maxpt_maxDeta++;
    if ( (flag_dRandpt_maxpt ==1  ) && (flag_dRandpt_maxMjj==1 ) )  cont_maxpt_maxMjj++;
    if ( (flag_dRandpt_maxMjj ==1 ) && (flag_dRandpt_maxDeta==1) )  cont_maxMjj_maxDeta++;
   
    int numjets = numPFCorJets + numPFCorVBFTagJets;
    
    cont_nj.at(numjets)++;

    if ( (dR1_maxpt   <= dR_treshold ) && (dR2_maxpt  <= dR_treshold) && (fabs(ptratio1_maxpt-1.)   <= pt_treshold) && (fabs(ptratio2_maxpt-1.)   <= pt_treshold) )  
     cont_nj_match_maxpt.at(numjets)++;
 
    if ( (dR1_maxDeta <= dR_treshold) && (dR2_maxDeta <= dR_treshold) && (fabs(ptratio1_maxDeta-1.) <= pt_treshold) && (fabs(ptratio2_maxDeta-1.) <= pt_treshold) ) 
     cont_nj_match_maxDeta.at(numjets)++;

    if ( (dR1_maxMjj  <= dR_treshold) && (dR2_maxMjj  <= dR_treshold) && (fabs(ptratio1_maxMjj-1.)  <= pt_treshold) && (fabs(ptratio2_maxMjj-1.)  <= pt_treshold) ) 
      cont_nj_match_maxMjj.at(numjets)++;

    if ( (dR1_maxpt   <= dR_treshold) && (dR2_maxpt   <= dR_treshold) )  cont_nj_match_dR_maxpt.at(numjets)++;
    if ( (dR1_maxDeta <= dR_treshold) && (dR2_maxDeta <= dR_treshold) )  cont_nj_match_dR_maxDeta.at(numjets)++;
    if ( (dR1_maxMjj  <= dR_treshold) && (dR2_maxMjj  <= dR_treshold) )  cont_nj_match_dR_maxMjj.at(numjets)++;

    //matching with pileUp conditions

    if ( (vbf_maxpt_j1_isPileUpTight==true) && (vbf_maxpt_j2_isPileUpTight ==true) ) {
      cont_tight_maxpt++;

      if ((dR1_maxpt <= dR_treshold) && (dR2_maxpt <= dR_treshold) )   cont_tight_match_dR_maxpt++;
      if ((dR1_maxpt <= dR_treshold) && (dR2_maxpt <= dR_treshold) && ( fabs(ptratio1_maxpt-1.) <= pt_treshold) && ( fabs(ptratio2_maxpt-1.) <= pt_treshold))   cont_tight_match_maxpt++;
    }

    if ( (vbf_maxDeta_j1_isPileUpTight==true) && (vbf_maxDeta_j2_isPileUpTight==true)) {

      cont_tight_maxDeta++;

      if ((dR1_maxDeta <= dR_treshold) && (dR2_maxDeta <= dR_treshold) )   cont_tight_match_dR_maxDeta++;
      if ((dR1_maxDeta <= dR_treshold) && (dR2_maxDeta <= dR_treshold) && (fabs(ptratio1_maxDeta-1.) <= pt_treshold) && (fabs(ptratio2_maxDeta-1. )<= pt_treshold)) cont_tight_match_maxDeta++;
    }

    if ( (vbf_maxMjj_j1_isPileUpTight==true) && (vbf_maxMjj_j2_isPileUpTight==true)) {
      cont_tight_maxMjj++;

      if ((dR1_maxMjj <= dR_treshold) && (dR2_maxMjj <= dR_treshold) )   cont_tight_match_dR_maxMjj++;
      if ((dR1_maxMjj <= dR_treshold) && (dR2_maxMjj <= dR_treshold) && (fabs(ptratio1_maxMjj-1.) <= pt_treshold) && (fabs(ptratio2_maxMjj-1.) <= pt_treshold)) cont_tight_match_maxMjj++;
    }

    if ( (vbf_maxpt_j1_isPileUpMedium==true) && (vbf_maxpt_j2_isPileUpMedium==true)) {
      cont_medium_maxpt++;

      if ( (dR1_maxpt <= dR_treshold) && (dR2_maxpt <= dR_treshold))   cont_medium_match_dR_maxpt++;
      if ( (dR1_maxpt <= dR_treshold) && (dR2_maxpt <= dR_treshold) && (fabs(ptratio1_maxpt-1.) <= pt_treshold) && (fabs(ptratio2_maxpt-1.) <= pt_treshold)) cont_medium_match_maxpt++;
    }

    if ( (vbf_maxDeta_j1_isPileUpMedium==true) && (vbf_maxDeta_j2_isPileUpMedium==true)) {
      cont_medium_maxDeta++;

      if ((dR1_maxDeta <= dR_treshold) && (dR2_maxDeta <= dR_treshold))   cont_medium_match_dR_maxDeta++;
      if ((dR1_maxDeta <= dR_treshold) && (dR2_maxDeta <= dR_treshold) && (fabs(ptratio1_maxDeta-1.) <= pt_treshold) && (fabs(ptratio2_maxDeta-1.) <= pt_treshold)) cont_medium_match_maxDeta++;
    }

    if ( (vbf_maxMjj_j1_isPileUpMedium==true) && (vbf_maxMjj_j2_isPileUpMedium==true)) {
      cont_medium_maxMjj++;

      if ((dR1_maxMjj <= dR_treshold) && (dR2_maxMjj <= dR_treshold))   cont_medium_match_dR_maxMjj++;
      if ((dR1_maxMjj <= dR_treshold) && (dR2_maxMjj <= dR_treshold) && (fabs(ptratio1_maxMjj-1.) <= pt_treshold) && (fabs(ptratio2_maxMjj-1.) <= pt_treshold)) cont_medium_match_maxMjj++;
    }

  }

  std::cout<<"    "<<std::endl;
  std::cout<<" Input Events : "<<nEntry<<" Passed Jet Met Selection "<<nPassingEvents<<std::endl;  
  std::cout<<"    "<<std::endl;

  //create canvas and draw the plots  --> maxpt

  TCanvas *histo_dR1_maxpt_c = new TCanvas("histo_dR1_maxpt_c");
  histo_dR1_maxpt_c->cd();
  histo_dR1_maxpt_c->SetTickx();
  histo_dR1_maxpt_c->SetTicky();
  histo_dR1_maxpt->Draw();
  histo_dR1_maxpt_c->Print ((OutputDirectory+"/DeltaR1_maxpt.png").c_str(), "png");
  histo_dR1_maxpt->Write();

  TCanvas *histo_dR2_maxpt_c = new TCanvas("histo_dR2_maxpt_c");
  histo_dR2_maxpt_c->cd();
  histo_dR2_maxpt_c->SetTickx();
  histo_dR2_maxpt_c->SetTicky();
  histo_dR2_maxpt->Draw();
  histo_dR2_maxpt_c->Print ((OutputDirectory+"/DeltaR2_maxpt.png").c_str(), "png");
  histo_dR2_maxpt->Write();

  TCanvas *histo_deta1_maxpt_c = new TCanvas("histo_deta1_maxpt_c");
  histo_deta1_maxpt_c->cd();
  histo_deta1_maxpt_c->SetTickx();
  histo_deta1_maxpt_c->SetTicky();
  histo_deta1_maxpt->Draw();
  histo_deta1_maxpt_c->Print ((OutputDirectory+"/Deta1_maxpt.png").c_str(), "png");
  histo_deta1_maxpt->Write();

  TCanvas *histo_deta2_maxpt_c = new TCanvas("histo_deta2_maxpt_c");
  histo_deta2_maxpt_c->cd();
  histo_deta2_maxpt_c->SetTickx();
  histo_deta2_maxpt_c->SetTicky();
  histo_deta2_maxpt->Draw();
  histo_deta2_maxpt_c->Print ((OutputDirectory+"/Deta2_maxpt.png").c_str(), "png");
  histo_deta2_maxpt->Write();

  TCanvas *histo_dphi1_maxpt_c = new TCanvas("histo_dphi1_maxpt_c");
  histo_dphi1_maxpt_c->cd();
  histo_dphi1_maxpt_c->SetTickx();
  histo_dphi1_maxpt_c->SetTicky();
  histo_dphi1_maxpt->Draw();
  histo_dphi1_maxpt_c->Print ((OutputDirectory+"/Dphi1_maxpt.png").c_str(), "png");
  histo_dphi1_maxpt->Write();

  TCanvas *histo_dphi2_maxpt_c = new TCanvas("histo_dphi2_maxpt_c");
  histo_dphi2_maxpt_c->cd();
  histo_dphi2_maxpt_c->SetTickx();
  histo_dphi2_maxpt_c->SetTicky();
  histo_dphi2_maxpt->Draw();
  histo_dphi2_maxpt_c->Print ((OutputDirectory+"/Dphi2_maxpt.png").c_str(), "png");
  histo_dphi2_maxpt->Write();

  TCanvas *histo_dpt1_maxpt_c = new TCanvas("histo_dpt1_maxpt_c");
  histo_dpt1_maxpt_c->cd();
  histo_dpt1_maxpt_c->SetTickx();
  histo_dpt1_maxpt_c->SetTicky();
  histo_dpt1_maxpt->Draw();
  histo_dpt1_maxpt_c->Print ((OutputDirectory+"/Dpt1_maxpt.png").c_str(), "png");
  histo_dpt1_maxpt->Write();


  TCanvas *histo_dpt2_maxpt_c = new TCanvas("histo_dpt2_maxpt_c");
  histo_dpt2_maxpt_c->cd();
  histo_dpt2_maxpt_c->SetTickx();
  histo_dpt2_maxpt_c->SetTicky();
  histo_dpt2_maxpt->Draw();
  histo_dpt2_maxpt_c->Print ((OutputDirectory+"/Dpt2_maxpt.png").c_str(), "png");
  histo_dpt2_maxpt->Write();

  TCanvas *histo_dR1_dpt1_maxpt_c = new TCanvas("histo_dR1_dpt1_maxpt_c");
  histo_dR1_dpt1_maxpt_c->cd();
  histo_dR1_dpt1_maxpt_c->SetTickx();
  histo_dR1_dpt1_maxpt_c->SetTicky();
  histo_dR1_dpt1_maxpt->Draw("colz");
  histo_dR1_dpt1_maxpt_c->Print ((OutputDirectory+"/DR1_Dpt1_maxpt.png").c_str(), "png");
  histo_dR1_dpt1_maxpt->Write();

  TCanvas *histo_dR2_dpt2_maxpt_c = new TCanvas("histo_dR2_dpt2_maxpt_c");
  histo_dR2_dpt2_maxpt_c->cd();
  histo_dR2_dpt2_maxpt_c->SetTickx();
  histo_dR2_dpt2_maxpt_c->SetTicky();
  histo_dR2_dpt2_maxpt->Draw("colz");
  histo_dR2_dpt2_maxpt_c->Print ((OutputDirectory+"/DR2_Dpt2_maxpt.png").c_str(), "png");
  histo_dR2_dpt2_maxpt->Write();

  //create canvas and draw the plots  --> maxdeta

  TCanvas *histo_dR1_maxDeta_c = new TCanvas("histo_dR1_maxDeta_c");
  histo_dR1_maxDeta_c->cd();
  histo_dR1_maxDeta_c->SetTickx();
  histo_dR1_maxDeta_c->SetTicky();
  histo_dR1_maxDeta->Draw();
  histo_dR1_maxDeta_c->Print ((OutputDirectory+"/DeltaR1_maxDeta.png").c_str(), "png");
  histo_dR1_maxDeta->Write();

  TCanvas *histo_dR2_maxDeta_c = new TCanvas("histo_dR2_maxDeta_c");
  histo_dR2_maxDeta_c->cd();
  histo_dR2_maxDeta_c->SetTickx();
  histo_dR2_maxDeta_c->SetTicky();
  histo_dR2_maxDeta->Draw();
  histo_dR2_maxDeta_c->Print ((OutputDirectory+"/DeltaR2_maxDeta.png").c_str(), "png");
  histo_dR2_maxDeta->Write();

  TCanvas *histo_deta1_maxDeta_c = new TCanvas("histo_deta1_maxDeta_c");
  histo_deta1_maxDeta_c->cd();
  histo_deta1_maxDeta_c->SetTickx();
  histo_deta1_maxDeta_c->SetTicky();
  histo_deta1_maxDeta->Draw();
  histo_deta1_maxDeta_c->Print ((OutputDirectory+"/Deta1_maxDeta.png").c_str(), "png");
  histo_deta1_maxDeta->Write();

  TCanvas *histo_deta2_maxDeta_c = new TCanvas("histo_deta2_maxDeta_c");
  histo_deta2_maxDeta_c->cd();
  histo_deta2_maxDeta_c->SetTickx();
  histo_deta2_maxDeta_c->SetTicky();
  histo_deta2_maxDeta->Draw();
  histo_deta2_maxDeta_c->Print ((OutputDirectory+"/Deta2_maxDeta.png").c_str(), "png");
  histo_deta2_maxDeta->Write();

  TCanvas *histo_dphi1_maxDeta_c = new TCanvas("histo_dphi1_maxDeta_c");
  histo_dphi1_maxDeta_c->cd();
  histo_dphi1_maxDeta_c->SetTickx();
  histo_dphi1_maxDeta_c->SetTicky();
  histo_dphi1_maxDeta->Draw();
  histo_dphi1_maxDeta_c->Print ((OutputDirectory+"/Dphi1_maxDeta.png").c_str(), "png");
  histo_dphi1_maxDeta->Write();

  TCanvas *histo_dphi2_maxDeta_c = new TCanvas("histo_dphi2_maxDeta_c");
  histo_dphi2_maxDeta_c->cd();
  histo_dphi2_maxDeta_c->SetTickx();
  histo_dphi2_maxDeta_c->SetTicky();
  histo_dphi2_maxDeta->Draw();
  histo_dphi2_maxDeta_c->Print ((OutputDirectory+"/Dphi2_maxDeta.png").c_str(), "png");
  histo_dphi2_maxDeta->Write();

  TCanvas *histo_dpt1_maxDeta_c = new TCanvas("histo_dpt1_maxDeta_c");
  histo_dpt1_maxDeta_c->cd();
  histo_dpt1_maxDeta_c->SetTickx();
  histo_dpt1_maxDeta_c->SetTicky();
  histo_dpt1_maxDeta->Draw();
  histo_dpt1_maxDeta_c->Print ((OutputDirectory+"/Dpt1_maxDeta.png").c_str(), "png");
  histo_dpt1_maxDeta->Write();


  TCanvas *histo_dpt2_maxDeta_c = new TCanvas("histo_dpt2_maxDeta_c");
  histo_dpt2_maxDeta_c->cd();
  histo_dpt2_maxDeta_c->SetTickx();
  histo_dpt2_maxDeta_c->SetTicky();
  histo_dpt2_maxDeta->Draw();
  histo_dpt2_maxDeta_c->Print ((OutputDirectory+"/Dpt2_maxDeta.png").c_str(), "png");
  histo_dpt2_maxDeta->Write();

  TCanvas *histo_dR1_dpt1_maxDeta_c = new TCanvas("histo_dR1_dpt1_maxDeta_c");
  histo_dR1_dpt1_maxDeta_c->cd();
  histo_dR1_dpt1_maxDeta_c->SetTickx();
  histo_dR1_dpt1_maxDeta_c->SetTicky();
  histo_dR1_dpt1_maxDeta->Draw("colz");
  histo_dR1_dpt1_maxDeta_c->Print ((OutputDirectory+"/DR1_Dpt1_maxDeta.png").c_str(), "png");
  histo_dR1_dpt1_maxDeta->Write();

  TCanvas *histo_dR2_dpt2_maxDeta_c = new TCanvas("histo_dR2_dpt2_maxDeta_c");
  histo_dR2_dpt2_maxDeta_c->cd();
  histo_dR2_dpt2_maxDeta_c->SetTickx();
  histo_dR2_dpt2_maxDeta_c->SetTicky();
  histo_dR2_dpt2_maxDeta->Draw("colz");
  histo_dR2_dpt2_maxDeta_c->Print ((OutputDirectory+"/DR2_Dpt2_maxDeta.png").c_str(), "png");
  histo_dR2_dpt2_maxDeta->Write();

  //create canvas and draw the plots  --> maxMjj

  TCanvas *histo_dR1_maxMjj_c = new TCanvas("histo_dR1_maxMjj_c");
  histo_dR1_maxMjj_c->cd();
  histo_dR1_maxMjj_c->SetTickx();
  histo_dR1_maxMjj_c->SetTicky();
  histo_dR1_maxMjj->Draw();
  histo_dR1_maxMjj_c->Print ((OutputDirectory+"/DeltaR1_maxMjj.png").c_str(), "png");
  histo_dR1_maxMjj->Write();

  TCanvas *histo_dR2_maxMjj_c = new TCanvas("histo_dR2_maxMjj_c");
  histo_dR2_maxMjj_c->cd();
  histo_dR2_maxMjj_c->SetTickx();
  histo_dR2_maxMjj_c->SetTicky();
  histo_dR2_maxMjj->Draw();
  histo_dR2_maxMjj_c->Print ((OutputDirectory+"/DeltaR2_maxMjj.png").c_str(), "png");
  histo_dR2_maxMjj->Write();

  TCanvas *histo_deta1_maxMjj_c = new TCanvas("histo_deta1_maxMjj_c");
  histo_deta1_maxMjj_c->cd();
  histo_deta1_maxMjj_c->SetTickx();
  histo_deta1_maxMjj_c->SetTicky();
  histo_deta1_maxMjj->Draw();
  histo_deta1_maxMjj_c->Print ((OutputDirectory+"/Deta1_maxMjj.png").c_str(), "png");
  histo_deta1_maxMjj->Write();

  TCanvas *histo_deta2_maxMjj_c = new TCanvas("histo_deta2_maxMjj_c");
  histo_deta2_maxMjj_c->cd();
  histo_deta2_maxMjj_c->SetTickx();
  histo_deta2_maxMjj_c->SetTicky();
  histo_deta2_maxMjj->Draw();
  histo_deta2_maxMjj_c->Print ((OutputDirectory+"/Deta2_maxMjj.png").c_str(), "png");
  histo_deta2_maxMjj->Write();

  TCanvas *histo_dphi1_maxMjj_c = new TCanvas("histo_dphi1_maxMjj_c");
  histo_dphi1_maxMjj_c->cd();
  histo_dphi1_maxMjj_c->SetTickx();
  histo_dphi1_maxMjj_c->SetTicky();
  histo_dphi1_maxMjj->Draw();
  histo_dphi1_maxMjj_c->Print ((OutputDirectory+"/Dphi1_maxMjj.png").c_str(), "png");
  histo_dphi1_maxMjj->Write();

  TCanvas *histo_dphi2_maxMjj_c = new TCanvas("histo_dphi2_maxMjj_c");
  histo_dphi2_maxMjj_c->cd();
  histo_dphi2_maxMjj_c->SetTickx();
  histo_dphi2_maxMjj_c->SetTicky();
  histo_dphi2_maxMjj->Draw();
  histo_dphi2_maxMjj_c->Print ((OutputDirectory+"/Dphi2_maxMjj.png").c_str(), "png");
  histo_dphi2_maxMjj->Write();

  TCanvas *histo_dpt1_maxMjj_c = new TCanvas("histo_dpt1_maxMjj_c");
  histo_dpt1_maxMjj_c->cd();
  histo_dpt1_maxMjj_c->SetTickx();
  histo_dpt1_maxMjj_c->SetTicky();
  histo_dpt1_maxMjj->Draw();
  histo_dpt1_maxMjj_c->Print ((OutputDirectory+"/Dpt1_maxMjj.png").c_str(), "png");
  histo_dpt1_maxMjj->Write();


  TCanvas *histo_dpt2_maxMjj_c = new TCanvas("histo_dpt2_maxMjj_c");
  histo_dpt2_maxMjj_c->cd();
  histo_dpt2_maxMjj_c->SetTickx();
  histo_dpt2_maxMjj_c->SetTicky();
  histo_dpt2_maxMjj->Draw();
  histo_dpt2_maxMjj_c->Print ((OutputDirectory+"/Dpt2_maxMjj.png").c_str(), "png");
  histo_dpt2_maxMjj->Write();

  TCanvas *histo_dR1_dpt1_maxMjj_c = new TCanvas("histo_dR1_dpt1_maxMjj_c");
  histo_dR1_dpt1_maxMjj_c->cd();
  histo_dR1_dpt1_maxMjj_c->SetTickx();
  histo_dR1_dpt1_maxMjj_c->SetTicky();
  histo_dR1_dpt1_maxMjj->Draw("colz");
  histo_dR1_dpt1_maxMjj_c->Print ((OutputDirectory+"/DR1_Dpt1_maxMjj.png").c_str(), "png");
  histo_dR1_dpt1_maxMjj->Write();

  TCanvas *histo_dR2_dpt2_maxMjj_c = new TCanvas("histo_dR2_dpt2_maxMjj_c");
  histo_dR2_dpt2_maxMjj_c->cd();
  histo_dR2_dpt2_maxMjj_c->SetTickx();
  histo_dR2_dpt2_maxMjj_c->SetTicky();
  histo_dR2_dpt2_maxMjj->Draw("colz");
  histo_dR2_dpt2_maxMjj_c->Print ((OutputDirectory+"/DR2_Dpt2_maxMjj.png").c_str(), "png");
  histo_dR2_dpt2_maxMjj->Write();

  OutputFile->Close();
 
  //print information
  std::cout<<"                                       "<<std::endl;
  std::cout<<"Event Fraction: (Total)     (thresholds dR=0.3,  ptratio=0.3)"<<std::endl;

  std::cout<<"######################################################################"<<std::endl;

  std::cout<<"Pt max: both tag jets in dR  = "<<1.*cont_maxpt2/(1.*nPassingEvents)<<std::endl;
  std::cout<<"Pt max: only 1 jet in dR     = "<<1.*cont_maxpt1/(1.*nPassingEvents)<<std::endl;
  std::cout<<"Pt max: no matching in dR    = "<<1.*cont_maxpt0/(1.*nPassingEvents)<<std::endl;
  std::cout<<"Pt max: Both in dR and dPt   = "<<1.*cont_dRandpt_maxpt/(1.*nPassingEvents)<<std::endl;

  std::cout<<"#######################################################################"<<std::endl;
 
  std::cout<<"Deta max: both tag jets in dR =  "<<1.*cont_maxDeta2/(1.*nPassingEvents)<<std::endl;
  std::cout<<"Deta max: only one jets in dR =  "<<1.*cont_maxDeta1/(1.*nPassingEvents)<<std::endl;
  std::cout<<"Deta max: no matching         =  "<<1.*cont_maxDeta0/(1.*nPassingEvents)<<std::endl;
  std::cout<<"Deta max: Both in dR and dPt  =  "<<1.*cont_dRandpt_maxDeta/(1.*nPassingEvents)<<std::endl;

  std::cout<<"#######################################################################"<<std::endl;

  std::cout<<"Mjj max: both tag jets in dR =  "<<1.*cont_maxMjj2/(1.*nPassingEvents)<<std::endl;
  std::cout<<"Mjj max: only one jets in dR =  "<<1.*cont_maxMjj1/(1.*nPassingEvents)<<std::endl;
  std::cout<<"Mjj max: no matching         =  "<<1.*cont_maxMjj0/(1.*nPassingEvents)<<std::endl;
  std::cout<<"Mjj max: Both in dR and dPt  =  "<<1.*cont_dRandpt_maxMjj/(1.*nPassingEvents)<<std::endl;


  std::cout<<"#######################################################################"<<std::endl;

  std::cout<<"Maxpt matching in dR and dPt --> There is a matching also in maxDeta "<<1.*cont_maxpt_maxDeta/(1.*cont_dRandpt_maxpt)<<std::endl; 
  std::cout<<"Maxpt matching in dR and dPt --> There is a matching also in maxMjj  "<<1.*cont_maxpt_maxMjj/(1.*cont_dRandpt_maxpt)<<std::endl; 
  std::cout<<"MaxDeta matching in dR and dPt --> There is a matching also in maxpt  "<<1.*cont_maxpt_maxDeta/(1.*cont_dRandpt_maxDeta)<<std::endl;
  std::cout<<"MaxDeta matching in dR and dPt --> There is a matching also in maxpt  "<<1.*cont_maxMjj_maxDeta/(1.*cont_dRandpt_maxDeta)<<std::endl; 
  std::cout<<"MaxMjj matching in dR and dPt --> There is a matching also in maxpt  "<<1.*cont_maxpt_maxMjj/(1.*cont_dRandpt_maxMjj)<<std::endl; 
  std::cout<<"MaxMjj matching in dR and dPt --> There is a matching also in maxDeta   "<<1.*cont_maxMjj_maxDeta/(1.*cont_dRandpt_maxMjj)<<std::endl; 

  std::cout<<"#######################################################################"<<std::endl;
 
  std::cout<<"Events with 3 jets: "<<std::endl;

  std::cout<<"#######################################################################"<<std::endl;

  std::cout<<"Pt max: both in dR   "<<1.*cont_nj_match_dR_maxpt[3]/(1.*cont_nj[3])<<std::endl;
  std::cout<<"Deta max: both in dR "<<1.*cont_nj_match_dR_maxDeta[3]/(1.*cont_nj[3])<<std::endl;
  std::cout<<"Mjj max: both in dR  "<<1.*cont_nj_match_dR_maxMjj[3]/(1.*cont_nj[3])<<std::endl;

  std::cout<<"Pt max: both in dR and dpt   "<<1.*cont_nj_match_maxpt[3]/(1.*cont_nj[3])<<std::endl;
  std::cout<<"Deta max: both in dR and dpt "<<1.*cont_nj_match_maxDeta[3]/(1.*cont_nj[3])<<std::endl;
  std::cout<<"Mjj max: both in dR and dpt  "<<1.*cont_nj_match_maxMjj[3]/(1.*cont_nj[3])<<std::endl;


  std::cout<<"#######################################################################"<<std::endl;
 
  std::cout<<"Events with 4 jets: "<<std::endl;

  std::cout<<"#######################################################################"<<std::endl;

  std::cout<<"Pt max: both in dR   "<<1.*cont_nj_match_dR_maxpt[4]/(1.*cont_nj[4])<<std::endl;
  std::cout<<"Deta max: both in dR "<<1.*cont_nj_match_dR_maxDeta[4]/(1.*cont_nj[4])<<std::endl;
  std::cout<<"Mjj max: both in dR  "<<1.*cont_nj_match_dR_maxMjj[4]/(1.*cont_nj[4])<<std::endl;

  std::cout<<"Pt max: both in dR and dpt   "<<1.*cont_nj_match_maxpt[4]/(1.*cont_nj[4])<<std::endl;
  std::cout<<"Deta max: both in dR and dpt "<<1.*cont_nj_match_maxDeta[4]/(1.*cont_nj[4])<<std::endl;
  std::cout<<"Mjj max: both in dR and dpt  "<<1.*cont_nj_match_maxMjj[4]/(1.*cont_nj[5])<<std::endl;


  std::cout<<"#######################################################################"<<std::endl;
 
  std::cout<<"Events with 5 jets: "<<std::endl;

  std::cout<<"#######################################################################"<<std::endl;

  std::cout<<"Pt max: both in dR   "<<1.*cont_nj_match_dR_maxpt[5]/(1.*cont_nj[5])<<std::endl;
  std::cout<<"Deta max: both in dR "<<1.*cont_nj_match_dR_maxDeta[5]/(1.*cont_nj[5])<<std::endl;
  std::cout<<"Mjj max: both in dR  "<<1.*cont_nj_match_dR_maxMjj[5]/(1.*cont_nj[5])<<std::endl;

  std::cout<<"Pt max: both in dR and dpt   "<<1.*cont_nj_match_maxpt[5]/(1.*cont_nj[5])<<std::endl;
  std::cout<<"Deta max: both in dR and dpt "<<1.*cont_nj_match_maxDeta[5]/(1.*cont_nj[5])<<std::endl;
  std::cout<<"Mjj max: both in dR and dpt  "<<1.*cont_nj_match_maxMjj[5]/(1.*cont_nj[5])<<std::endl;


  std::cout<<"#######################################################################"<<std::endl;

  std::cout<<"Evens with both vbf jets satisify medium pile-up : "<<std::endl;

  std::cout<<"#######################################################################"<<std::endl;

  std::cout<<"Pt max: both in dR   "<<1.*cont_medium_match_dR_maxpt/(1.*cont_medium_maxpt)<<std::endl;
  std::cout<<"Deta max: both in dR "<<1.*cont_medium_match_dR_maxDeta/(1.*cont_medium_maxpt)<<std::endl;
  std::cout<<"Mjj max: both in dR  "<<1.*cont_medium_match_dR_maxMjj/(1.*cont_medium_maxpt)<<std::endl;

  std::cout<<"Pt max: both in dR and dpt   "<<1.*cont_medium_match_maxpt/(1.*cont_medium_maxpt)<<std::endl;
  std::cout<<"Deta max: both in dR and dpt "<<1.*cont_medium_match_maxDeta/(1.*cont_medium_maxpt)<<std::endl;
  std::cout<<"Mjj max: both in dR  and dpt "<<1.*cont_medium_match_maxMjj/(1.*cont_medium_maxpt)<<std::endl;

  std::cout<<"#######################################################################"<<std::endl;

  std::cout<<"Evens with both vbf jets satisify tight pile-up : "<<std::endl;

  std::cout<<"#######################################################################"<<std::endl;

  std::cout<<"Pt max: both in dR   "<<1.*cont_tight_match_dR_maxpt/(1.*cont_tight_maxpt)<<std::endl;
  std::cout<<"Deta max: both in dR "<<1.*cont_tight_match_dR_maxDeta/(1.*cont_tight_maxpt)<<std::endl;
  std::cout<<"Mjj max: both in dR  "<<1.*cont_tight_match_dR_maxMjj/(1.*cont_tight_maxpt)<<std::endl;

  std::cout<<"Pt max: both in dR and dpt   "<<1.*cont_tight_match_maxpt/(1.*cont_tight_maxpt)<<std::endl;
  std::cout<<"Deta max: both in dR and dpt "<<1.*cont_tight_match_maxDeta/(1.*cont_tight_maxpt)<<std::endl;
  std::cout<<"Mjj max: both in dR  and dpt "<<1.*cont_tight_match_maxMjj/(1.*cont_tight_maxpt)<<std::endl;


  std::cout<<"#######################################################################"<<std::endl;

  std::cout<<"Fraction of events with 2 jets  pileup medium: "<<std::endl;

  std::cout<<"#######################################################################"<<std::endl;

  std::cout<<"Pt max: "<<1.*cont_medium_maxpt/(1.*nPassingEvents)<<std::endl;
  std::cout<<"Deta max: "<<1.*cont_medium_maxDeta/(1.*nPassingEvents)<<std::endl;
  std::cout<<"Mjj max: "<<1.*cont_medium_maxMjj/(1.*nPassingEvents)<<std::endl;

  std::cout<<"#######################################################################"<<std::endl;

  std::cout<<"Fraction of events with 2 jets  pileup medium: "<<std::endl;

  std::cout<<"#######################################################################"<<std::endl;

  std::cout<<"Pt max: "<<1.*cont_tight_maxpt/(1.*nPassingEvents)<<std::endl;
  std::cout<<"Deta max: "<<1.*cont_tight_maxDeta/(1.*nPassingEvents)<<std::endl;
  std::cout<<"Mjj max: "<<cont_tight_maxMjj/(1.*nPassingEvents)<<std::endl;


  std::cout<<"                                       "<<std::endl;
  std::cout<<"Finish : exit from Programme          "<<std::endl;
  std::cout<<"                                       "<<std::endl;
  
  return(0);

}
										       

//  LocalWords:  endl
