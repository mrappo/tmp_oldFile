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
  std::string InputFileName      = gConfigParser -> readStringOption("Input::inputFileName");

  std::cout<<"         "<<std::endl;
  std::cout<<" InputDirectory :    "<<InputDirectory<<std::endl;
  std::cout<<"         "<<std::endl;

  std::cout<<" InputFileName :    "<<InputFileName<<std::endl;
  std::cout<<"         "<<std::endl;

  
  std::string TreeName ;
  try{ TreeName = gConfigParser -> readStringOption("Input::TreeName"); }
  catch( char const* exceptionString){ TreeName = "WJet" ;
                                       std::cerr<<" TreeName set by default to --> WJet" <<std::endl;
  }

  float BoostedPtCut =  gConfigParser -> readFloatOption("Input::BoostedPtCut");
  float MetCut       =  gConfigParser -> readFloatOption("Input::MetCut");

  std::cout<<" TreeName :    "<<TreeName<<std::endl;
  std::cout<<"         "<<std::endl;

  float dR_treshold ;
  try{ dR_treshold = gConfigParser -> readFloatOption("Option::DRCut"); }
  catch(char const* exceptionString){ dR_treshold = 0.3 ;
                                      std::cerr<<" DR Cut Set To --> 0.3"<<std::endl;
  }

  std::cout<<" dR_treshold :    "<<dR_treshold<<std::endl;
  std::cout<<"         "<<std::endl;

  float pt_treshold ;
  try{ pt_treshold = gConfigParser -> readFloatOption("Input::PtCut"); }
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

  std::cout<<" DRMin :    "<<DRMin<<" DRMax : "<<DRMax<<std::endl;
  std::cout<<"         "<<std::endl;

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


  std::cout<<" DPtMin :    "<<DPtMin<<" DPtMax"<<DPtMax<<std::endl;
  std::cout<<"         "<<std::endl;




  // Open Input File and Get Input Tree

  TFile * inputFile  = new TFile((InputDirectory+"/"+InputFileName).c_str(),"READ");

  if(inputFile==0){ std::cerr<<" No Input File found --> exit "<<std::endl; return -1 ;}
  
  TTree * fTree      = (TTree*) inputFile -> Get(TreeName.c_str());

  if(fTree==0){ std::cerr<<" No Input Tree found --> exit "<<std::endl;  return -1; }


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

  std::vector<int> cont_nj (10,0) ;

  std::vector<int> cont_nj_match_maxpt (10,0); 
  std::vector<int> cont_nj_match_maxDeta (10,0); 
  std::vector<int> cont_nj_match_maxMjj (10,0);  
  std::vector<int> cont_nj_match_dR_maxpt (10,0);
  std::vector<int> cont_nj_match_dR_maxDeta (10,0);
  std::vector<int> cont_nj_match_dR_maxMjj (10,0) ;
/*
  int cont_tight_maxpt=0, cont_tight_maxDeta=0, cont_tight_maxMjj=0;
  int cont_tight_match_maxpt=0, cont_tight_match_maxDeta=0, cont_tight_match_maxMjj=0;
  int cont_medium_maxpt=0, cont_medium_maxDeta=0, cont_medium_maxMjj=0;
  int cont_medium_match_maxpt=0, cont_medium_match_maxDeta=0, cont_medium_match_maxMjj=0;
  */

  int nPassingEvents =0;
  
  //start reading tree

  for(int iEntry =0; iEntry<nEntry ; iEntry++){


    dR11_maxpt   = 0., dR12_maxpt   = 0.; dR11_maxDeta = 0., dR12_maxDeta = 0.; dR11_maxMjj  = 0., dR12_maxMjj  = 0.; 
    dR21_maxpt   = 0., dR22_maxpt   = 0.; dR21_maxDeta = 0., dR22_maxDeta = 0.; dR21_maxMjj  = 0., dR22_maxMjj  = 0.;

    cont_maxpt0 =0, cont_maxDeta0 =0, cont_maxMjj0 =0; cont_maxpt1 =0, cont_maxDeta1 =0, cont_maxMjj1 =0;
    cont_maxpt2 =0, cont_maxDeta2 =0, cont_maxMjj2 =0; dR1_maxpt =0., dR1_maxDeta =0., dR1_maxMjj =0.;
    dR2_maxpt =0., dR2_maxDeta =0., dR2_maxMjj =0.;
  
    ptratio1_maxpt=0., ptratio1_maxDeta=0., ptratio1_maxMjj=0.; ptratio2_maxpt=0., ptratio2_maxDeta=0., ptratio2_maxMjj=0.;
    cont_dRandpt_maxpt = 0, cont_dRandpt_maxDeta=0, cont_dRandpt_maxMjj=0; cont_maxpt_maxDeta = 0, cont_maxpt_maxMjj=0, cont_maxMjj_maxDeta =0;

 
    if( iEntry%1000 == 0 ) std::cout << "reading saved entry " << iEntry << "\r" << std::endl;
 
    fTree->GetEntry(iEntry);
    
    if(GroomedJet_CA8_pt[0]<BoostedPtCut || W_pt < BoostedPtCut || event_met_pfmet < MetCut ) continue;
    nPassingEvents++;

    //calculation of dR for max pt
        
    dR11_maxpt = deltaR(W_TagQuark_phi[0],vbf_maxpt_j1_phi,W_TagQuark_eta[0],vbf_maxpt_j1_eta); 
    dR12_maxpt = deltaR(W_TagQuark_phi[0],vbf_maxpt_j2_phi,W_TagQuark_eta[0],vbf_maxpt_j2_eta); 
    dR21_maxpt = deltaR(W_TagQuark_phi[1],vbf_maxpt_j1_phi,W_TagQuark_eta[1],vbf_maxpt_j1_eta); 
    dR22_maxpt = deltaR(W_TagQuark_phi[1],vbf_maxpt_j2_phi,W_TagQuark_eta[1],vbf_maxpt_j2_eta); 

    if(dR11_maxpt < dR12_maxpt){
   
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
    else{

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

    if ( ((dR1_maxpt   <= dR_treshold) && (dR2_maxpt   > dR_treshold)) || ((dR1_maxpt   >  dR_treshold) && (dR2_maxpt   <= dR_treshold)) )   cont_maxpt1++;
    if ( ((dR1_maxDeta <= dR_treshold) && (dR2_maxDeta > dR_treshold)) || ((dR1_maxDeta <= dR_treshold) && (dR2_maxDeta <= dR_treshold)) )   cont_maxDeta1++;
    if ( ((dR1_maxMjj  <= dR_treshold) && (dR2_maxMjj  > dR_treshold)) || ((dR1_maxMjj  <= dR_treshold) && (dR2_maxMjj  <= dR_treshold)) )   cont_maxMjj1++;

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



    /*
    int arejetsflipped_maxpt=0, arejetsflipped_maxDeta=0, arejetsflipped_maxMjj=0;





    float deta1_maxpt, deta2_maxpt, dpt1_maxpt, dpt2_maxpt;

    //flip dR if they are exchanged

 
    //if conditions are respected, then increases the counters

    //only 1 jets has matching in dR


    //matching in dR and ptratio for events with 3,4,5,6,7.. jets

    //matching with pileUp conditions
    if ( (vbf_maxpt_j1_isPileUpTight==true) && (vbf_maxpt_j2_isPileUpTight==true)) {
      cont_tight_maxpt++;
      if ((dR1_maxpt<=dR_treshold)&&(dR2_maxpt<=dR_treshold)&&(ptratio1_maxpt<=pt_treshold)&&(ptratio2_maxpt<=pt_treshold))   cont_tight_match_maxpt++;
  }

    if ( (vbf_maxDeta_j1_isPileUpTight==true) && (vbf_maxDeta_j2_isPileUpTight==true)) {
      cont_tight_maxDeta++;
      if ((dR1_maxDeta<=dR_treshold)&&(dR2_maxDeta<=dR_treshold)&&(ptratio1_maxDeta<=pt_treshold)&&(ptratio2_maxDeta<=pt_treshold))   cont_tight_match_maxDeta++;
  }

    if ( (vbf_maxMjj_j1_isPileUpTight==true) && (vbf_maxMjj_j2_isPileUpTight==true)) {
      cont_tight_maxMjj++;
      if ((dR1_maxMjj<=dR_treshold)&&(dR2_maxMjj<=dR_treshold)&&(ptratio1_maxMjj<=pt_treshold)&&(ptratio2_maxMjj<=pt_treshold))   cont_tight_match_maxMjj++;
  }

    if ( (vbf_maxpt_j1_isPileUpMedium==true) && (vbf_maxpt_j2_isPileUpMedium==true)) {
      cont_medium_maxpt++;
      if ((dR1_maxpt<=dR_treshold)&&(dR2_maxpt<=dR_treshold)&&(ptratio1_maxpt<=pt_treshold)&&(ptratio2_maxpt<=pt_treshold))   cont_medium_match_maxpt++;
  }

    if ( (vbf_maxDeta_j1_isPileUpMedium==true) && (vbf_maxDeta_j2_isPileUpMedium==true)) {
      cont_medium_maxDeta++;
      if ((dR1_maxDeta<=dR_treshold)&&(dR2_maxDeta<=dR_treshold)&&(ptratio1_maxDeta<=pt_treshold)&&(ptratio2_maxDeta<=pt_treshold))   cont_medium_match_maxDeta++;
  }

    if ( (vbf_maxMjj_j1_isPileUpMedium==true) && (vbf_maxMjj_j2_isPileUpMedium==true)) {
      cont_medium_maxMjj++;
      if ((dR1_maxMjj<=dR_treshold)&&(dR2_maxMjj<=dR_treshold)&&(ptratio1_maxMjj<=pt_treshold)&&(ptratio2_maxMjj<=pt_treshold))   cont_medium_match_maxMjj++;
  }

    if ( (flag_dRandpt_maxpt==1)&&(flag_dRandpt_maxDeta==1) )  cont_maxpt_maxDeta++;
    if ( (flag_dRandpt_maxpt==1)&&(flag_dRandpt_maxMjj==1) )   cont_maxpt_maxMjj++;
    if ( (flag_dRandpt_maxMjj==1)&&(flag_dRandpt_maxDeta==1) ) cont_maxMjj_maxDeta++;
    */    
  }
  
  std::cout<<" Input Events : "<<nEntry<<" Passed Jet Met Selection "<<nPassingEvents<<std::endl;  
  //create canvas and draw
  /*
  TCanvas *c1 = new TCanvas("c1","Graph Draw Options",200,10,600,400);
  c1->cd();
  c1->SetTickx();
  c1->SetTicky();
  histo_dR1_maxpt->GetXaxis()->SetTitle("#DeltaR_(j_{1},q_{1})");
  histo_dR1_maxpt->Draw();
  c1->Print ("DeltaR1_maxpt.png", "png");
  delete c1;

  TCanvas *c2 = new TCanvas("c2","Graph Draw Options",200,10,600,400);
  c2->cd();
  c2->SetTickx();
  c2->SetTicky();
  histo_dR2_maxpt->GetXaxis()->SetTitle("#DeltaR_(j_{2},q_{2})");
  histo_dR2_maxpt->Draw();
  c2->Print ("DeltaR2_maxpt.png", "png");
  delete c2;

  TCanvas *c3 = new TCanvas("c3","Graph Draw Options",200,10,600,400);
  c3->cd();
  c3->SetTickx();
  c3->SetTicky();
  histo_deta1_maxpt->GetXaxis()->SetTitle("#Delta#eta_(j_{1},q_{1})");
  histo_deta1_maxpt->Draw();
  c3->Print ("DeltaEta1_maxpt.png", "png");
  delete c3;

  TCanvas *c4 = new TCanvas("c4","Graph Draw Options",200,10,600,400);
  c4->cd();
  c4->SetTickx();
  c4->SetTicky();
  histo_deta2_maxpt->GetXaxis()->SetTitle("#Delta#eta_(j_{2},q_{2})");
  histo_deta2_maxpt->Draw();
  c4->Print ("DeltaEta2_maxpt.png", "png");
  delete c4;

  TCanvas *c5 = new TCanvas("c5","Graph Draw Options",200,10,600,400);
  c5->cd();
  c5->SetTickx();
  c5->SetTicky();
  histo_dpt1_maxpt->GetXaxis()->SetTitle("#Deltap_{T}_(j_{1},q_{1})");
  histo_dpt1_maxpt->Draw();
  c5->Print ("DeltaPt1_maxpt.png", "png");
  delete c5;

  TCanvas *c6 = new TCanvas("c6","Graph Draw Options",200,10,600,400);
  c6->cd();
  c6->SetTickx();
  c6->SetTicky();
  histo_dpt2_maxpt->GetXaxis()->SetTitle("#Deltap_T_(j_{2},q_{2})");
  histo_dpt2_maxpt->Draw();
  c6->Print ("DeltaPt2_maxpt.png", "png");
  delete c6;
  
  //print information

  cout<<"Frazioni di eventi in cui: (totale)     (soglia dR=0.3, soglia pt=0.3)"<<endl;

  cout<<"######################################################################"<<endl;

  cout<<"Pt max: entrambi i jet matchano in dR "<<1.*cont_maxpt2/(1.*nEntry)<<endl;
  cout<<"Pt max: un solo jet matcha in dR "<<1.*cont_maxpt1/(1.*nEntry)<<endl;
  cout<<"Pt max: nessun jet matcha in dR "<<1.*cont_maxpt0/(1.*nEntry)<<endl;
  cout<<"Pt max: entrambi i jet matchano sia in dR sia in pt: "<<1.*cont_dRandpt_maxpt/(1.*nEntry)<<endl;

  cout<<"#######################################################################"<<endl;

  cout<<"Deta max: entrambi i jet matchano in dR "<<1.*cont_maxDeta2/(1.*nEntry)<<endl;
  cout<<"Deta max: un solo jet matcha in dR "<<1.*cont_maxDeta1/(1.*nEntry)<<endl;
  cout<<"Deta max: nessun jet matcha in dR "<<1.*cont_maxDeta0/(1.*nEntry)<<endl;
  cout<<"Deta max: entrambi i jet matchano sia in dR sia in pt: "<<1.*cont_dRandpt_maxDeta/(1.*nEntry)<<endl;

  cout<<"#######################################################################"<<endl;

  cout<<"Mjj max: entrambi i jet matchano in dR "<<1.*cont_maxMjj2/(1.*nEntry)<<endl; 
  cout<<"Mjj max: un solo jet matcha in dR "<<1.*cont_maxMjj1/(1.*nEntry)<<endl; 
  cout<<"Mjj max: nessun jet matcha dR "<<1.*cont_maxMjj0/(1.*nEntry)<<endl;
  cout<<"Mjj max: entrambi i jet matchano sia in dR sia in pt: "<<1.*cont_dRandpt_maxMjj/(1.*nEntry)<<endl;


  cout<<"#######################################################################"<<endl;

  cout<<"Degli eventi in cui matchano (in dR e pt) i 2 jet a maxpt, c'è matching anche per quelli in Deta "<<1.*cont_maxpt_maxDeta/(1.*cont_dRandpt_maxpt)<<endl; 
  cout<<"Degli eventi in cui matchano (in dR e pt) i 2 jet a maxpt, c'è matching anche per quelli in Mjj "<<1.*cont_maxpt_maxMjj/(1.*cont_dRandpt_maxpt)<<endl; 
  cout<<"Degli eventi in cui matchano (in dR e pt) i 2 jet a maxDeta, c'è matching anche per quelli in pt  "<<1.*cont_maxpt_maxDeta/(1.*cont_dRandpt_maxDeta)<<endl;
  cout<<"Degli eventi in cui matchano (in dR e pt) i 2 jet a maxDeta, c'è matching anche per quelli in Mjj "<<1.*cont_maxMjj_maxDeta/(1.*cont_dRandpt_maxDeta)<<endl; 
  cout<<"Degli eventi in cui matchano (in dR e pt) i 2 jet a maxMjj, c'è matching anche per quelli in pt "<<1.*cont_maxpt_maxMjj/(1.*cont_dRandpt_maxMjj)<<endl; 
  cout<<"Degli eventi in cui matchano (in dR e pt) i 2 jet a maxMjj, c'è matching anche per quelli in Deta  "<<1.*cont_maxMjj_maxDeta/(1.*cont_dRandpt_maxMjj)<<endl; 

  cout<<"#######################################################################"<<endl;
 
  cout<<"Eventi a 3 jets: "<<endl;

  cout<<"#######################################################################"<<endl;

  cout<<"Pt max: entrambi i jets matchano in dR e pt "<<1.*cont_nj_match_maxpt[3]/(1.*cont_nj[3])<<endl;
  cout<<"Deta max: entrambi i jets matchano in dR e pt "<<1.*cont_nj_match_maxDeta[3]/(1.*cont_nj[3])<<endl;
  cout<<"Mjj max: entrambi i jets matchano in dR e pt "<<1.*cont_nj_match_maxMjj[3]/(1.*cont_nj[3])<<endl;

  cout<<"#######################################################################"<<endl;

  cout<<"#######################################################################"<<endl;
 
  cout<<"Eventi a 3 jets: "<<endl;

  cout<<"#######################################################################"<<endl;

  cout<<"Pt max: entrambi i jets matchano in dR "<<1.*cont_nj_match_dR_maxpt[3]/(1.*cont_nj[3])<<endl;
  cout<<"Deta max: entrambi i jets matchano in dR "<<1.*cont_nj_match_dR_maxDeta[3]/(1.*cont_nj[3])<<endl;
  cout<<"Mjj max: entrambi i jets matchano in dR "<<1.*cont_nj_match_dR_maxMjj[3]/(1.*cont_nj[3])<<endl;

  cout<<"#######################################################################"<<endl;

  cout<<"Eventi a 4 jets: "<<endl;

  cout<<"#######################################################################"<<endl;

  cout<<"Pt max: entrambi i jets matchano in dR e pt "<<1.*cont_nj_match_maxpt[4]/(1.*cont_nj[4])<<endl;
  cout<<"Deta max: entrambi i jets matchano in dR e pt "<<1.*cont_nj_match_maxDeta[4]/(1.*cont_nj[4])<<endl;
  cout<<"Mjj max: entrambi i jets matchano in dR e pt "<<1.*cont_nj_match_maxMjj[4]/(1.*cont_nj[4])<<endl;

  cout<<"#######################################################################"<<endl;

  cout<<"Eventi a 5 jets: "<<endl;

  cout<<"#######################################################################"<<endl;

  cout<<"Pt max: entrambi i jets matchano in dR e pt "<<1.*cont_nj_match_maxpt[5]/(1.*cont_nj[5])<<endl;
  cout<<"Deta max: entrambi i jets matchano in dR e pt "<<1.*cont_nj_match_maxDeta[5]/(1.*cont_nj[5])<<endl;
  cout<<"Mjj max: entrambi i jets matchano in dR e pt "<<1.*cont_nj_match_maxMjj[5]/(1.*cont_nj[5])<<endl;

  cout<<"#######################################################################"<<endl;

  cout<<"Eventi a 6 jets: "<<endl;

  cout<<"#######################################################################"<<endl;

  cout<<"Pt max: entrambi i jets matchano in dR e pt "<<1.*cont_nj_match_maxpt[6]/(1.*cont_nj[6])<<endl;
  cout<<"Deta max: entrambi i jets matchano in dR e pt "<<1.*cont_nj_match_maxDeta[6]/(1.*cont_nj[6])<<endl;
  cout<<"Mjj max: entrambi i jets matchano in dR e pt "<<1.*cont_nj_match_maxMjj[6]/(1.*cont_nj[6])<<endl;

  cout<<"#######################################################################"<<endl;

  cout<<"Eventi a 7 jets: "<<endl;

  cout<<"#######################################################################"<<endl;

  cout<<"Pt max: entrambi i jets matchano in dR e pt "<<1.*cont_nj_match_maxpt[7]/(1.*cont_nj[7])<<endl;
  cout<<"Deta max: entrambi i jets matchano in dR e pt "<<1.*cont_nj_match_maxDeta[7]/(1.*cont_nj[7])<<endl;
  cout<<"Mjj max: entrambi i jets matchano in dR e pt "<<1.*cont_nj_match_maxMjj[7]/(1.*cont_nj[7])<<endl;

  cout<<"#######################################################################"<<endl;

  cout<<"Eventi in cui entrambi i jets hanno passato la medium pileup: "<<endl;

  cout<<"#######################################################################"<<endl;

  cout<<"Pt max: entrambi i jets matchano in dR e pt "<<1.*cont_medium_match_maxpt/(1.*cont_medium_maxpt)<<endl;
  cout<<"Deta max: entrambi i jets matchano in dR e pt "<<1.*cont_medium_match_maxDeta/(1.*cont_medium_maxDeta)<<endl;
  cout<<"Mjj max: entrambi i jets matchano in dR e pt "<<1.*cont_medium_match_maxMjj/(1.*cont_medium_maxMjj)<<endl;

  cout<<"#######################################################################"<<endl;

  cout<<"Eventi in cui entrambi i jets hanno passato la tight pileup: "<<endl;

  cout<<"#######################################################################"<<endl;

  cout<<"Pt max: entrambi i jets matchano in dR e pt "<<1.*cont_tight_match_maxpt/(1.*cont_tight_maxpt)<<endl;
  cout<<"Deta max: entrambi i jets matchano in dR e pt "<<1.*cont_tight_match_maxDeta/(1.*cont_tight_maxDeta)<<endl;
  cout<<"Mjj max: entrambi i jets matchano in dR e pt "<<1.*cont_tight_match_maxMjj/(1.*cont_tight_maxMjj)<<endl;

  cout<<"#######################################################################"<<endl;

  cout<<"Frazioni di eventi sul totale in cui i 2 jet han passato la pileup medium: "<<endl;

  cout<<"#######################################################################"<<endl;

  cout<<"Pt max: "<<1.*cont_medium_maxpt/(1.*nEntry)<<endl;
  cout<<"Deta max: "<<1.*cont_medium_maxDeta/(1.*nEntry)<<endl;
  cout<<"Mjj max: "<<1.*cont_medium_maxMjj/(1.*nEntry)<<endl;

  cout<<"#######################################################################"<<endl;

  cout<<"Frazioni di eventi sul totale in cui i 2 jet han passato la pileup tight: "<<endl;

  cout<<"#######################################################################"<<endl;

  cout<<"Pt max: "<<1.*cont_tight_maxpt/(1.*nEntry)<<endl;
  cout<<"Deta max: "<<1.*cont_tight_maxDeta/(1.*nEntry)<<endl;
  cout<<"Mjj max: "<<cont_tight_maxMjj/(1.*nEntry)<<endl;

  */
  return(0);}
										       
