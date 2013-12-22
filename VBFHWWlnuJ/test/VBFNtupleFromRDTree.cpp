
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <functional>
#include <assert.h>
#include <map>

#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TCanvas.h"

#include "ntpleUtils.h"
#include "ConfigParser.h"

#include "JetCollectionSorting.h"
#include "METzCalculator.h"
#include "VBFAnalysisUtils.h"

#include "VBFMuonClass.h"
#include "VBFElectronClass.h"

/// Main Programme

int main (int argc, char** argv){

  if (argc != 2){

      std::cerr << ">>> Usage:   " << argv[1] << "   treeFile.root" << std::endl;
      return -1;
    }

  float Wmass = 80.385;

  // parse config file parameter                                                                                                                                                                
  parseConfigFile(argv[1]);

  std::string InputDirectory     = gConfigParser -> readStringOption("Input::InputDirectory");
  std::string InputRootFile      = gConfigParser -> readStringOption("Input::InputRootFile"); 
  std::string TreeName           = gConfigParser -> readStringOption("Input::TreeName");
  std::string LeptonType         = gConfigParser -> readStringOption("Input::LeptonType");
  
  std::cout<<"                   "<<std::endl;
  std::cout<<" Input Directory   "<<InputDirectory<<std::endl;
  std::cout<<" Input Root File   "<<InputRootFile<<std::endl;
  std::cout<<" Input TreeName    "<<TreeName<<std::endl;
  std::cout<<" Input Lepton Type "<<LeptonType<<std::endl;
  std::cout<<"                   "<<std::endl;

  double JetPtWboostedMin          = gConfigParser -> readDoubleOption("Input::JetPtWboostedMin");
  std::vector<double> JetPtCutMin  = gConfigParser -> readDoubleListOption("Input::JetPtCutMin");
  double JetEtaCutMax              = gConfigParser -> readDoubleOption("Input::JetEtaCutMax");
  double CleaningTreshold          = gConfigParser -> readDoubleOption("Input::CleaningTreshold");

  std::cout<<"                   "<<std::endl;
  std::cout<<" Input JetPtWboostedMin   "<<JetPtWboostedMin<<std::endl;
  std::cout<<" Input JetEtaCutMax       "<<JetEtaCutMax<<std::endl;
  std::cout<<" Input CleaningTreshold   "<<CleaningTreshold<<std::endl;
  //  std::cout<<" JetPtCutMin              "<<JetPtCutMin.size()<<std::endl;
  std::cout<<"                          "<<std::endl;

  int JetCollectionDimension  = gConfigParser -> readIntOption("Input::JetCollectionDimension");
  int NumJetMin               = gConfigParser -> readIntOption("Input::NumJetMin");

  std::cout<<"                   "<<std::endl;
  std::cout<<" JetCollectionDimension  "<<JetCollectionDimension<<std::endl;
  std::cout<<" NumJetMin   "<<NumJetMin<<std::endl;
  std::cout<<"                   "<<std::endl;

  std::string OutputRootDirectory   = gConfigParser -> readStringOption("Output::OutputRootDirectory");
  std::string OutputRootFile        = gConfigParser -> readStringOption("Output::OutputRootFile");

  std::cout<<"                   "<<std::endl;
  std::cout<<" Output OutputRootDirectory "<<OutputRootDirectory<<std::endl;
  std::cout<<" Output OutputRootFile "<<OutputRootFile<<std::endl;
  std::cout<<"                   "<<std::endl;

  std::string command = "if [ ! -e "+OutputRootDirectory+" ] ; then mkdir "+OutputRootDirectory+" ; fi";
  system(command.c_str());
  
  // create and open the input file 

  TFile * inputFile = new TFile((InputDirectory+"/"+InputRootFile).c_str(),"READ");

  // Efficiency Histo
  TH1F* SelectionEvents = new TH1F ("SelectionEvents","Selection Events",100,0,100);
  TH1F* SelectionEfficiency = new TH1F ("SelectionEfficiency","Selection Efficiency",100,0,100);
  int nstepEvents [100];
  for ( int iEvents = 0 ; iEvents < 100 ; iEvents ++) nstepEvents[iEvents] = 0 ;
  int nStep ;

  // Muon Sample Processing 
  if(LeptonType == "Muon"){

   std::cout<<" Enter in the Muon Category "<<std::endl;
   std::cout<<"                            "<<std::endl;

   // create and open the output file 

   TFile *outputFile = new TFile((OutputRootDirectory+"/"+OutputRootFile).c_str(),"RECREATE");
   outputFile->cd();

   std::cout<<" Open Input File : "<<inputFile->GetName()<<" TreeName  "<<TreeName<<std::endl;

   VBFMuonClass* MuonTree = new VBFMuonClass (inputFile,TreeName);
   
   MuonTree->SetReader(MuonTree->fTree);

   std::cout<<"                                   "<<std::endl;
   std::cout<<" Clone Tree  "<<std::endl;
   
   TTree *newtree = MuonTree->fTree->CloneTree(0);

   // Add new Branches 
   VBFMuonClass* NewMuonTree = new VBFMuonClass(newtree);
   
   NewMuonTree->SetNewBranches(NewMuonTree->fTree);
  
   // Loop on the events 

   std::cout<<"                                   "<<std::endl;
   std::cout << "Input Tree Number of Entries : " <<  MuonTree->fTree->GetEntries ()  << std::endl ;
   std::cout<<"                                   "<<std::endl;

   for(int iEntry = 0 ; iEntry <  MuonTree->fTree->GetEntries ()  ; iEntry++){

    nStep = 1 ;
    if(std::string(SelectionEvents->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEvents->GetXaxis()->SetBinLabel(nStep,"All Events");
    if(std::string(SelectionEfficiency->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEfficiency->GetXaxis()->SetBinLabel(nStep,"All Events");
    nstepEvents[nStep-1]++;
    nStep = 2;
    
    MuonTree->fTree->GetEntry(iEntry); 

    NewMuonTree->InitializateVariables(); 

    if (iEntry % 10000 == 0) std::cout << "reading event " << iEntry << std::endl ;

    // Basic Selections for boosted region 
    if(!(MuonTree->fReader->getInt("isgengdboostedWevt")[0]) || (MuonTree->fReader->getFloat("GroomedJet_CA8_deltaR_lca8jet")[0]) < TMath::Pi()/ 2.0  || 
        (MuonTree->fReader->getFloat("GroomedJet_CA8_pt")[0])< JetPtWboostedMin || (MuonTree->fReader->getFloat("W_pt")[0]) < JetPtWboostedMin ) continue ;
    
    if(std::string(SelectionEvents->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEvents->GetXaxis()->SetBinLabel(nStep,"Base Boosted W");    
    if(std::string(SelectionEfficiency->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEfficiency->GetXaxis()->SetBinLabel(nStep,"Base Boosted W");
    nstepEvents[nStep-1]++;
    nStep = 3;

    // Join forward and central PFJetCor Collection    
    std::vector<std::vector<TLorentzVector> > GroomedJet_CA8_Collection (JetPtCutMin.size());

    std::vector< std::vector <JetAK5> > JetPFCor_AK5_Collection (JetPtCutMin.size()) ;
    std::vector< std::vector <JetAK5> > CleanedJetPFCor_AK5_Collection(JetPtCutMin.size()) ;
    std::vector< std::vector <JetAK5> > HadronicW_AK5_Collection(JetPtCutMin.size()) ;

    std::vector<std::vector<TLorentzVector> > GenGroomedJet_CA8_Collection(JetPtCutMin.size()); 
    std::vector< std::vector <JetAK5> > GenJetPFCor_AK5_Collection(JetPtCutMin.size()); 
    std::vector< std::vector <JetAK5> > GenCleanedJetPFCor_AK5_Collection(JetPtCutMin.size());
    std::vector< std::vector <JetAK5> > GenHadronicW_AK5_Collection(JetPtCutMin.size());
   
    for(int iJetPtCutMin = 0 ; iJetPtCutMin < int(JetPtCutMin.size()); iJetPtCutMin++){

     for(int iJet = 0 ; iJet < JetCollectionDimension ; iJet++) { // run over the whole CA8 jet collection and took the 4vetcor 
      TLorentzVector JetTemp , GenJetTemp ;
      std::string nameCollection ; 
      JetTemp.SetPtEtaPhiE(MuonTree->fReader->getFloat("GroomedJet_CA8_pt")[iJet],MuonTree->fReader->getFloat("GroomedJet_CA8_eta")[iJet], 
                           MuonTree->fReader->getFloat("GroomedJet_CA8_phi")[iJet],MuonTree->fReader->getFloat("GroomedJet_CA8_e")[iJet]);

      if(NewMuonTree->fTree->FindBranch("GenGroomedJet_CA8_pt") && NewMuonTree->fTree->FindBranch("JetGen_Pt")){
       GenJetTemp.SetPtEtaPhiE(MuonTree->fReader->getFloat("GenGroomedJet_CA8_pt")[iJet],MuonTree->fReader->getFloat("GenGroomedJet_CA8_eta")[iJet], 
                              MuonTree->fReader->getFloat("GenGroomedJet_CA8_phi")[iJet],MuonTree->fReader->getFloat("GenGroomedJet_CA8_e")[iJet]);
      }

     // Selection on CA8 Jets -> pt cut on each jet over the threshold and acceptance
      if(fabs(JetTemp.Eta())<JetEtaCutMax && JetTemp.Pt()>JetPtCutMin.at(iJetPtCutMin))
		GroomedJet_CA8_Collection.at(iJetPtCutMin).push_back(JetTemp);
               

      if(NewMuonTree->fTree->FindBranch("GenGroomedJet_CA8_pt") && NewMuonTree->fTree->FindBranch("JetGen_Pt")){
       if(fabs(GenJetTemp.Eta())<JetEtaCutMax && GenJetTemp.Pt()>JetPtCutMin.at(iJetPtCutMin))
		GenGroomedJet_CA8_Collection.at(iJetPtCutMin).push_back(GenJetTemp);
      }
      // take the central AK5 jets 
      JetTemp.SetPtEtaPhiE(MuonTree->fReader->getFloat("JetPFCor_Pt")[iJet],MuonTree->fReader->getFloat("JetPFCor_Eta")[iJet],
 			  MuonTree->fReader->getFloat("JetPFCor_Phi")[iJet],MuonTree->fReader->getFloat("JetPFCor_E")[iJet]);

      if(NewMuonTree->fTree->FindBranch("GenGroomedJet_CA8_pt") && NewMuonTree->fTree->FindBranch("JetGen_Pt"))
       GenJetTemp.SetPtEtaPhiE(MuonTree->fReader->getFloat("JetGen_Pt")[iJet],MuonTree->fReader->getFloat("JetGen_Eta")[iJet],
   	   		       MuonTree->fReader->getFloat("JetGen_Phi")[iJet],MuonTree->fReader->getFloat("JetGen_E")[iJet]);
      
     // Selection on PF Cor Central Jets --> AK5
     if(fabs(JetTemp.Eta())<JetEtaCutMax && JetTemp.Pt()>JetPtCutMin.at(iJetPtCutMin)){
       JetAK5 tempJetAK5 (iJet,"JetPFCor",JetTemp);
       JetPFCor_AK5_Collection.at(iJetPtCutMin).push_back(tempJetAK5);
     }

     if(NewMuonTree->fTree->FindBranch("GenGroomedJet_CA8_pt") && NewMuonTree->fTree->FindBranch("JetGen_Pt")){
      if(fabs(GenJetTemp.Eta())<JetEtaCutMax && GenJetTemp.Pt()>JetPtCutMin.at(iJetPtCutMin)){
       JetAK5 tempJetAK5 (iJet,"GenJet",GenJetTemp);
       GenJetPFCor_AK5_Collection.at(iJetPtCutMin).push_back(tempJetAK5);
      }
      }
    }
   
    for(int iJet = 0 ; iJet < JetCollectionDimension ; iJet++) { //only AK5 forward jet over the pt threshold 
     TLorentzVector JetTemp ; 
     JetTemp.SetPtEtaPhiE(MuonTree->fReader->getFloat("JetPFCorVBFTag_Pt")[iJet],MuonTree->fReader->getFloat("JetPFCorVBFTag_Eta")[iJet], 
                          MuonTree->fReader->getFloat("JetPFCorVBFTag_Phi")[iJet],MuonTree->fReader->getFloat("JetPFCorVBFTag_E")[iJet]);

     // Selection on PF Cor Forward Jets --> AK5
     if(fabs(JetTemp.Eta())<JetEtaCutMax && JetTemp.Pt()>JetPtCutMin.at(iJetPtCutMin)){
       JetAK5 tempJetAK5 (iJet,"JetPFCorVBFTag",JetTemp);
       JetPFCor_AK5_Collection.at(iJetPtCutMin).push_back(tempJetAK5);
     }
    }
   
    //choose the jet corresponding to the hadronic W and fill new branches with its variables --> apply the pT cut on the whole CA8 jet and select the one with mass closer
    // to the W mass just to have another solution 
    float difference = 1000.;
    int iWHadronic = 0;

    for(size_t i = 0; i < GroomedJet_CA8_Collection.size() ; i ++){
      if ( MuonTree->fReader->getFloat("GroomedJet_CA8_pt")[i]>JetPtWboostedMin ) {
	if ( fabs (MuonTree->fReader->getFloat("GroomedJet_CA8_mass_pr")[i] - Wmass) < difference ) {
	  difference = fabs (MuonTree->fReader->getFloat("GroomedJet_CA8_mass_pr")[i] - Wmass);
	  iWHadronic = i;
	}
	}
    }
    if(iJetPtCutMin == 0){

     NewMuonTree -> WHadposition = iWHadronic;   //position of the hadronic W in the CA8Jet collection with pt over threshold + mass closer to the W one

     NewMuonTree->Hadronic_W_Jet_mass_uncorr    = MuonTree->fReader->getFloat("GroomedJet_CA8_mass_uncorr")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_mass_tr_uncorr = MuonTree->fReader->getFloat("GroomedJet_CA8_mass_tr_uncorr")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_mass_ft_uncorr = MuonTree->fReader->getFloat("GroomedJet_CA8_mass_ft_uncorr")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_mass_pr_uncorr = MuonTree->fReader->getFloat("GroomedJet_CA8_mass_pr_uncorr")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_massdrop_pr_uncorr = MuonTree->fReader->getFloat("GroomedJet_CA8_massdrop_pr_uncorr")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_tau2tau1           = MuonTree->fReader->getFloat("GroomedJet_CA8_tau2tau1")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_tau1 = MuonTree->fReader->getFloat("GroomedJet_CA8_tau1")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_tau2 = MuonTree->fReader->getFloat("GroomedJet_CA8_tau2")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_tau3 = MuonTree->fReader->getFloat("GroomedJet_CA8_tau3")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_tau4 = MuonTree->fReader->getFloat("GroomedJet_CA8_tau4")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_pt   = MuonTree->fReader->getFloat("GroomedJet_CA8_pt")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_eta  = MuonTree->fReader->getFloat("GroomedJet_CA8_eta")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_phi  = MuonTree->fReader->getFloat("GroomedJet_CA8_phi")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_e    = MuonTree->fReader->getFloat("GroomedJet_CA8_e")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_pt_tr_uncorr = MuonTree->fReader->getFloat("GroomedJet_CA8_pt_tr_uncorr")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_pt_tr  = MuonTree->fReader->getFloat("GroomedJet_CA8_pt_tr")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_eta_tr = MuonTree->fReader->getFloat("GroomedJet_CA8_eta_tr")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_phi_tr = MuonTree->fReader->getFloat("GroomedJet_CA8_phi_tr")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_e_tr   = MuonTree->fReader->getFloat("GroomedJet_CA8_e_tr")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_pt_ft_uncorr = MuonTree->fReader->getFloat("GroomedJet_CA8_pt_ft_uncorr")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_pt_ft  = MuonTree->fReader->getFloat("GroomedJet_CA8_pt_ft")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_eta_ft = MuonTree->fReader->getFloat("GroomedJet_CA8_eta_ft")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_phi_ft = MuonTree->fReader->getFloat("GroomedJet_CA8_phi_ft")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_e_ft   = MuonTree->fReader->getFloat("GroomedJet_CA8_e_ft")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_pt_pr_uncorr = MuonTree->fReader->getFloat("GroomedJet_CA8_pt_pr_uncorr")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_pt_pr  = MuonTree->fReader->getFloat("GroomedJet_CA8_pt_tr")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_eta_pr = MuonTree->fReader->getFloat("GroomedJet_CA8_eta_pr")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_phi_pr = MuonTree->fReader->getFloat("GroomedJet_CA8_phi_pr")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_e_pr   = MuonTree->fReader->getFloat("GroomedJet_CA8_e_pr")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_prsubjet1_px = MuonTree->fReader->getFloat("GroomedJet_CA8_prsubjet1_px")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_prsubjet1_py = MuonTree->fReader->getFloat("GroomedJet_CA8_prsubjet1_py")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_prsubjet1_pz = MuonTree->fReader->getFloat("GroomedJet_CA8_prsubjet1_pz")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_prsubjet1_e  = MuonTree->fReader->getFloat("GroomedJet_CA8_prsubjet1_e")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_prsubjet2_px = MuonTree->fReader->getFloat("GroomedJet_CA8_prsubjet2_px")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_prsubjet2_py = MuonTree->fReader->getFloat("GroomedJet_CA8_prsubjet2_py")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_prsubjet2_pz = MuonTree->fReader->getFloat("GroomedJet_CA8_prsubjet2_pz")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_prsubjet2_e  = MuonTree->fReader->getFloat("GroomedJet_CA8_prsubjet2_e")[iWHadronic];  
     NewMuonTree->Hadronic_W_Jet_mass    = MuonTree->fReader->getFloat("GroomedJet_CA8_mass")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_mass_tr = MuonTree->fReader->getFloat("GroomedJet_CA8_mass_tr")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_mass_ft = MuonTree->fReader->getFloat("GroomedJet_CA8_mass_ft")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_mass_pr = MuonTree->fReader->getFloat("GroomedJet_CA8_mass_pr")[iWHadronic];  
     NewMuonTree->Hadronic_W_Jet_massdrop = MuonTree->fReader->getFloat("GroomedJet_CA8_massdrop_pr")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_area    = MuonTree->fReader->getFloat("GroomedJet_CA8_area")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_area_tr = MuonTree->fReader->getFloat("GroomedJet_CA8_area_tr")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_area_ft = MuonTree->fReader->getFloat("GroomedJet_CA8_area_ft")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_area_pr = MuonTree->fReader->getFloat("GroomedJet_CA8_area_pr")[iWHadronic]; 
     NewMuonTree->Hadronic_W_Jet_jetconsituents = MuonTree->fReader->getFloat("GroomedJet_CA8_jetconstituents")[iWHadronic]; 
     NewMuonTree->Hadronic_W_Jet_jetcharge = MuonTree->fReader->getFloat("GroomedJet_CA8_jetcharge")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_rcores = MuonTree->fReader->getFloat("GroomedJet_CA8_rcores")[iWHadronic];  
     NewMuonTree->Hadronic_W_Jet_ptcores = MuonTree->fReader->getFloat("GroomedJet_CA8_ptcores")[iWHadronic];  
     NewMuonTree->Hadronic_W_Jet_planarflow = MuonTree->fReader->getFloat("GroomedJet_CA8_planarflow")[iWHadronic];  
     NewMuonTree->Hadronic_W_Jet_qjetmass = MuonTree->fReader->getFloat("GroomedJet_CA8_qjetmass")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_qjetmassdrop = MuonTree->fReader->getFloat("GroomedJet_CA8_qjetmassdrop")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_deltaR_ljet = MuonTree->fReader->getFloat("GroomedJet_CA8_deltaR_lca8jet")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_deltaphi_METjet = MuonTree->fReader->getFloat("GroomedJet_CA8_deltaphi_METca8jet")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_deltaphi_Vca8jet = MuonTree->fReader->getFloat("GroomedJet_CA8_deltaphi_Vca8jet")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_rcores01 = MuonTree->fReader->getFloat("GroomedJet_CA8_rcores01")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_rcores02 = MuonTree->fReader->getFloat("GroomedJet_CA8_rcores02")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_rcores03 = MuonTree->fReader->getFloat("GroomedJet_CA8_rcores03")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_rcores04 = MuonTree->fReader->getFloat("GroomedJet_CA8_rcores04")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_rcores05 = MuonTree->fReader->getFloat("GroomedJet_CA8_rcores05")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_rcores06 = MuonTree->fReader->getFloat("GroomedJet_CA8_rcores06")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_rcores07 = MuonTree->fReader->getFloat("GroomedJet_CA8_rcores07")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_rcores08 = MuonTree->fReader->getFloat("GroomedJet_CA8_rcores08")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_rcores09 = MuonTree->fReader->getFloat("GroomedJet_CA8_rcores09")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_rcores10 = MuonTree->fReader->getFloat("GroomedJet_CA8_rcores10")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_rcores11 = MuonTree->fReader->getFloat("GroomedJet_CA8_rcores11")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_ptcores01 = MuonTree->fReader->getFloat("GroomedJet_CA8_ptcores01")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_ptcores02 = MuonTree->fReader->getFloat("GroomedJet_CA8_ptcores02")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_ptcores03 = MuonTree->fReader->getFloat("GroomedJet_CA8_ptcores03")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_ptcores04 = MuonTree->fReader->getFloat("GroomedJet_CA8_ptcores04")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_ptcores05 = MuonTree->fReader->getFloat("GroomedJet_CA8_ptcores05")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_ptcores06 = MuonTree->fReader->getFloat("GroomedJet_CA8_ptcores06")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_ptcores07 = MuonTree->fReader->getFloat("GroomedJet_CA8_ptcores07")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_ptcores08 = MuonTree->fReader->getFloat("GroomedJet_CA8_ptcores08")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_ptcores09 = MuonTree->fReader->getFloat("GroomedJet_CA8_ptcores09")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_ptcores10 = MuonTree->fReader->getFloat("GroomedJet_CA8_ptcores10")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_ptcores11 = MuonTree->fReader->getFloat("GroomedJet_CA8_ptcores11")[iWHadronic]; 
     NewMuonTree->Hadronic_W_Jet_planarflow01 = MuonTree->fReader->getFloat("GroomedJet_CA8_planarflow01")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_planarflow02 = MuonTree->fReader->getFloat("GroomedJet_CA8_planarflow02")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_planarflow03 = MuonTree->fReader->getFloat("GroomedJet_CA8_planarflow03")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_planarflow04 = MuonTree->fReader->getFloat("GroomedJet_CA8_planarflow04")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_planarflow05 = MuonTree->fReader->getFloat("GroomedJet_CA8_planarflow05")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_planarflow06 = MuonTree->fReader->getFloat("GroomedJet_CA8_planarflow06")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_planarflow07 = MuonTree->fReader->getFloat("GroomedJet_CA8_planarflow07")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_planarflow08 = MuonTree->fReader->getFloat("GroomedJet_CA8_planarflow08")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_planarflow09 = MuonTree->fReader->getFloat("GroomedJet_CA8_planarflow09")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_planarflow10 = MuonTree->fReader->getFloat("GroomedJet_CA8_planarflow10")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_planarflow11 = MuonTree->fReader->getFloat("GroomedJet_CA8_planarflow11")[iWHadronic];
     NewMuonTree->Hadronic_W_Jet_mass_sensi_tr = MuonTree->fReader->getFloat("GroomedJet_CA8_mass_sensi_tr")[iWHadronic]; 
     NewMuonTree->Hadronic_W_Jet_mass_sensi_ft = MuonTree->fReader->getFloat("GroomedJet_CA8_mass_sensi_ft")[iWHadronic];     
     NewMuonTree->Hadronic_W_Jet_mass_sensi_pr = MuonTree->fReader->getFloat("GroomedJet_CA8_mass_sensi_pr")[iWHadronic]; 
     NewMuonTree->Hadronic_W_Jet_qjetmassvolatility = MuonTree->fReader->getFloat("GroomedJet_CA8_qjetmassvolatility")[iWHadronic]; 
     NewMuonTree->Hadronic_W_Jet_prsubjet1ptoverjetpt = MuonTree->fReader->getFloat("GroomedJet_CA8_prsubjet1ptoverjetpt")[iWHadronic]; 
     NewMuonTree->Hadronic_W_Jet_prsubjet2ptoverjetpt = MuonTree->fReader->getFloat("GroomedJet_CA8_prsubjet2ptoverjetpt")[iWHadronic]; 
     NewMuonTree->Hadronic_W_Jet_prsubjet1subjet2_deltaR = MuonTree->fReader->getFloat("GroomedJet_CA8_prsubjet1subjet2_deltaR")[iWHadronic]; 

     }
   }
    // Calculate Neutrino Pz using all the possible choices : type0 -> if real roots, pick the one nearest to the lepton Pz except when the Pz so chosen 
    //                                                                 is greater than 300 GeV in which case pick the most central root.               
    //                                                        type1 -> type = 1: if real roots, choose the one closest to the lepton Pz                                                                                                                                    if complex roots, use only the real part.     
    //                                                        type = 2: if real roots, choose the most central solution.                                                                                                                                          if complex roots, use only the real part.                                                                                                                                       type = 3: if real roots, pick the largest value of the cosine*                         

    TLorentzVector W_mu, W_Met;
   
    W_mu.SetPxPyPzE(MuonTree->fReader->getFloat("W_muon_px")[0],MuonTree->fReader->getFloat("W_muon_py")[0],
                    MuonTree->fReader->getFloat("W_muon_pz")[0],MuonTree->fReader->getFloat("W_muon_e")[0]);
    W_Met.SetPxPyPzE(MuonTree->fReader->getFloat("event_met_pfmet")[0] * TMath::Cos(MuonTree->fReader->getFloat("event_met_pfmetPhi")[0]),
                     MuonTree->fReader->getFloat("event_met_pfmet")[0] * TMath::Sin(MuonTree->fReader->getFloat("event_met_pfmetPhi")[0]),0.,
                     fabs(MuonTree->fReader->getFloat("event_met_pfmet")[0]));

    if(W_mu.Pt()<=0 || W_Met.Pt() <= 0 ){ std::cerr<<" Negative Lepton - Neutrino Pt "<<std::endl; continue ; }

    if(std::string(SelectionEvents->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEvents->GetXaxis()->SetBinLabel(nStep,"Leptonic 4V");
    if(std::string(SelectionEfficiency->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEfficiency->GetXaxis()->SetBinLabel(nStep,"Leptonic 4V");
    nstepEvents[nStep-1]++;
    nStep = 4;

    // type0 calculation of neutrino pZ
    METzCalculator<TLorentzVector> NeutrinoPz_type0;
    NeutrinoPz_type0.SetMET(W_Met);
    NeutrinoPz_type0.SetLepton(W_mu);
    NeutrinoPz_type0.SetLeptonType("muon");
    double pz1_type0 = NeutrinoPz_type0.Calculate(); // Default one -> according to type0
    double pz2_type0 = NeutrinoPz_type0.getOther(); // Default one

    // don't touch the neutrino pT
    TLorentzVector W_neutrino_type0_met; W_neutrino_type0_met.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type0,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type0*pz1_type0));
    NewMuonTree->W_mass_type0_met   = (W_neutrino_type0_met+W_mu).M();  
    NewMuonTree->W_pz_type0_met     = (W_neutrino_type0_met+W_mu).Pz();   
    NewMuonTree->W_nu1_pz_type0_met = pz1_type0; 
    NewMuonTree->W_nu2_pz_type0_met = pz2_type0;

    // chenge the neutrino pT in case of complex solution in order to make it real
    TLorentzVector W_neutrino_type0;  W_neutrino_type0.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type0,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type0*pz1_type0));

    if (NeutrinoPz_type0.IsComplex()) {// if this is a complix, change MET                                                                                                                            
     double nu_pt1 = NeutrinoPz_type0.getPtneutrino(1);
     double nu_pt2 = NeutrinoPz_type0.getPtneutrino(2);
   
     TLorentzVector W_neutrino_1;
     W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(MuonTree->fReader->getFloat("event_met_pfmetPhi")[0]), 
                             nu_pt1 * TMath::Sin(MuonTree->fReader->getFloat("event_met_pfmetPhi")[0]), pz1_type0, sqrt(nu_pt1*nu_pt1 + pz1_type0*pz1_type0) );
     TLorentzVector W_neutrino_2;
     W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(MuonTree->fReader->getFloat("event_met_pfmetPhi")[0]),
                             nu_pt2 * TMath::Sin(MuonTree->fReader->getFloat("event_met_pfmetPhi")[0]), pz1_type0, sqrt(nu_pt2*nu_pt2 + pz1_type0*pz1_type0) );

     if ( fabs((W_mu+W_neutrino_1).M()-Wmass) < fabs((W_mu+W_neutrino_2).M()-Wmass) )  W_neutrino_type0 = W_neutrino_1;
     else W_neutrino_type0 = W_neutrino_2;

    }

    NewMuonTree->W_mass_type0 = (W_mu+W_neutrino_type0).M();  
    NewMuonTree->W_pz_type0   = (W_mu+W_neutrino_type0).Pz();  
    NewMuonTree->W_nu1_pz_type0 = pz1_type0;  
    NewMuonTree->W_nu2_pz_type0 = pz2_type0;

  

    // type2 calculation of neutrino pZ
    METzCalculator<TLorentzVector> NeutrinoPz_type2;
    NeutrinoPz_type2.SetMET(W_Met);
    NeutrinoPz_type2.SetLepton(W_mu);
    NeutrinoPz_type2.SetLeptonType("muon");
    double pz1_type2 = NeutrinoPz_type2.Calculate(2); // Default one -> according to type2
    double pz2_type2 = NeutrinoPz_type2.getOther(); // Default one

    // don't touch the neutrino pT
    TLorentzVector W_neutrino_type2_met; W_neutrino_type2_met.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type2,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type2*pz1_type2));
    NewMuonTree->W_mass_type2_met   = (W_neutrino_type2_met+W_mu).M();  
    NewMuonTree->W_pz_type2_met     = (W_neutrino_type2_met+W_mu).Pz();   
    NewMuonTree->W_nu1_pz_type2_met = pz1_type2; 
    NewMuonTree->W_nu2_pz_type2_met = pz2_type2;

    // chenge the neutrino pT in case of complex solution in order to make it real
    TLorentzVector W_neutrino_type2;  W_neutrino_type2.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type2,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type2*pz1_type2));

    if (NeutrinoPz_type2.IsComplex()) {// if this is a complix, change MET                                                                                                                            
     double nu_pt1 = NeutrinoPz_type2.getPtneutrino(1);
     double nu_pt2 = NeutrinoPz_type2.getPtneutrino(2);
   
     TLorentzVector W_neutrino_1;
     W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(MuonTree->fReader->getFloat("event_met_pfmetPhi")[0]), 
                             nu_pt1 * TMath::Sin(MuonTree->fReader->getFloat("event_met_pfmetPhi")[0]), pz1_type2, sqrt(nu_pt1*nu_pt1 + pz1_type2*pz1_type2) );
     TLorentzVector W_neutrino_2;
     W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(MuonTree->fReader->getFloat("event_met_pfmetPhi")[0]),
                             nu_pt2 * TMath::Sin(MuonTree->fReader->getFloat("event_met_pfmetPhi")[0]), pz1_type2, sqrt(nu_pt2*nu_pt2 + pz1_type2*pz1_type2) );

     if ( fabs((W_mu+W_neutrino_1).M()-Wmass) < fabs((W_mu+W_neutrino_2).M()-Wmass) )  W_neutrino_type2 = W_neutrino_1;
     else W_neutrino_type2 = W_neutrino_2;

    }

    NewMuonTree->W_mass_type2 = (W_mu+W_neutrino_type2).M();  
    NewMuonTree->W_pz_type2   = (W_mu+W_neutrino_type2).Pz();  
    NewMuonTree->W_nu1_pz_type2 = pz1_type2;  
    NewMuonTree->W_nu2_pz_type2 = pz2_type2;


    //////////////////////////////
    
    TLorentzVector W_subjet1, W_subjet2 ;  // take the two subjet of the hardest CA8 after pruning
   
    W_subjet1.SetPxPyPzE(MuonTree->fReader->getFloat("GroomedJet_CA8_prsubjet1_px")[0],MuonTree->fReader->getFloat("GroomedJet_CA8_prsubjet1_py")[0],
         		 MuonTree->fReader->getFloat("GroomedJet_CA8_prsubjet1_pz")[0],MuonTree->fReader->getFloat("GroomedJet_CA8_prsubjet1_e")[0] );
    W_subjet2.SetPxPyPzE(MuonTree->fReader->getFloat("GroomedJet_CA8_prsubjet2_px")[0],MuonTree->fReader->getFloat("GroomedJet_CA8_prsubjet2_py")[0],
  		  	 MuonTree->fReader->getFloat("GroomedJet_CA8_prsubjet2_pz")[0],MuonTree->fReader->getFloat("GroomedJet_CA8_prsubjet2_e")[0] );

    if(W_subjet1.Pt() <= 0 || W_subjet2.Pt() <= 0){ std::cerr<<" Problem with subjets "<<std::endl; continue ; }

    if(std::string(SelectionEvents->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEvents->GetXaxis()->SetBinLabel(nStep,"Subjets 4V");
    if(std::string(SelectionEfficiency->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEfficiency->GetXaxis()->SetBinLabel(nStep,"Subjets 4V");
    nstepEvents[nStep-1]++;
    nStep = 5;


    TLorentzVector W_GroomedJet_CA8; 
    W_GroomedJet_CA8.SetPtEtaPhiE(MuonTree->fReader->getFloat("GroomedJet_CA8_pt")[0], MuonTree->fReader->getFloat("GroomedJet_CA8_eta")[0],
                                     MuonTree->fReader->getFloat("GroomedJet_CA8_phi")[0], MuonTree->fReader->getFloat("GroomedJet_CA8_e")[0]);

    if(W_GroomedJet_CA8.Pt() <=0){ std::cerr<<" Problem with pruned CA8 "<<std::endl; continue ;}

    if(std::string(SelectionEvents->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEvents->GetXaxis()->SetBinLabel(nStep,"Pruned CA8 4V");
    if(std::string(SelectionEfficiency->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEfficiency->GetXaxis()->SetBinLabel(nStep,"Pruned CA8 4V");
    nstepEvents[nStep-1]++;
    nStep = 6;

   
    // Kinematic Fit  for each neutrino type --> just type 0  and type 2                                                                                                                
    TLorentzVector fit_muon_type0(0,0,0,0), fit_neutrino_type0(0,0,0,0), fit_W_subjet1_type0(0,0,0,0), fit_W_subjet2_type0(0,0,0,0) ;
    TLorentzVector fit_muon_type2(0,0,0,0), fit_neutrino_type2(0,0,0,0), fit_W_subjet1_type2(0,0,0,0), fit_W_subjet2_type2(0,0,0,0) ;
    TLorentzVector fit_muon_type0_met(0,0,0,0), fit_neutrino_type0_met(0,0,0,0), fit_W_subjet1_type0_met(0,0,0,0), fit_W_subjet2_type0_met(0,0,0,0) ;
    TLorentzVector fit_muon_type2_met(0,0,0,0), fit_neutrino_type2_met(0,0,0,0), fit_W_subjet1_type2_met(0,0,0,0), fit_W_subjet2_type2_met(0,0,0,0) ;

    doKinematicFit(1, W_mu, W_neutrino_type0, W_subjet1, W_subjet2,  fit_muon_type0, fit_neutrino_type0, fit_W_subjet1_type0, fit_W_subjet2_type0, NewMuonTree->fit_chi2_type0, 
                   NewMuonTree->fit_NDF_type0, NewMuonTree->fit_status_type0, LeptonType);
    doKinematicFit(1, W_mu, W_neutrino_type2, W_subjet1, W_subjet2,  fit_muon_type2, fit_neutrino_type2, fit_W_subjet1_type2, fit_W_subjet2_type2, NewMuonTree->fit_chi2_type2, 
                   NewMuonTree->fit_NDF_type2, NewMuonTree->fit_status_type2, LeptonType);
    doKinematicFit(1, W_mu, W_neutrino_type0_met, W_subjet1, W_subjet2,  fit_muon_type0_met, fit_neutrino_type0_met, fit_W_subjet1_type0_met, fit_W_subjet2_type0_met, 
                   NewMuonTree->fit_chi2_type0_met, NewMuonTree->fit_NDF_type0_met, NewMuonTree->fit_status_type0_met, LeptonType);
    doKinematicFit(1, W_mu, W_neutrino_type2_met, W_subjet1, W_subjet2,  fit_muon_type2_met, fit_neutrino_type2_met, fit_W_subjet1_type2_met, fit_W_subjet2_type2_met,
                   NewMuonTree->fit_chi2_type2_met, NewMuonTree->fit_NDF_type2_met, NewMuonTree->fit_status_type2_met, LeptonType);
    
    if(fit_muon_type0.Pt() >0 && fit_neutrino_type0.Pt()>0 && fit_W_subjet1_type0.Pt()>0 && fit_W_subjet2_type0.Pt()>0){

     NewMuonTree->fit_mu_px_type0 = fit_muon_type0.Px();
     NewMuonTree->fit_mu_py_type0 = fit_muon_type0.Py(); 
     NewMuonTree->fit_mu_pz_type0 = fit_muon_type0.Pz(); 
     NewMuonTree->fit_mu_e_type0  = fit_muon_type0.E();
     NewMuonTree->fit_nv_px_type0 = fit_neutrino_type0.Px(); 
     NewMuonTree->fit_nv_py_type0 = fit_neutrino_type0.Py(); 
     NewMuonTree->fit_nv_pz_type0 = fit_neutrino_type0.Pz(); 
     NewMuonTree->fit_nv_e_type0  = fit_neutrino_type0.E();

     NewMuonTree->fit_subjet1_px_type0 = fit_W_subjet1_type0.Px();  NewMuonTree->fit_subjet2_px_type0 = fit_W_subjet2_type0.Px();
     NewMuonTree->fit_subjet1_py_type0 = fit_W_subjet1_type0.Py();  NewMuonTree->fit_subjet2_py_type0 = fit_W_subjet2_type0.Py();  
     NewMuonTree->fit_subjet1_pz_type0 = fit_W_subjet1_type0.Pz();  NewMuonTree->fit_subjet2_pz_type0 = fit_W_subjet2_type0.Pz();
     NewMuonTree->fit_subjet1_e_type0  = fit_W_subjet1_type0.E();   NewMuonTree->fit_subjet2_e_type0  = fit_W_subjet2_type0.E();
     NewMuonTree->fit_subjet1_m_type0  = fit_W_subjet1_type0.M();   NewMuonTree->fit_subjet2_m_type0  = fit_W_subjet2_type0.M();

     NewMuonTree->fit_lvj_m_type0   = (fit_muon_type0+fit_neutrino_type0+fit_W_subjet1_type0+fit_W_subjet2_type0).M();
     NewMuonTree->fit_lv_m_type0    = (fit_muon_type0+fit_neutrino_type0).M();
     NewMuonTree->fit_j_m_type0     = (fit_W_subjet1_type0+fit_W_subjet2_type0).M();
     NewMuonTree->fit_lvj_pt_type0  = (fit_muon_type0+fit_neutrino_type0+fit_W_subjet1_type0+fit_W_subjet2_type0).M();
     NewMuonTree->fit_lvj_eta_type0 = (fit_muon_type0+fit_neutrino_type0+fit_W_subjet1_type0+fit_W_subjet2_type0).Eta();
     NewMuonTree->fit_lvj_phi_type0 = (fit_muon_type0+fit_neutrino_type0+fit_W_subjet1_type0+fit_W_subjet2_type0).Phi();
     NewMuonTree->fit_lvj_e_type0   = (fit_muon_type0+fit_neutrino_type0+fit_W_subjet1_type0+fit_W_subjet2_type0).E();
    }

    if(fit_muon_type2.Pt() >0 && fit_neutrino_type2.Pt()>0 && fit_W_subjet1_type2.Pt()>0 && fit_W_subjet2_type2.Pt()>0){

     NewMuonTree->fit_mu_px_type2 = fit_muon_type2.Px();
     NewMuonTree->fit_mu_py_type2 = fit_muon_type2.Py(); 
     NewMuonTree->fit_mu_pz_type2 = fit_muon_type2.Pz(); 
     NewMuonTree->fit_mu_e_type2  = fit_muon_type2.E();
     NewMuonTree->fit_nv_px_type2 = fit_neutrino_type2.Px(); 
     NewMuonTree->fit_nv_py_type2 = fit_neutrino_type2.Py(); 
     NewMuonTree->fit_nv_pz_type2 = fit_neutrino_type2.Pz(); 
     NewMuonTree->fit_nv_e_type2  = fit_neutrino_type2.E();

     NewMuonTree->fit_subjet1_px_type2 = fit_W_subjet1_type2.Px();  NewMuonTree->fit_subjet2_px_type2 = fit_W_subjet2_type2.Px();
     NewMuonTree->fit_subjet1_py_type2 = fit_W_subjet1_type2.Py();  NewMuonTree->fit_subjet2_py_type2 = fit_W_subjet2_type2.Py();  
     NewMuonTree->fit_subjet1_pz_type2 = fit_W_subjet1_type2.Pz();  NewMuonTree->fit_subjet2_pz_type2 = fit_W_subjet2_type2.Pz();
     NewMuonTree->fit_subjet1_e_type2  = fit_W_subjet1_type2.E();   NewMuonTree->fit_subjet2_e_type2  = fit_W_subjet2_type2.E();
     NewMuonTree->fit_subjet1_m_type2  = fit_W_subjet1_type2.M();   NewMuonTree->fit_subjet2_m_type2  = fit_W_subjet2_type2.M();

     NewMuonTree->fit_lvj_m_type2   = (fit_muon_type2+fit_neutrino_type2+fit_W_subjet1_type2+fit_W_subjet2_type2).M();
     NewMuonTree->fit_lv_m_type2    = (fit_muon_type2+fit_neutrino_type2).M();
     NewMuonTree->fit_j_m_type2     = (fit_W_subjet1_type2+fit_W_subjet2_type2).M();
     NewMuonTree->fit_lvj_pt_type2  = (fit_muon_type2+fit_neutrino_type2+fit_W_subjet1_type2+fit_W_subjet2_type2).M();
     NewMuonTree->fit_lvj_eta_type2 = (fit_muon_type2+fit_neutrino_type2+fit_W_subjet1_type2+fit_W_subjet2_type2).Eta();
     NewMuonTree->fit_lvj_phi_type2 = (fit_muon_type2+fit_neutrino_type2+fit_W_subjet1_type2+fit_W_subjet2_type2).Phi();
     NewMuonTree->fit_lvj_e_type2   = (fit_muon_type2+fit_neutrino_type2+fit_W_subjet1_type2+fit_W_subjet2_type2).E();
    }

    if(fit_muon_type0_met.Pt() >0 && fit_neutrino_type0_met.Pt()>0 && fit_W_subjet1_type0_met.Pt()>0 && fit_W_subjet2_type0_met.Pt()>0){

     NewMuonTree->fit_mu_px_type0_met = fit_muon_type0_met.Px();
     NewMuonTree->fit_mu_py_type0_met = fit_muon_type0_met.Py(); 
     NewMuonTree->fit_mu_pz_type0_met = fit_muon_type0_met.Pz(); 
     NewMuonTree->fit_mu_e_type0_met  = fit_muon_type0_met.E();
     NewMuonTree->fit_nv_px_type0_met = fit_neutrino_type0_met.Px(); 
     NewMuonTree->fit_nv_py_type0_met = fit_neutrino_type0_met.Py(); 
     NewMuonTree->fit_nv_pz_type0_met = fit_neutrino_type0_met.Pz(); 
     NewMuonTree->fit_nv_e_type0_met  = fit_neutrino_type0_met.E();

     NewMuonTree->fit_subjet1_px_type0_met = fit_W_subjet1_type0_met.Px();  NewMuonTree->fit_subjet2_px_type0_met = fit_W_subjet2_type0_met.Px();
     NewMuonTree->fit_subjet1_py_type0_met = fit_W_subjet1_type0_met.Py();  NewMuonTree->fit_subjet2_py_type0_met = fit_W_subjet2_type0_met.Py();  
     NewMuonTree->fit_subjet1_pz_type0_met = fit_W_subjet1_type0_met.Pz();  NewMuonTree->fit_subjet2_pz_type0_met = fit_W_subjet2_type0_met.Pz();
     NewMuonTree->fit_subjet1_e_type0_met  = fit_W_subjet1_type0_met.E();   NewMuonTree->fit_subjet2_e_type0_met  = fit_W_subjet2_type0_met.E();
     NewMuonTree->fit_subjet1_m_type0_met  = fit_W_subjet1_type0_met.M();   NewMuonTree->fit_subjet2_m_type0_met  = fit_W_subjet2_type0_met.M();

     NewMuonTree->fit_lvj_m_type0_met   = (fit_muon_type0_met+fit_neutrino_type0_met+fit_W_subjet1_type0_met+fit_W_subjet2_type0_met).M();
     NewMuonTree->fit_lv_m_type0_met    = (fit_muon_type0_met+fit_neutrino_type0_met).M();
     NewMuonTree->fit_j_m_type0_met     = (fit_W_subjet1_type0_met+fit_W_subjet2_type0_met).M();
     NewMuonTree->fit_lvj_pt_type0_met  = (fit_muon_type0_met+fit_neutrino_type0_met+fit_W_subjet1_type0_met+fit_W_subjet2_type0_met).M();
     NewMuonTree->fit_lvj_eta_type0_met = (fit_muon_type0_met+fit_neutrino_type0_met+fit_W_subjet1_type0_met+fit_W_subjet2_type0_met).Eta();
     NewMuonTree->fit_lvj_phi_type0_met = (fit_muon_type0_met+fit_neutrino_type0_met+fit_W_subjet1_type0_met+fit_W_subjet2_type0_met).Phi();
     NewMuonTree->fit_lvj_e_type0_met   = (fit_muon_type0_met+fit_neutrino_type0_met+fit_W_subjet1_type0_met+fit_W_subjet2_type0_met).E();
    }


    if(fit_muon_type2_met.Pt() >0 && fit_neutrino_type2_met.Pt()>0 && fit_W_subjet1_type2_met.Pt()>0 && fit_W_subjet2_type2_met.Pt()>0){

     NewMuonTree->fit_mu_px_type2_met = fit_muon_type2_met.Px();
     NewMuonTree->fit_mu_py_type2_met = fit_muon_type2_met.Py(); 
     NewMuonTree->fit_mu_pz_type2_met = fit_muon_type2_met.Pz(); 
     NewMuonTree->fit_mu_e_type2_met  = fit_muon_type2_met.E();
     NewMuonTree->fit_nv_px_type2_met = fit_neutrino_type2_met.Px(); 
     NewMuonTree->fit_nv_py_type2_met = fit_neutrino_type2_met.Py(); 
     NewMuonTree->fit_nv_pz_type2_met = fit_neutrino_type2_met.Pz(); 
     NewMuonTree->fit_nv_e_type2_met  = fit_neutrino_type2_met.E();

     NewMuonTree->fit_subjet1_px_type2_met = fit_W_subjet1_type2_met.Px();  NewMuonTree->fit_subjet2_px_type2_met = fit_W_subjet2_type2_met.Px();
     NewMuonTree->fit_subjet1_py_type2_met = fit_W_subjet1_type2_met.Py();  NewMuonTree->fit_subjet2_py_type2_met = fit_W_subjet2_type2_met.Py();  
     NewMuonTree->fit_subjet1_pz_type2_met = fit_W_subjet1_type2_met.Pz();  NewMuonTree->fit_subjet2_pz_type2_met = fit_W_subjet2_type2_met.Pz();
     NewMuonTree->fit_subjet1_e_type2_met  = fit_W_subjet1_type2_met.E();   NewMuonTree->fit_subjet2_e_type2_met  = fit_W_subjet2_type2_met.E();
     NewMuonTree->fit_subjet1_m_type2_met  = fit_W_subjet1_type2_met.M();   NewMuonTree->fit_subjet2_m_type2_met  = fit_W_subjet2_type2_met.M();

     NewMuonTree->fit_lvj_m_type2_met   = (fit_muon_type2_met+fit_neutrino_type2_met+fit_W_subjet1_type2_met+fit_W_subjet2_type2_met).M();
     NewMuonTree->fit_lv_m_type2_met    = (fit_muon_type2_met+fit_neutrino_type2_met).M();
     NewMuonTree->fit_j_m_type2_met     = (fit_W_subjet1_type2_met+fit_W_subjet2_type2_met).M();
     NewMuonTree->fit_lvj_pt_type2_met  = (fit_muon_type2_met+fit_neutrino_type2_met+fit_W_subjet1_type2_met+fit_W_subjet2_type2_met).M();
     NewMuonTree->fit_lvj_eta_type2_met = (fit_muon_type2_met+fit_neutrino_type2_met+fit_W_subjet1_type2_met+fit_W_subjet2_type2_met).Eta();
     NewMuonTree->fit_lvj_phi_type2_met = (fit_muon_type2_met+fit_neutrino_type2_met+fit_W_subjet1_type2_met+fit_W_subjet2_type2_met).Phi();
     NewMuonTree->fit_lvj_e_type2_met   = (fit_muon_type2_met+fit_neutrino_type2_met+fit_W_subjet1_type2_met+fit_W_subjet2_type2_met).E();
    }

    
    NewMuonTree->boosted_lvj_m_type0   = (W_mu+W_neutrino_type0+W_subjet1+W_subjet2).M();
    NewMuonTree->boosted_lv_m_type0    = (W_mu+W_neutrino_type0).M();
    NewMuonTree->boosted_j_m_type0     = (W_subjet1+W_subjet2).M();
    NewMuonTree->boosted_lvj_pt_type0  = (W_mu+W_neutrino_type0+W_subjet1+W_subjet2).Pt();
    NewMuonTree->boosted_lvj_eta_type0 = (W_mu+W_neutrino_type0+W_subjet1+W_subjet2).Eta();
    NewMuonTree->boosted_lvj_phi_type0 = (W_mu+W_neutrino_type0+W_subjet1+W_subjet2).Phi();
    NewMuonTree->boosted_lvj_e_type0   = (W_mu+W_neutrino_type0+W_subjet1+W_subjet2).E();
 
    NewMuonTree->boostedW_lvj_m_type0   = (W_mu+W_neutrino_type0+W_GroomedJet_CA8).M();
    NewMuonTree->boostedW_lv_m_type0    = (W_mu+W_neutrino_type0).M();
    NewMuonTree->boostedW_j_m_type0     = MuonTree->fReader->getFloat("GroomedJet_CA8_mass_pr")[0];
    NewMuonTree->boostedW_lvj_pt_type0  = (W_mu+W_neutrino_type0+W_GroomedJet_CA8).Pt();
    NewMuonTree->boostedW_lvj_eta_type0 = (W_mu+W_neutrino_type0+W_GroomedJet_CA8).Eta();
    NewMuonTree->boostedW_lvj_phi_type0 = (W_mu+W_neutrino_type0+W_GroomedJet_CA8).Phi();
    NewMuonTree->boostedW_lvj_e_type0   = (W_mu+W_neutrino_type0+W_GroomedJet_CA8).E();

    NewMuonTree->boosted_lvj_m_type2   = (W_mu+W_neutrino_type2+W_subjet1+W_subjet2).M();
    NewMuonTree->boosted_lv_m_type2    = (W_mu+W_neutrino_type2).M();
    NewMuonTree->boosted_j_m_type2     = (W_subjet1+W_subjet2).M();
    NewMuonTree->boosted_lvj_pt_type2  = (W_mu+W_neutrino_type2+W_subjet1+W_subjet2).Pt();
    NewMuonTree->boosted_lvj_eta_type2 = (W_mu+W_neutrino_type2+W_subjet1+W_subjet2).Eta();
    NewMuonTree->boosted_lvj_phi_type2 = (W_mu+W_neutrino_type2+W_subjet1+W_subjet2).Phi();
    NewMuonTree->boosted_lvj_e_type2   = (W_mu+W_neutrino_type2+W_subjet1+W_subjet2).E();
 
    NewMuonTree->boostedW_lvj_m_type2   = (W_mu+W_neutrino_type2+W_GroomedJet_CA8).M();
    NewMuonTree->boostedW_lv_m_type2    = (W_mu+W_neutrino_type2).M();
    NewMuonTree->boostedW_j_m_type2     = MuonTree->fReader->getFloat("GroomedJet_CA8_mass_pr")[0];
    NewMuonTree->boostedW_lvj_pt_type2  = (W_mu+W_neutrino_type2+W_GroomedJet_CA8).Pt();
    NewMuonTree->boostedW_lvj_eta_type2 = (W_mu+W_neutrino_type2+W_GroomedJet_CA8).Eta();
    NewMuonTree->boostedW_lvj_phi_type2 = (W_mu+W_neutrino_type2+W_GroomedJet_CA8).Phi();
    NewMuonTree->boostedW_lvj_e_type2   = (W_mu+W_neutrino_type2+W_GroomedJet_CA8).E();

    NewMuonTree->boosted_lvj_m_type0_met   = (W_mu+W_neutrino_type0_met+W_subjet1+W_subjet2).M();
    NewMuonTree->boosted_lv_m_type0_met    = (W_mu+W_neutrino_type0_met).M();
    NewMuonTree->boosted_j_m_type0_met     = (W_subjet1+W_subjet2).M();
    NewMuonTree->boosted_lvj_pt_type0_met  = (W_mu+W_neutrino_type0_met+W_subjet1+W_subjet2).Pt();
    NewMuonTree->boosted_lvj_eta_type0_met = (W_mu+W_neutrino_type0_met+W_subjet1+W_subjet2).Eta();
    NewMuonTree->boosted_lvj_phi_type0_met = (W_mu+W_neutrino_type0_met+W_subjet1+W_subjet2).Phi();
    NewMuonTree->boosted_lvj_e_type0_met   = (W_mu+W_neutrino_type0_met+W_subjet1+W_subjet2).E();
 
    NewMuonTree->boostedW_lvj_m_type0_met   = (W_mu+W_neutrino_type0_met+W_GroomedJet_CA8).M();
    NewMuonTree->boostedW_lv_m_type0_met    = (W_mu+W_neutrino_type0_met).M();
    NewMuonTree->boostedW_j_m_type0_met     = MuonTree->fReader->getFloat("GroomedJet_CA8_mass_pr")[0];
    NewMuonTree->boostedW_lvj_pt_type0_met  = (W_mu+W_neutrino_type0_met+W_GroomedJet_CA8).Pt();
    NewMuonTree->boostedW_lvj_eta_type0_met = (W_mu+W_neutrino_type0_met+W_GroomedJet_CA8).Eta();
    NewMuonTree->boostedW_lvj_phi_type0_met = (W_mu+W_neutrino_type0_met+W_GroomedJet_CA8).Phi();
    NewMuonTree->boostedW_lvj_e_type0_met   = (W_mu+W_neutrino_type0_met+W_GroomedJet_CA8).E();

    NewMuonTree->boosted_lvj_m_type2_met   = (W_mu+W_neutrino_type2_met+W_subjet1+W_subjet2).M();
    NewMuonTree->boosted_lv_m_type2_met    = (W_mu+W_neutrino_type2_met).M();
    NewMuonTree->boosted_j_m_type2_met     = (W_subjet1+W_subjet2).M();
    NewMuonTree->boosted_lvj_pt_type2_met  = (W_mu+W_neutrino_type2_met+W_subjet1+W_subjet2).Pt();
    NewMuonTree->boosted_lvj_eta_type2_met = (W_mu+W_neutrino_type2_met+W_subjet1+W_subjet2).Eta();
    NewMuonTree->boosted_lvj_phi_type2_met = (W_mu+W_neutrino_type2_met+W_subjet1+W_subjet2).Phi();
    NewMuonTree->boosted_lvj_e_type2_met   = (W_mu+W_neutrino_type2_met+W_subjet1+W_subjet2).E();
 
    NewMuonTree->boostedW_lvj_m_type2_met   = (W_mu+W_neutrino_type2_met+W_GroomedJet_CA8).M();
    NewMuonTree->boostedW_lv_m_type2_met    = (W_mu+W_neutrino_type2_met).M();
    NewMuonTree->boostedW_j_m_type2_met     = MuonTree->fReader->getFloat("GroomedJet_CA8_mass_pr")[0];
    NewMuonTree->boostedW_lvj_pt_type2_met  = (W_mu+W_neutrino_type2_met+W_GroomedJet_CA8).Pt();
    NewMuonTree->boostedW_lvj_eta_type2_met = (W_mu+W_neutrino_type2_met+W_GroomedJet_CA8).Eta();
    NewMuonTree->boostedW_lvj_phi_type2_met = (W_mu+W_neutrino_type2_met+W_GroomedJet_CA8).Phi();
    
    
    // Angles for the central Higgs Kinematics
    double costheta1, costheta2, phi, costhetastar, phistar1, phistar2;

    //Use the Subjet in the Boosted W Analyisis                                                                                                                                             
    if (MuonTree->fReader->getFloat("W_muon_charge")[0] < 0) calculateAngles(W_mu, W_neutrino_type0,W_subjet1,W_subjet2,costheta1,costheta2,phi,costhetastar,phistar1,phistar2);
    else calculateAngles(W_neutrino_type0, W_mu, W_subjet1, W_subjet2, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);
   
    NewMuonTree->boosted_wjj_ang_ha_type0   = costheta1;
    NewMuonTree->boosted_wjj_ang_hb_type0   = fabs(costheta2); 
    NewMuonTree->boosted_wjj_ang_hs_type0   = costhetastar;
    NewMuonTree->boosted_wjj_ang_phi_type0  = phi;
    NewMuonTree->boosted_wjj_ang_phia_type0 = phistar1;																
    NewMuonTree->boosted_wjj_ang_phib_type0 = phistar2;

    if (MuonTree->fReader->getFloat("W_muon_charge")[0] < 0) calculateAngles(W_mu, W_neutrino_type2,W_subjet1,W_subjet2,costheta1,costheta2,phi,costhetastar,phistar1,phistar2);
    else calculateAngles(W_neutrino_type2, W_mu, W_subjet1, W_subjet2, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);
   
    NewMuonTree->boosted_wjj_ang_ha_type2   = costheta1;
    NewMuonTree->boosted_wjj_ang_hb_type2   = fabs(costheta2); 
    NewMuonTree->boosted_wjj_ang_hs_type2   = costhetastar;
    NewMuonTree->boosted_wjj_ang_phi_type2  = phi;
    NewMuonTree->boosted_wjj_ang_phia_type2 = phistar1;															     
    NewMuonTree->boosted_wjj_ang_phib_type2 = phistar2;

    if (MuonTree->fReader->getFloat("W_muon_charge")[0] < 0) calculateAngles(W_mu, W_neutrino_type0_met,W_subjet1,W_subjet2,costheta1,costheta2,phi,costhetastar,phistar1,phistar2);
    else calculateAngles(W_neutrino_type0_met, W_mu, W_subjet1, W_subjet2, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);
   
    NewMuonTree->boosted_wjj_ang_ha_type0_met   = costheta1;
    NewMuonTree->boosted_wjj_ang_hb_type0_met   = fabs(costheta2); 
    NewMuonTree->boosted_wjj_ang_hs_type0_met   = costhetastar;
    NewMuonTree->boosted_wjj_ang_phi_type0_met  = phi;
    NewMuonTree->boosted_wjj_ang_phia_type0_met = phistar1;															
    NewMuonTree->boosted_wjj_ang_phib_type0_met = phistar2;

    if (MuonTree->fReader->getFloat("W_muon_charge")[0] < 0) calculateAngles(W_mu, W_neutrino_type2_met,W_subjet1,W_subjet2,costheta1,costheta2,phi,costhetastar,phistar1,phistar2);
    else calculateAngles(W_neutrino_type2_met, W_mu, W_subjet1, W_subjet2, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);
   
    NewMuonTree->boosted_wjj_ang_ha_type2_met   = costheta1;
    NewMuonTree->boosted_wjj_ang_hb_type2_met   = fabs(costheta2); 
    NewMuonTree->boosted_wjj_ang_hs_type2_met   = costhetastar;
    NewMuonTree->boosted_wjj_ang_phi_type2_met  = phi;
    NewMuonTree->boosted_wjj_ang_phia_type2_met = phistar1;														
    NewMuonTree->boosted_wjj_ang_phib_type2_met = phistar2;

    // Clean AK5 Jet Collection from the hadronic W and sotre the jet binning
    std::vector<int> numberJetBin ;
    for( size_t iJetPtCutMin = 0; iJetPtCutMin < JetPtCutMin.size(); iJetPtCutMin++){
      for(size_t iJet = 0; iJet < JetPFCor_AK5_Collection.at(iJetPtCutMin).size() ; iJet ++){       
	if(deltaR(JetPFCor_AK5_Collection.at(iJetPtCutMin).at(iJet).Momentum_.Phi(),GroomedJet_CA8_Collection.at(iJetPtCutMin).at(0).Phi(),
		  JetPFCor_AK5_Collection.at(iJetPtCutMin).at(iJet).Momentum_.Eta(),GroomedJet_CA8_Collection.at(iJetPtCutMin).at(0).Eta()) < CleaningTreshold ){
	  HadronicW_AK5_Collection.at(iJetPtCutMin).push_back(JetPFCor_AK5_Collection.at(iJetPtCutMin).at(iJet)); continue ;}

	CleanedJetPFCor_AK5_Collection.at(iJetPtCutMin).push_back(JetPFCor_AK5_Collection.at(iJetPtCutMin).at(iJet));

      }
 
      if(!CleanedJetPFCor_AK5_Collection.at(iJetPtCutMin).empty())
       numberJetBin.push_back(CleanedJetPFCor_AK5_Collection.at(iJetPtCutMin).size());
      else
       numberJetBin.push_back(0);

    }

    std::vector<int> numberJetBinGen ;
    if(NewMuonTree->fTree->FindBranch("GenGroomedJet_CA8_pt") && NewMuonTree->fTree->FindBranch("JetGen_Pt")){
     for( size_t iJetPtCutMin = 0; iJetPtCutMin < JetPtCutMin.size(); iJetPtCutMin++){
      for(size_t iJet = 0; iJet < GenJetPFCor_AK5_Collection.at(iJetPtCutMin).size() ; iJet ++){
	if(GenGroomedJet_CA8_Collection.at(iJetPtCutMin).empty()) break;
	if(deltaR(GenJetPFCor_AK5_Collection.at(iJetPtCutMin).at(iJet).Momentum_.Phi(),GenGroomedJet_CA8_Collection.at(iJetPtCutMin).at(0).Phi(),
		  GenJetPFCor_AK5_Collection.at(iJetPtCutMin).at(iJet).Momentum_.Eta(),GenGroomedJet_CA8_Collection.at(iJetPtCutMin).at(0).Eta()) < CleaningTreshold ){
	  GenHadronicW_AK5_Collection.at(iJetPtCutMin).push_back(GenJetPFCor_AK5_Collection.at(iJetPtCutMin).at(iJet)); continue ;}

	GenCleanedJetPFCor_AK5_Collection.at(iJetPtCutMin).push_back(GenJetPFCor_AK5_Collection.at(iJetPtCutMin).at(iJet));
        
      }
      
      if(!GenCleanedJetPFCor_AK5_Collection.at(iJetPtCutMin).empty())
	numberJetBinGen.push_back(GenCleanedJetPFCor_AK5_Collection.at(iJetPtCutMin).size());
      else
       numberJetBinGen.push_back(0);
     
     }
    }

    NewMuonTree -> numberJetBin = numberJetBin;
    NewMuonTree -> numberJetBinGen = numberJetBinGen;

    if (NewMuonTree -> numberJetBin.at(0) == 1 && CleanedJetPFCor_AK5_Collection.at(0).size() == 1){

      std::sort(CleanedJetPFCor_AK5_Collection.at(0).begin(),CleanedJetPFCor_AK5_Collection.at(0).end(),TLVP_PtSort());
      std::sort(GenCleanedJetPFCor_AK5_Collection.at(0).begin(),GenCleanedJetPFCor_AK5_Collection.at(0).end(),TLVP_PtSort());

      NewMuonTree->vbf_maxpt_j1_e   = CleanedJetPFCor_AK5_Collection.at(0).at(0).Momentum_.E();
      NewMuonTree->vbf_maxpt_j1_pt  = CleanedJetPFCor_AK5_Collection.at(0).at(0).Momentum_.Pt();
      NewMuonTree->vbf_maxpt_j1_eta = CleanedJetPFCor_AK5_Collection.at(0).at(0).Momentum_.Eta();
      NewMuonTree->vbf_maxpt_j1_phi = CleanedJetPFCor_AK5_Collection.at(0).at(0).Momentum_.Phi();
      NewMuonTree->vbf_maxpt_j1_m   = CleanedJetPFCor_AK5_Collection.at(0).at(0).Momentum_.M(); 
      
      if(!numberJetBinGen.empty()){
	if(numberJetBinGen.at(0)==1 && GenCleanedJetPFCor_AK5_Collection.at(0).size()==1){
        NewMuonTree->vbf_maxpt_j1_e_gen   = GenCleanedJetPFCor_AK5_Collection.at(0).at(0).Momentum_.E();
        NewMuonTree->vbf_maxpt_j1_pt_gen  = GenCleanedJetPFCor_AK5_Collection.at(0).at(0).Momentum_.Pt();
        NewMuonTree->vbf_maxpt_j1_eta_gen = GenCleanedJetPFCor_AK5_Collection.at(0).at(0).Momentum_.Eta();
        NewMuonTree->vbf_maxpt_j1_phi_gen = GenCleanedJetPFCor_AK5_Collection.at(0).at(0).Momentum_.Phi();
        NewMuonTree->vbf_maxpt_j1_m_gen   = GenCleanedJetPFCor_AK5_Collection.at(0).at(0).Momentum_.M(); 
        NewMuonTree->vbf_maxpt_j1_bDiscriminatorSSVHE_gen = MuonTree->fReader->getFloat("JetGen_bDiscriminatorSSVHE")[GenCleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
        NewMuonTree->vbf_maxpt_j1_bDiscriminatorTCHE_gen  = MuonTree->fReader->getFloat("JetGen_bDiscriminatorTCHE")[GenCleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
        NewMuonTree->vbf_maxpt_j1_bDiscriminatorCSV_gen   = MuonTree->fReader->getFloat("JetGen_bDiscriminatorCSV")[GenCleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
        NewMuonTree->vbf_maxpt_j1_bDiscriminatorSSVHP_gen = MuonTree->fReader->getFloat("JetGen_bDiscriminatorSSVHP")[GenCleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
        NewMuonTree->vbf_maxpt_j1_bDiscriminatorTCHP_gen  = MuonTree->fReader->getFloat("JetGen_bDiscriminatorTCHP")[GenCleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
       }
      }
      
      NewMuonTree->vbf_maxpt_j1_QGLikelihood = MuonTree->fReader->getFloat("JetPFCor_QGLikelihood")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;

      NewMuonTree->vbf_maxpt_j1_isPileUpLoose  = MuonTree->fReader->getBool("JetPFCor_isPileUpJetLoose")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_isPileUpMedium = MuonTree->fReader->getBool("JetPFCor_isPileUpJetMedium")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_isPileUpTight  = MuonTree->fReader->getBool("JetPFCor_isPileUpJetTight")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;

      NewMuonTree->vbf_maxpt_j1_bDiscriminatorSSVHE = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorSSVHE")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_bDiscriminatorTCHE  = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorTCHE")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_bDiscriminatorCSV   = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorCSV")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_bDiscriminatorSSVHP = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorSSVHP")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_bDiscriminatorTCHP  = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorTCHP")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;

      NewMuonTree->vbf_maxpt_j1_ChargedHadronEnergy      = MuonTree->fReader->getFloat("JetPFCor_ChargedHadronEnergy")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_ChargedHadronEnergyFrac  = MuonTree->fReader->getFloat("JetPFCor_ChargedHadronEnergyFrac")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_NeutralHadronEnergy      = MuonTree->fReader->getFloat("JetPFCor_NeutralHadronEnergy")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_NeutralHadronEnergyFrac  = MuonTree->fReader->getFloat("JetPFCor_NeutralHadronEnergyFrac")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_ChargedEmEnergy          = MuonTree->fReader->getFloat("JetPFCor_ChargedEmEnergy")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_ChargedEmEnergyFrac      = MuonTree->fReader->getFloat("JetPFCor_ChargedEmEnergyFrac")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_ChargedMuEnergy          = MuonTree->fReader->getFloat("JetPFCor_ChargedMuEnergy")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_ChargedMuEnergyFrac      = MuonTree->fReader->getFloat("JetPFCor_ChargedMuEnergyFrac")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_NeutralEmEnergy          = MuonTree->fReader->getFloat("JetPFCor_NeutralEmEnergy")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_NeutralEmEnergyFrac      = MuonTree->fReader->getFloat("JetPFCor_NeutralEmEnergyFrac")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_PhotonEnergy             = MuonTree->fReader->getFloat("JetPFCor_PhotonEnergy")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_PhotonEnergyFraction     = MuonTree->fReader->getFloat("JetPFCor_PhotonEnergyFraction")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_ElectronEnergy           = MuonTree->fReader->getFloat("JetPFCor_ElectronEnergy")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_ElectronEnergyFraction   = MuonTree->fReader->getFloat("JetPFCor_ElectronEnergyFraction")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_HFHadronEnergy           = MuonTree->fReader->getFloat("JetPFCor_HFHadronEnergy")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_HFHadronEnergyFraction   = MuonTree->fReader->getFloat("JetPFCor_HFHadronEnergyFraction")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_HFEMEnergy               = MuonTree->fReader->getFloat("JetPFCor_HFEMEnergy")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_HFEMEnergyFraction       = MuonTree->fReader->getFloat("JetPFCor_HFEMEnergyFraction")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;

      NewMuonTree->vbf_maxpt_j1_ChargedMultiplicity        = MuonTree->fReader->getFloat("JetPFCor_ChargedMultiplicity")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_NeutralMultiplicity        = MuonTree->fReader->getFloat("JetPFCor_NeutralMultiplicity")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_MuonMultiplicity           = MuonTree->fReader->getFloat("JetPFCor_MuonMultiplicity")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_ChargedHadronMultiplicity  = MuonTree->fReader->getFloat("JetPFCor_ChargedHadronMultiplicity")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_NeutralHadronMultiplicity  = MuonTree->fReader->getFloat("JetPFCor_NeutralHadronMultiplicity")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_PhotonMultiplicity         = MuonTree->fReader->getFloat("JetPFCor_PhotonMultiplicity")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ; 
      NewMuonTree->vbf_maxpt_j1_ElectronMultiplicity       = MuonTree->fReader->getFloat("JetPFCor_ElectronMultiplicity")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_HFHadronMultiplicity       = MuonTree->fReader->getFloat("JetPFCor_HFHadronMultiplicity")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;


    }
    /// store info only for the VBF case 
    else if (NewMuonTree -> numberJetBin.at(0) >= 2){
    
     // vbf Tag Jet Selection
    
     std::vector<JetAK5> outputAK5_PtSorted;
     std::vector<JetAK5> outputAK5_DEtaSorted;
     std::vector<JetAK5> outputAK5_MjjSorted;

     std::vector<JetAK5> outputGenAK5_PtSorted;
     std::vector<JetAK5> outputGenAK5_DEtaSorted;
     std::vector<JetAK5> outputGenAK5_MjjSorted;

     // Sorting of AK5 Cleaned Collection in Pt

     std::sort(CleanedJetPFCor_AK5_Collection.at(0).begin(),CleanedJetPFCor_AK5_Collection.at(0).end(),TLVP_PtSort());
     outputAK5_PtSorted.push_back(CleanedJetPFCor_AK5_Collection.at(0).at(0));
     outputAK5_PtSorted.push_back(CleanedJetPFCor_AK5_Collection.at(0).at(1));
     if(outputAK5_PtSorted.size() < 2) continue ;
     
     if(!numberJetBinGen.empty() && !GenCleanedJetPFCor_AK5_Collection.at(0).empty()){
      if(numberJetBinGen.at(0)>=2 && GenCleanedJetPFCor_AK5_Collection.at(0).size()>=2){
       std::sort(GenCleanedJetPFCor_AK5_Collection.at(0).begin(),GenCleanedJetPFCor_AK5_Collection.at(0).end(),TLVP_PtSort());
       outputGenAK5_PtSorted.push_back(GenCleanedJetPFCor_AK5_Collection.at(0).at(0));
       outputGenAK5_PtSorted.push_back(GenCleanedJetPFCor_AK5_Collection.at(0).at(1));
      }
     }
         
     // Sorting of AK5 Cleaned Collection in DeltaEta

     std::sort(CleanedJetPFCor_AK5_Collection.at(0).begin(),CleanedJetPFCor_AK5_Collection.at(0).end(),TLVP_EtaSort());
     outputAK5_DEtaSorted.push_back(CleanedJetPFCor_AK5_Collection.at(0).front());
     outputAK5_DEtaSorted.push_back(CleanedJetPFCor_AK5_Collection.at(0).back());
     if(outputAK5_DEtaSorted.size() < 2) continue ;
     
     if(!numberJetBinGen.empty() && !GenCleanedJetPFCor_AK5_Collection.at(0).empty()){
      if(numberJetBinGen.at(0)>=2 && GenCleanedJetPFCor_AK5_Collection.at(0).size()>=2){
       std::sort(GenCleanedJetPFCor_AK5_Collection.at(0).begin(),GenCleanedJetPFCor_AK5_Collection.at(0).end(),TLVP_EtaSort());
       outputGenAK5_DEtaSorted.push_back(GenCleanedJetPFCor_AK5_Collection.at(0).front());
       outputGenAK5_DEtaSorted.push_back(GenCleanedJetPFCor_AK5_Collection.at(0).back());
      }
     }          
     // Sorting of AK5 Cleaned Collection in Mjj
     float maxMjj = 0. ;
     int iJ1 = 0 ;
     int iJ2 = 0 ;

     for (size_t iJet = 0 ; iJet < CleanedJetPFCor_AK5_Collection.at(0).size()-1 ; ++iJet){
       for (size_t jJet = iJet + 1 ; jJet < CleanedJetPFCor_AK5_Collection.at(0).size() ; ++jJet){

        TLorentzVector SumMomentum = CleanedJetPFCor_AK5_Collection.at(0).at(iJet).Momentum_ + CleanedJetPFCor_AK5_Collection.at(0).at(jJet).Momentum_ ;
        float Mjj = SumMomentum.M();
        if(Mjj > maxMjj){
         
          iJ1 = iJet ;
          iJ2 = jJet ;
          maxMjj = Mjj ;
        }
      }
     }

     outputAK5_MjjSorted.push_back (CleanedJetPFCor_AK5_Collection.at(0).at (iJ1)) ;
     outputAK5_MjjSorted.push_back (CleanedJetPFCor_AK5_Collection.at(0).at (iJ2)) ;
     if(outputAK5_MjjSorted.size() < 2) continue ;
    
     maxMjj = 0. ;
     iJ1 = 0 ; iJ2 = 0;

     if(!numberJetBinGen.empty() && !GenCleanedJetPFCor_AK5_Collection.at(0).empty()){
       if(numberJetBinGen.at(0)>=2 && GenCleanedJetPFCor_AK5_Collection.at(0).size()>=2){
       for (size_t iJet = 0 ; iJet < GenCleanedJetPFCor_AK5_Collection.at(0).size()-1 ; ++iJet){
        for (size_t jJet = iJet + 1 ; jJet < GenCleanedJetPFCor_AK5_Collection.at(0).size() ; ++jJet){

         TLorentzVector SumMomentum = GenCleanedJetPFCor_AK5_Collection.at(0).at(iJet).Momentum_ + GenCleanedJetPFCor_AK5_Collection.at(0).at(jJet).Momentum_ ;
         float Mjj = SumMomentum.M();
         if(Mjj > maxMjj){         
          iJ1 = iJet ;
          iJ2 = jJet ;
          maxMjj = Mjj ;
        }
       }
      }
      outputGenAK5_MjjSorted.push_back (GenCleanedJetPFCor_AK5_Collection.at(0).at (iJ1)) ;
      outputGenAK5_MjjSorted.push_back (GenCleanedJetPFCor_AK5_Collection.at(0).at (iJ2)) ;
      }     
     }
      
     //////////////////////////////////////////////////////////////////////////////////////////////////
     // Fill Information for Max Pt Pair of vbf tag jets
     //////////////////////////////////////////////////////////////////////////////////////////////////

     TLorentzVector Total4VMaxPt = outputAK5_PtSorted.at(0).Momentum_ + outputAK5_PtSorted.at(1).Momentum_ ;
    
     NewMuonTree->vbf_maxpt_jj_e   = Total4VMaxPt.E(); 
     NewMuonTree->vbf_maxpt_jj_pt  = Total4VMaxPt.Pt(); 
     NewMuonTree->vbf_maxpt_jj_eta = Total4VMaxPt.Eta(); 
     NewMuonTree->vbf_maxpt_jj_phi = Total4VMaxPt.Phi(); 
     NewMuonTree->vbf_maxpt_jj_m   = Total4VMaxPt.M(); 

     NewMuonTree->vbf_maxpt_j1_e   = outputAK5_PtSorted.at(0).Momentum_.E();
     NewMuonTree->vbf_maxpt_j1_pt  = outputAK5_PtSorted.at(0).Momentum_.Pt();
     NewMuonTree->vbf_maxpt_j1_eta = outputAK5_PtSorted.at(0).Momentum_.Eta();
     NewMuonTree->vbf_maxpt_j1_phi = outputAK5_PtSorted.at(0).Momentum_.Phi();
     NewMuonTree->vbf_maxpt_j1_m   = outputAK5_PtSorted.at(0).Momentum_.M(); 
 
     NewMuonTree->vbf_maxpt_j2_e   = outputAK5_PtSorted.at(1).Momentum_.E();
     NewMuonTree->vbf_maxpt_j2_pt  = outputAK5_PtSorted.at(1).Momentum_.Pt();
     NewMuonTree->vbf_maxpt_j2_eta = outputAK5_PtSorted.at(1).Momentum_.Eta();
     NewMuonTree->vbf_maxpt_j2_phi = outputAK5_PtSorted.at(1).Momentum_.Phi();
     NewMuonTree->vbf_maxpt_j2_m   = outputAK5_PtSorted.at(1).Momentum_.M();

   
     NewMuonTree->vbf_maxpt_jj_deta = fabs(outputAK5_PtSorted.at(0).Momentum_.Eta() - outputAK5_PtSorted.at(1).Momentum_.Eta()) ;
     if (fabs(outputAK5_PtSorted.at(0).Momentum_.Phi() - outputAK5_PtSorted.at(1).Momentum_.Phi()) < TMath::Pi())
       NewMuonTree->vbf_maxpt_jj_dphi = fabs(outputAK5_PtSorted.at(0).Momentum_.Phi() - outputAK5_PtSorted.at(1).Momentum_.Phi()) ;
     else 
       NewMuonTree->vbf_maxpt_jj_dphi = 2*TMath::Pi() - fabs(outputAK5_PtSorted.at(0).Momentum_.Phi() - outputAK5_PtSorted.at(1).Momentum_.Phi()) ;


     if(outputAK5_PtSorted.at(0).NameCollection_ == "JetPFCor" && outputAK5_PtSorted.at(1).NameCollection_ == "JetPFCor"){
       NewMuonTree->vbf_maxpt_jj_type = 1 ; /// both central 

      int nexcj = 0 , nexfj = 0; 
      std::vector<JetAK5>::const_iterator itVec = outputAK5_PtSorted.begin();
      for( ; itVec != outputAK5_PtSorted.end() ; itVec++){

	if(itVec->NameCollection_ == "JetPFCor") nexcj ++ ;
   	if(itVec->NameCollection_ == "JetPFCorVBFTag") nexfj ++ ;
   
      }

      NewMuonTree->vbf_maxpt_n_excj = nexcj ; // number of central jets 
      NewMuonTree->vbf_maxpt_n_exfj = nexfj ; // number of forward jets
     }

     else if(outputAK5_PtSorted.at(0).NameCollection_ == "JetPFCor" && outputAK5_PtSorted.at(1).NameCollection_ == "JetPFCorVBFTag"){
      
      NewMuonTree->vbf_maxpt_jj_type = 2 ;
      int nexcj = 0 , nexfj = 0;
      std::vector<JetAK5>::const_iterator itVec = outputAK5_PtSorted.begin();
      for( ; itVec != outputAK5_PtSorted.end() ; itVec++){

	if(itVec->NameCollection_ == "JetPFCor") nexcj ++ ;
   	if(itVec->NameCollection_ == "JetPFCorVBFTag") nexfj ++ ;
   
      }

      NewMuonTree->vbf_maxpt_n_excj = nexcj ;
      NewMuonTree->vbf_maxpt_n_exfj = nexfj ;
     }

     else if(outputAK5_PtSorted.at(0).NameCollection_ == "JetPFCorVBFTag" && outputAK5_PtSorted.at(1).NameCollection_ == "JetPFCor"){
      
      NewMuonTree->vbf_maxpt_jj_type = 3 ;
      int nexcj = 0 , nexfj = 0;
      std::vector<JetAK5>::const_iterator itVec = outputAK5_PtSorted.begin();
      for( ; itVec != outputAK5_PtSorted.end() ; itVec++){

	if(itVec->NameCollection_ == "JetPFCor") nexcj ++ ;
   	if(itVec->NameCollection_ == "JetPFCorVBFTag") nexfj ++ ;
   
      }

      NewMuonTree->vbf_maxpt_n_excj = nexcj ;
      NewMuonTree->vbf_maxpt_n_exfj = nexfj ;
     }

     else if(outputAK5_PtSorted.at(0).NameCollection_ == "JetPFCorVBFTag" && outputAK5_PtSorted.at(1).NameCollection_ == "JetPFCorVBFTag"){
      
      NewMuonTree->vbf_maxpt_jj_type = 4 ;
      int nexcj = 0 , nexfj = 0;
      std::vector<JetAK5>::const_iterator itVec = outputAK5_PtSorted.begin();
      for( ; itVec != outputAK5_PtSorted.end() ; itVec++){

	if(itVec->NameCollection_ == "JetPFCor") nexcj ++ ;
   	if(itVec->NameCollection_ == "JetPFCorVBFTag") nexfj ++ ;
   
      }

      NewMuonTree->vbf_maxpt_n_excj = nexcj ;
      NewMuonTree->vbf_maxpt_n_exfj = nexfj ;
     }
     else{ std::cerr<<" Something Wrong in MaxPt Jet Categorization "<<std::endl; continue ;}

     if(std::string(SelectionEvents->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEvents->GetXaxis()->SetBinLabel(nStep,"Jet MaxPt Category");
     if(std::string(SelectionEfficiency->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEfficiency->GetXaxis()->SetBinLabel(nStep,"Jet MaxPt Category");
     nstepEvents[nStep-1]++;
     nStep = 7;

     if( NewMuonTree->vbf_maxpt_jj_type < 0 || NewMuonTree->vbf_maxpt_n_excj < 0 || NewMuonTree->vbf_maxpt_n_exfj < 0 ) continue ;

     
     if(!numberJetBinGen.empty() && !outputGenAK5_PtSorted.empty()){
      if(numberJetBinGen.at(0)>=2 && outputGenAK5_PtSorted.size()>=2){

      TLorentzVector Total4VMaxPtGen = outputGenAK5_PtSorted.at(0).Momentum_ + outputGenAK5_PtSorted.at(1).Momentum_ ;
      NewMuonTree->vbf_maxpt_jj_e_gen   = Total4VMaxPtGen.E(); 
      NewMuonTree->vbf_maxpt_jj_pt_gen  = Total4VMaxPtGen.Pt(); 
      NewMuonTree->vbf_maxpt_jj_eta_gen = Total4VMaxPtGen.Eta(); 
      NewMuonTree->vbf_maxpt_jj_phi_gen = Total4VMaxPtGen.Phi(); 
      NewMuonTree->vbf_maxpt_jj_m_gen   = Total4VMaxPtGen.M(); 
      NewMuonTree->vbf_maxpt_j1_e_gen   = outputGenAK5_PtSorted.at(0).Momentum_.E();
      NewMuonTree->vbf_maxpt_j1_pt_gen  = outputGenAK5_PtSorted.at(0).Momentum_.Pt();
      NewMuonTree->vbf_maxpt_j1_eta_gen = outputGenAK5_PtSorted.at(0).Momentum_.Eta();
      NewMuonTree->vbf_maxpt_j1_phi_gen = outputGenAK5_PtSorted.at(0).Momentum_.Phi();
      NewMuonTree->vbf_maxpt_j1_m_gen   = outputGenAK5_PtSorted.at(0).Momentum_.M(); 
      NewMuonTree->vbf_maxpt_j2_e_gen   = outputGenAK5_PtSorted.at(1).Momentum_.E();
      NewMuonTree->vbf_maxpt_j2_pt_gen  = outputGenAK5_PtSorted.at(1).Momentum_.Pt();
      NewMuonTree->vbf_maxpt_j2_eta_gen = outputGenAK5_PtSorted.at(1).Momentum_.Eta();
      NewMuonTree->vbf_maxpt_j2_phi_gen = outputGenAK5_PtSorted.at(1).Momentum_.Phi();
      NewMuonTree->vbf_maxpt_j2_m_gen   = outputGenAK5_PtSorted.at(1).Momentum_.M();

      NewMuonTree->vbf_maxpt_jj_deta_gen = fabs(outputGenAK5_PtSorted.at(0).Momentum_.Eta() - outputGenAK5_PtSorted.at(1).Momentum_.Eta()) ;
      if (fabs(outputGenAK5_PtSorted.at(0).Momentum_.Phi() - outputGenAK5_PtSorted.at(1).Momentum_.Phi()) < TMath::Pi())
       NewMuonTree->vbf_maxpt_jj_dphi_gen = fabs(outputGenAK5_PtSorted.at(0).Momentum_.Phi() - outputGenAK5_PtSorted.at(1).Momentum_.Phi()) ;
      else 
       NewMuonTree->vbf_maxpt_jj_dphi_gen = 2*TMath::Pi() - fabs(outputGenAK5_PtSorted.at(0).Momentum_.Phi() - outputGenAK5_PtSorted.at(1).Momentum_.Phi()) ;

      NewMuonTree->vbf_maxpt_j1_bDiscriminatorSSVHE_gen = MuonTree->fReader->getFloat("JetGen_bDiscriminatorSSVHE")[outputGenAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_bDiscriminatorTCHE_gen  = MuonTree->fReader->getFloat("JetGen_bDiscriminatorTCHE")[outputGenAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_bDiscriminatorCSV_gen   = MuonTree->fReader->getFloat("JetGen_bDiscriminatorCSV")[outputGenAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_bDiscriminatorSSVHP_gen = MuonTree->fReader->getFloat("JetGen_bDiscriminatorSSVHP")[outputGenAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_bDiscriminatorTCHP_gen  = MuonTree->fReader->getFloat("JetGen_bDiscriminatorTCHP")[outputGenAK5_PtSorted.at(0).position_] ;

      NewMuonTree->vbf_maxpt_j2_bDiscriminatorSSVHE_gen = MuonTree->fReader->getFloat("JetGen_bDiscriminatorSSVHE")[outputGenAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_bDiscriminatorTCHE_gen  = MuonTree->fReader->getFloat("JetGen_bDiscriminatorTCHE")[outputGenAK5_PtSorted.at(1).position_] ; 
      NewMuonTree->vbf_maxpt_j2_bDiscriminatorCSV_gen   = MuonTree->fReader->getFloat("JetGen_bDiscriminatorCSV")[outputGenAK5_PtSorted.at(1).position_]  ;
      NewMuonTree->vbf_maxpt_j2_bDiscriminatorSSVHP_gen = MuonTree->fReader->getFloat("JetGen_bDiscriminatorSSVHP")[outputGenAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_bDiscriminatorTCHP_gen  = MuonTree->fReader->getFloat("JetGen_bDiscriminatorTCHP")[outputGenAK5_PtSorted.at(1).position_] ;
      }
      }    
     if(outputAK5_PtSorted.at(0).NameCollection_ == "JetPFCor") {
    
      NewMuonTree->vbf_maxpt_j1_QGLikelihood = MuonTree->fReader->getFloat("JetPFCor_QGLikelihood")[outputAK5_PtSorted.at(0).position_] ;

      NewMuonTree->vbf_maxpt_j1_isPileUpLoose  = MuonTree->fReader->getBool("JetPFCor_isPileUpJetLoose")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_isPileUpMedium = MuonTree->fReader->getBool("JetPFCor_isPileUpJetMedium")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_isPileUpTight  = MuonTree->fReader->getBool("JetPFCor_isPileUpJetTight")[outputAK5_PtSorted.at(0).position_] ;


      NewMuonTree->vbf_maxpt_j1_bDiscriminatorSSVHE = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorSSVHE")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_bDiscriminatorTCHE  = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorTCHE")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_bDiscriminatorCSV   = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorCSV")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_bDiscriminatorSSVHP = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorSSVHP")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_bDiscriminatorTCHP  = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorTCHP")[outputAK5_PtSorted.at(0).position_] ;

      NewMuonTree->vbf_maxpt_j1_ChargedHadronEnergy      = MuonTree->fReader->getFloat("JetPFCor_ChargedHadronEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_ChargedHadronEnergyFrac  = MuonTree->fReader->getFloat("JetPFCor_ChargedHadronEnergyFrac")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_NeutralHadronEnergy      = MuonTree->fReader->getFloat("JetPFCor_NeutralHadronEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_NeutralHadronEnergyFrac  = MuonTree->fReader->getFloat("JetPFCor_NeutralHadronEnergyFrac")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_ChargedEmEnergy          = MuonTree->fReader->getFloat("JetPFCor_ChargedEmEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_ChargedEmEnergyFrac      = MuonTree->fReader->getFloat("JetPFCor_ChargedEmEnergyFrac")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_ChargedMuEnergy          = MuonTree->fReader->getFloat("JetPFCor_ChargedMuEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_ChargedMuEnergyFrac      = MuonTree->fReader->getFloat("JetPFCor_ChargedMuEnergyFrac")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_NeutralEmEnergy          = MuonTree->fReader->getFloat("JetPFCor_NeutralEmEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_NeutralEmEnergyFrac      = MuonTree->fReader->getFloat("JetPFCor_NeutralEmEnergyFrac")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_PhotonEnergy             = MuonTree->fReader->getFloat("JetPFCor_PhotonEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_PhotonEnergyFraction     = MuonTree->fReader->getFloat("JetPFCor_PhotonEnergyFraction")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_ElectronEnergy           = MuonTree->fReader->getFloat("JetPFCor_ElectronEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_ElectronEnergyFraction   = MuonTree->fReader->getFloat("JetPFCor_ElectronEnergyFraction")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_HFHadronEnergy           = MuonTree->fReader->getFloat("JetPFCor_HFHadronEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_HFHadronEnergyFraction   = MuonTree->fReader->getFloat("JetPFCor_HFHadronEnergyFraction")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_HFEMEnergy               = MuonTree->fReader->getFloat("JetPFCor_HFEMEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_HFEMEnergyFraction       = MuonTree->fReader->getFloat("JetPFCor_HFEMEnergyFraction")[outputAK5_PtSorted.at(0).position_] ;

      NewMuonTree->vbf_maxpt_j1_ChargedMultiplicity        = MuonTree->fReader->getFloat("JetPFCor_ChargedMultiplicity")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_NeutralMultiplicity        = MuonTree->fReader->getFloat("JetPFCor_NeutralMultiplicity")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_MuonMultiplicity           = MuonTree->fReader->getFloat("JetPFCor_MuonMultiplicity")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_ChargedHadronMultiplicity  = MuonTree->fReader->getFloat("JetPFCor_ChargedHadronMultiplicity")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_NeutralHadronMultiplicity  = MuonTree->fReader->getFloat("JetPFCor_NeutralHadronMultiplicity")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_PhotonMultiplicity         = MuonTree->fReader->getFloat("JetPFCor_PhotonMultiplicity")[outputAK5_PtSorted.at(0).position_] ; 
      NewMuonTree->vbf_maxpt_j1_ElectronMultiplicity       = MuonTree->fReader->getFloat("JetPFCor_ElectronMultiplicity")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_HFHadronMultiplicity       = MuonTree->fReader->getFloat("JetPFCor_HFHadronMultiplicity")[outputAK5_PtSorted.at(0).position_] ;
     }
     else if(outputAK5_PtSorted.at(0).NameCollection_ == "JetPFCorVBFTag" ){

      NewMuonTree->vbf_maxpt_j1_QGLikelihood = MuonTree->fReader->getFloat("JetPFCorVBFTag_QGLikelihood")[outputAK5_PtSorted.at(0).position_] ;

      NewMuonTree->vbf_maxpt_j1_isPileUpLoose  = MuonTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetLoose")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_isPileUpMedium = MuonTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetMedium")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_isPileUpTight  = MuonTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetTight")[outputAK5_PtSorted.at(0).position_] ;


      NewMuonTree->vbf_maxpt_j1_bDiscriminatorSSVHE = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorSSVHE")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_bDiscriminatorTCHE  = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorTCHE")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_bDiscriminatorCSV   = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorCSV")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_bDiscriminatorSSVHP = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorSSVHP")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_bDiscriminatorTCHP  = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorTCHP")[outputAK5_PtSorted.at(0).position_] ;

      NewMuonTree->vbf_maxpt_j1_ChargedHadronEnergy      = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_ChargedHadronEnergyFrac  = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronEnergyFrac")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_NeutralHadronEnergy      = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_NeutralHadronEnergyFrac  = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronEnergyFrac")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_ChargedEmEnergy          = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedEmEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_ChargedEmEnergyFrac      = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedEmEnergyFrac")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_ChargedMuEnergy          = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedMuEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_ChargedMuEnergyFrac      = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedMuEnergyFrac")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_NeutralEmEnergy          = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralEmEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_NeutralEmEnergyFrac      = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralEmEnergyFrac")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_PhotonEnergy             = MuonTree->fReader->getFloat("JetPFCorVBFTag_PhotonEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_PhotonEnergyFraction     = MuonTree->fReader->getFloat("JetPFCorVBFTag_PhotonEnergyFraction")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_ElectronEnergy           = MuonTree->fReader->getFloat("JetPFCorVBFTag_ElectronEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_ElectronEnergyFraction   = MuonTree->fReader->getFloat("JetPFCorVBFTag_ElectronEnergyFraction")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_HFHadronEnergy           = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFHadronEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_HFHadronEnergyFraction   = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFHadronEnergyFraction")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_HFEMEnergy               = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFEMEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_HFEMEnergyFraction       = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFEMEnergyFraction")[outputAK5_PtSorted.at(0).position_] ;

      NewMuonTree->vbf_maxpt_j1_ChargedMultiplicity        = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedMultiplicity")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_NeutralMultiplicity        = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralMultiplicity")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_MuonMultiplicity           = MuonTree->fReader->getFloat("JetPFCorVBFTag_MuonMultiplicity")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_ChargedHadronMultiplicity  = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronMultiplicity")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_NeutralHadronMultiplicity  = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronMultiplicity")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_PhotonMultiplicity         = MuonTree->fReader->getFloat("JetPFCorVBFTag_PhotonMultiplicity")[outputAK5_PtSorted.at(0).position_] ; 
      NewMuonTree->vbf_maxpt_j1_ElectronMultiplicity       = MuonTree->fReader->getFloat("JetPFCorVBFTag_ElectronMultiplicity")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_HFHadronMultiplicity       = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFHadronMultiplicity")[outputAK5_PtSorted.at(0).position_] ;  

     }
     else { std::cerr<<" problem with High pT Jet Name Collection "<<std::endl; continue ; }

     if(outputAK5_PtSorted.at(1).NameCollection_ == "JetPFCor") {
    
      NewMuonTree->vbf_maxpt_j2_QGLikelihood = MuonTree->fReader->getFloat("JetPFCor_QGLikelihood")[outputAK5_PtSorted.at(1).position_] ;

      NewMuonTree->vbf_maxpt_j2_isPileUpLoose  = MuonTree->fReader->getBool("JetPFCor_isPileUpJetLoose")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_isPileUpMedium = MuonTree->fReader->getBool("JetPFCor_isPileUpJetMedium")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_isPileUpTight  = MuonTree->fReader->getBool("JetPFCor_isPileUpJetTight")[outputAK5_PtSorted.at(1).position_] ;


      NewMuonTree->vbf_maxpt_j2_bDiscriminatorSSVHE = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorSSVHE")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_bDiscriminatorTCHE  = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorTCHE")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_bDiscriminatorCSV   = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorCSV")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_bDiscriminatorSSVHP = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorSSVHP")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_bDiscriminatorTCHP  = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorTCHP")[outputAK5_PtSorted.at(1).position_] ;

      NewMuonTree->vbf_maxpt_j2_ChargedHadronEnergy      = MuonTree->fReader->getFloat("JetPFCor_ChargedHadronEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_ChargedHadronEnergyFrac  = MuonTree->fReader->getFloat("JetPFCor_ChargedHadronEnergyFrac")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_NeutralHadronEnergy      = MuonTree->fReader->getFloat("JetPFCor_NeutralHadronEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_NeutralHadronEnergyFrac  = MuonTree->fReader->getFloat("JetPFCor_NeutralHadronEnergyFrac")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_ChargedEmEnergy          = MuonTree->fReader->getFloat("JetPFCor_ChargedEmEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_ChargedEmEnergyFrac      = MuonTree->fReader->getFloat("JetPFCor_ChargedEmEnergyFrac")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_ChargedMuEnergy          = MuonTree->fReader->getFloat("JetPFCor_ChargedMuEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_ChargedMuEnergyFrac      = MuonTree->fReader->getFloat("JetPFCor_ChargedMuEnergyFrac")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_NeutralEmEnergy          = MuonTree->fReader->getFloat("JetPFCor_NeutralEmEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_NeutralEmEnergyFrac      = MuonTree->fReader->getFloat("JetPFCor_NeutralEmEnergyFrac")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_PhotonEnergy             = MuonTree->fReader->getFloat("JetPFCor_PhotonEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_PhotonEnergyFraction     = MuonTree->fReader->getFloat("JetPFCor_PhotonEnergyFraction")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_ElectronEnergy           = MuonTree->fReader->getFloat("JetPFCor_ElectronEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_ElectronEnergyFraction   = MuonTree->fReader->getFloat("JetPFCor_ElectronEnergyFraction")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_HFHadronEnergy           = MuonTree->fReader->getFloat("JetPFCor_HFHadronEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_HFHadronEnergyFraction   = MuonTree->fReader->getFloat("JetPFCor_HFHadronEnergyFraction")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_HFEMEnergy               = MuonTree->fReader->getFloat("JetPFCor_HFEMEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_HFEMEnergyFraction       = MuonTree->fReader->getFloat("JetPFCor_HFEMEnergyFraction")[outputAK5_PtSorted.at(1).position_] ;

      NewMuonTree->vbf_maxpt_j2_ChargedMultiplicity        = MuonTree->fReader->getFloat("JetPFCor_ChargedMultiplicity")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_NeutralMultiplicity        = MuonTree->fReader->getFloat("JetPFCor_NeutralMultiplicity")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_MuonMultiplicity           = MuonTree->fReader->getFloat("JetPFCor_MuonMultiplicity")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_ChargedHadronMultiplicity  = MuonTree->fReader->getFloat("JetPFCor_ChargedHadronMultiplicity")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_NeutralHadronMultiplicity  = MuonTree->fReader->getFloat("JetPFCor_NeutralHadronMultiplicity")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_PhotonMultiplicity         = MuonTree->fReader->getFloat("JetPFCor_PhotonMultiplicity")[outputAK5_PtSorted.at(1).position_] ; 
      NewMuonTree->vbf_maxpt_j2_ElectronMultiplicity       = MuonTree->fReader->getFloat("JetPFCor_ElectronMultiplicity")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_HFHadronMultiplicity       = MuonTree->fReader->getFloat("JetPFCor_HFHadronMultiplicity")[outputAK5_PtSorted.at(1).position_] ;
     }
     else if(outputAK5_PtSorted.at(1).NameCollection_ == "JetPFCorVBFTag" ){

      NewMuonTree->vbf_maxpt_j2_QGLikelihood = MuonTree->fReader->getFloat("JetPFCorVBFTag_QGLikelihood")[outputAK5_PtSorted.at(1).position_] ;

      NewMuonTree->vbf_maxpt_j2_isPileUpLoose  = MuonTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetLoose")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_isPileUpMedium = MuonTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetMedium")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_isPileUpTight  = MuonTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetTight")[outputAK5_PtSorted.at(1).position_] ;


      NewMuonTree->vbf_maxpt_j2_bDiscriminatorSSVHE = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorSSVHE")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_bDiscriminatorTCHE  = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorTCHE")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_bDiscriminatorCSV   = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorCSV")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_bDiscriminatorSSVHP = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorSSVHP")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_bDiscriminatorTCHP  = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorTCHP")[outputAK5_PtSorted.at(1).position_] ;

      NewMuonTree->vbf_maxpt_j2_ChargedHadronEnergy      = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_ChargedHadronEnergyFrac  = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronEnergyFrac")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_NeutralHadronEnergy      = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_NeutralHadronEnergyFrac  = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronEnergyFrac")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_ChargedEmEnergy          = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedEmEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_ChargedEmEnergyFrac      = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedEmEnergyFrac")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_ChargedMuEnergy          = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedMuEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_ChargedMuEnergyFrac      = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedMuEnergyFrac")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_NeutralEmEnergy          = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralEmEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_NeutralEmEnergyFrac      = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralEmEnergyFrac")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_PhotonEnergy             = MuonTree->fReader->getFloat("JetPFCorVBFTag_PhotonEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_PhotonEnergyFraction     = MuonTree->fReader->getFloat("JetPFCorVBFTag_PhotonEnergyFraction")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_ElectronEnergy           = MuonTree->fReader->getFloat("JetPFCorVBFTag_ElectronEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_ElectronEnergyFraction   = MuonTree->fReader->getFloat("JetPFCorVBFTag_ElectronEnergyFraction")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_HFHadronEnergy           = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFHadronEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_HFHadronEnergyFraction   = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFHadronEnergyFraction")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_HFEMEnergy               = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFEMEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_HFEMEnergyFraction       = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFEMEnergyFraction")[outputAK5_PtSorted.at(1).position_] ;

      NewMuonTree->vbf_maxpt_j2_ChargedMultiplicity        = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedMultiplicity")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_NeutralMultiplicity        = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralMultiplicity")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_MuonMultiplicity           = MuonTree->fReader->getFloat("JetPFCorVBFTag_MuonMultiplicity")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_ChargedHadronMultiplicity  = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronMultiplicity")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_NeutralHadronMultiplicity  = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronMultiplicity")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_PhotonMultiplicity         = MuonTree->fReader->getFloat("JetPFCorVBFTag_PhotonMultiplicity")[outputAK5_PtSorted.at(1).position_] ; 
      NewMuonTree->vbf_maxpt_j2_ElectronMultiplicity       = MuonTree->fReader->getFloat("JetPFCorVBFTag_ElectronMultiplicity")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_HFHadronMultiplicity       = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFHadronMultiplicity")[outputAK5_PtSorted.at(1).position_] ;  

     }
     else { std::cerr<<" problem with High pT Jet Name Collection "<<std::endl; continue ; }

     /////////////////////////////////////////////////////////////////////////////////////////////////
     // Fill Information for Max Deta Pair of vbf tag jets
     //////////////////////////////////////////////////////////////////////////////////////////////////
   
     TLorentzVector Total4VMaxDeta = outputAK5_DEtaSorted.at(0).Momentum_ + outputAK5_DEtaSorted.at(1).Momentum_ ;
    
     NewMuonTree->vbf_maxDeta_jj_e   = Total4VMaxDeta.E(); 
     NewMuonTree->vbf_maxDeta_jj_pt  = Total4VMaxDeta.Pt(); 
     NewMuonTree->vbf_maxDeta_jj_eta = Total4VMaxDeta.Eta(); 
     NewMuonTree->vbf_maxDeta_jj_phi = Total4VMaxDeta.Phi(); 
     NewMuonTree->vbf_maxDeta_jj_m   = Total4VMaxDeta.M(); 
 
     NewMuonTree->vbf_maxDeta_j1_e   = outputAK5_DEtaSorted.at(0).Momentum_.E();
     NewMuonTree->vbf_maxDeta_j1_pt  = outputAK5_DEtaSorted.at(0).Momentum_.Pt();
     NewMuonTree->vbf_maxDeta_j1_eta = outputAK5_DEtaSorted.at(0).Momentum_.Eta();
     NewMuonTree->vbf_maxDeta_j1_phi = outputAK5_DEtaSorted.at(0).Momentum_.Phi();
     NewMuonTree->vbf_maxDeta_j1_m   = outputAK5_DEtaSorted.at(0).Momentum_.M();

     NewMuonTree->vbf_maxDeta_j2_e   = outputAK5_DEtaSorted.at(1).Momentum_.E();
     NewMuonTree->vbf_maxDeta_j2_pt  = outputAK5_DEtaSorted.at(1).Momentum_.Pt();
     NewMuonTree->vbf_maxDeta_j2_eta = outputAK5_DEtaSorted.at(1).Momentum_.Eta();
     NewMuonTree->vbf_maxDeta_j2_phi = outputAK5_DEtaSorted.at(1).Momentum_.Phi();
     NewMuonTree->vbf_maxDeta_j2_m   = outputAK5_DEtaSorted.at(1).Momentum_.M();
   
     NewMuonTree->vbf_maxDeta_jj_deta = fabs(outputAK5_DEtaSorted.at(0).Momentum_.Eta() - outputAK5_DEtaSorted.at(1).Momentum_.Eta()) ;
     if (fabs(outputAK5_DEtaSorted.at(0).Momentum_.Phi() - outputAK5_DEtaSorted.at(1).Momentum_.Phi()) < 3.14)
      NewMuonTree->vbf_maxDeta_jj_dphi = fabs(outputAK5_DEtaSorted.at(0).Momentum_.Phi() - outputAK5_DEtaSorted.at(1).Momentum_.Phi()) ;
     else  
      NewMuonTree->vbf_maxDeta_jj_dphi = 6.28 - fabs(outputAK5_DEtaSorted.at(0).Momentum_.Phi() - outputAK5_DEtaSorted.at(1).Momentum_.Phi()) ;

     if(outputAK5_DEtaSorted.at(0).NameCollection_ == "JetPFCor" && outputAK5_DEtaSorted.at(1).NameCollection_ == "JetPFCor"){
      
      NewMuonTree->vbf_maxDeta_jj_type = 1 ;

      int nexcj = 0 , nexfj = 0; 
      std::vector<JetAK5>::const_iterator itVec = outputAK5_DEtaSorted.begin();
      for( ; itVec != outputAK5_DEtaSorted.end() ; itVec++){

	if(itVec->NameCollection_ == "JetPFCor") nexcj ++ ;
   	if(itVec->NameCollection_ == "JetPFCorVBFTag") nexfj ++ ;
   
      }

      NewMuonTree->vbf_maxDeta_n_excj = nexcj ;
      NewMuonTree->vbf_maxDeta_n_exfj = nexfj ;
     }

     else if(outputAK5_DEtaSorted.at(0).NameCollection_ == "JetPFCor" && outputAK5_DEtaSorted.at(1).NameCollection_ == "JetPFCorVBFTag"){
      
      NewMuonTree->vbf_maxDeta_jj_type = 2 ;

      int nexcj = 0 , nexfj = 0;
      std::vector<JetAK5>::const_iterator itVec = outputAK5_DEtaSorted.begin();
      for( ; itVec != outputAK5_DEtaSorted.end() ; itVec++){

	if(itVec->NameCollection_ == "JetPFCor") nexcj ++ ;
   	if(itVec->NameCollection_ == "JetPFCorVBFTag") nexfj ++ ;
   
      }

      NewMuonTree->vbf_maxDeta_n_excj = nexcj ;
      NewMuonTree->vbf_maxDeta_n_exfj = nexfj ;
     }

     else if(outputAK5_DEtaSorted.at(0).NameCollection_ == "JetPFCorVBFTag" && outputAK5_DEtaSorted.at(1).NameCollection_ == "JetPFCor"){
      
      NewMuonTree->vbf_maxDeta_jj_type = 3 ;

      int nexcj = 0 , nexfj = 0;
      std::vector<JetAK5>::const_iterator itVec = outputAK5_DEtaSorted.begin();
      for( ; itVec != outputAK5_DEtaSorted.end() ; itVec++){

	if(itVec->NameCollection_ == "JetPFCor") nexcj ++ ;
   	if(itVec->NameCollection_ == "JetPFCorVBFTag") nexfj ++ ;
   
      }

      NewMuonTree->vbf_maxDeta_n_excj = nexcj ;
      NewMuonTree->vbf_maxDeta_n_exfj = nexfj ;
     }

     else if(outputAK5_DEtaSorted.at(0).NameCollection_ == "JetPFCorVBFTag" && outputAK5_DEtaSorted.at(1).NameCollection_ == "JetPFCorVBFTag"){
      
      NewMuonTree->vbf_maxDeta_jj_type = 4 ;

      int nexcj = 0 , nexfj = 0;
      std::vector<JetAK5>::const_iterator itVec = outputAK5_DEtaSorted.begin();
      for( ; itVec != outputAK5_DEtaSorted.end() ; itVec++){

	if(itVec->NameCollection_ == "JetPFCor") nexcj ++ ;
   	if(itVec->NameCollection_ == "JetPFCorVBFTag") nexfj ++ ;
   
      }

      NewMuonTree->vbf_maxDeta_n_excj = nexcj ;
      NewMuonTree->vbf_maxDeta_n_exfj = nexfj ;
     }
     else{ std::cerr<<" Something Wrong in MaxDeta Jet Categorization "<<std::endl; continue ;}

     if(std::string(SelectionEvents->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEvents->GetXaxis()->SetBinLabel(nStep,"Jet MaxDeta Category");
     if(std::string(SelectionEfficiency->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEfficiency->GetXaxis()->SetBinLabel(nStep,"Jet MaxDeta Category");
     nstepEvents[nStep-1]++;
     nStep = 8;

     
     if(!numberJetBinGen.empty() && !outputGenAK5_DEtaSorted.empty()){
      if(numberJetBinGen.at(0)>=2 && outputGenAK5_DEtaSorted.size()>=2){

      TLorentzVector Total4VMaxDetaGen = outputGenAK5_DEtaSorted.at(0).Momentum_ + outputGenAK5_DEtaSorted.at(1).Momentum_ ;
      NewMuonTree->vbf_maxDeta_jj_e_gen   = Total4VMaxDetaGen.E(); 
      NewMuonTree->vbf_maxDeta_jj_pt_gen  = Total4VMaxDetaGen.Pt(); 
      NewMuonTree->vbf_maxDeta_jj_eta_gen = Total4VMaxDetaGen.Eta(); 
      NewMuonTree->vbf_maxDeta_jj_phi_gen = Total4VMaxDetaGen.Phi(); 
      NewMuonTree->vbf_maxDeta_jj_m_gen   = Total4VMaxDetaGen.M(); 
      NewMuonTree->vbf_maxDeta_j1_e_gen   = outputGenAK5_DEtaSorted.at(0).Momentum_.E();
      NewMuonTree->vbf_maxDeta_j1_pt_gen  = outputGenAK5_DEtaSorted.at(0).Momentum_.Pt();
      NewMuonTree->vbf_maxDeta_j1_eta_gen = outputGenAK5_DEtaSorted.at(0).Momentum_.Eta();
      NewMuonTree->vbf_maxDeta_j1_phi_gen = outputGenAK5_DEtaSorted.at(0).Momentum_.Phi();
      NewMuonTree->vbf_maxDeta_j1_m_gen   = outputGenAK5_DEtaSorted.at(0).Momentum_.M(); 
      NewMuonTree->vbf_maxDeta_j2_e_gen   = outputGenAK5_DEtaSorted.at(1).Momentum_.E();
      NewMuonTree->vbf_maxDeta_j2_pt_gen  = outputGenAK5_DEtaSorted.at(1).Momentum_.Pt();
      NewMuonTree->vbf_maxDeta_j2_eta_gen = outputGenAK5_DEtaSorted.at(1).Momentum_.Eta();
      NewMuonTree->vbf_maxDeta_j2_phi_gen = outputGenAK5_DEtaSorted.at(1).Momentum_.Phi();
      NewMuonTree->vbf_maxDeta_j2_m_gen   = outputGenAK5_DEtaSorted.at(1).Momentum_.M();

      NewMuonTree->vbf_maxDeta_jj_deta_gen = fabs(outputGenAK5_DEtaSorted.at(0).Momentum_.Eta() - outputGenAK5_DEtaSorted.at(1).Momentum_.Eta()) ;
      if (fabs(outputGenAK5_DEtaSorted.at(0).Momentum_.Phi() - outputGenAK5_DEtaSorted.at(1).Momentum_.Phi()) < TMath::Pi())
       NewMuonTree->vbf_maxDeta_jj_dphi_gen = fabs(outputGenAK5_DEtaSorted.at(0).Momentum_.Phi() - outputGenAK5_DEtaSorted.at(1).Momentum_.Phi()) ;
      else 
       NewMuonTree->vbf_maxDeta_jj_dphi_gen = 2*TMath::Pi() - fabs(outputGenAK5_DEtaSorted.at(0).Momentum_.Phi() - outputGenAK5_DEtaSorted.at(1).Momentum_.Phi()) ;

      NewMuonTree->vbf_maxDeta_j1_bDiscriminatorSSVHE_gen = MuonTree->fReader->getFloat("JetGen_bDiscriminatorSSVHE")[outputGenAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_bDiscriminatorTCHE_gen  = MuonTree->fReader->getFloat("JetGen_bDiscriminatorTCHE")[outputGenAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_bDiscriminatorCSV_gen   = MuonTree->fReader->getFloat("JetGen_bDiscriminatorCSV")[outputGenAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_bDiscriminatorSSVHP_gen = MuonTree->fReader->getFloat("JetGen_bDiscriminatorSSVHP")[outputGenAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_bDiscriminatorTCHP_gen  = MuonTree->fReader->getFloat("JetGen_bDiscriminatorTCHP")[outputGenAK5_DEtaSorted.at(0).position_] ;

      NewMuonTree->vbf_maxDeta_j2_bDiscriminatorSSVHE_gen = MuonTree->fReader->getFloat("JetGen_bDiscriminatorSSVHE")[outputGenAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_bDiscriminatorTCHE_gen  = MuonTree->fReader->getFloat("JetGen_bDiscriminatorTCHE")[outputGenAK5_DEtaSorted.at(1).position_] ; 
      NewMuonTree->vbf_maxDeta_j2_bDiscriminatorCSV_gen   = MuonTree->fReader->getFloat("JetGen_bDiscriminatorCSV")[outputGenAK5_DEtaSorted.at(1).position_]  ;
      NewMuonTree->vbf_maxDeta_j2_bDiscriminatorSSVHP_gen = MuonTree->fReader->getFloat("JetGen_bDiscriminatorSSVHP")[outputGenAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_bDiscriminatorTCHP_gen  = MuonTree->fReader->getFloat("JetGen_bDiscriminatorTCHP")[outputGenAK5_DEtaSorted.at(1).position_] ;
     }
    }

    if(outputAK5_DEtaSorted.at(0).NameCollection_ == "JetPFCor") {
    
      NewMuonTree->vbf_maxDeta_j1_QGLikelihood = MuonTree->fReader->getFloat("JetPFCor_QGLikelihood")[outputAK5_DEtaSorted.at(0).position_] ;

      NewMuonTree->vbf_maxDeta_j1_isPileUpLoose  = MuonTree->fReader->getBool("JetPFCor_isPileUpJetLoose")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_isPileUpMedium = MuonTree->fReader->getBool("JetPFCor_isPileUpJetMedium")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_isPileUpTight  = MuonTree->fReader->getBool("JetPFCor_isPileUpJetTight")[outputAK5_DEtaSorted.at(0).position_] ;


      NewMuonTree->vbf_maxDeta_j1_bDiscriminatorSSVHE = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorSSVHE")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_bDiscriminatorTCHE  = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorTCHE")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_bDiscriminatorCSV   = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorCSV")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_bDiscriminatorSSVHP = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorSSVHP")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_bDiscriminatorTCHP  = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorTCHP")[outputAK5_DEtaSorted.at(0).position_] ;

      NewMuonTree->vbf_maxDeta_j1_ChargedHadronEnergy      = MuonTree->fReader->getFloat("JetPFCor_ChargedHadronEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_ChargedHadronEnergyFrac  = MuonTree->fReader->getFloat("JetPFCor_ChargedHadronEnergyFrac")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_NeutralHadronEnergy      = MuonTree->fReader->getFloat("JetPFCor_NeutralHadronEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_NeutralHadronEnergyFrac  = MuonTree->fReader->getFloat("JetPFCor_NeutralHadronEnergyFrac")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_ChargedEmEnergy          = MuonTree->fReader->getFloat("JetPFCor_ChargedEmEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_ChargedEmEnergyFrac      = MuonTree->fReader->getFloat("JetPFCor_ChargedEmEnergyFrac")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_ChargedMuEnergy          = MuonTree->fReader->getFloat("JetPFCor_ChargedMuEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_ChargedMuEnergyFrac      = MuonTree->fReader->getFloat("JetPFCor_ChargedMuEnergyFrac")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_NeutralEmEnergy          = MuonTree->fReader->getFloat("JetPFCor_NeutralEmEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_NeutralEmEnergyFrac      = MuonTree->fReader->getFloat("JetPFCor_NeutralEmEnergyFrac")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_PhotonEnergy             = MuonTree->fReader->getFloat("JetPFCor_PhotonEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_PhotonEnergyFraction     = MuonTree->fReader->getFloat("JetPFCor_PhotonEnergyFraction")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_ElectronEnergy           = MuonTree->fReader->getFloat("JetPFCor_ElectronEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_ElectronEnergyFraction   = MuonTree->fReader->getFloat("JetPFCor_ElectronEnergyFraction")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_HFHadronEnergy           = MuonTree->fReader->getFloat("JetPFCor_HFHadronEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_HFHadronEnergyFraction   = MuonTree->fReader->getFloat("JetPFCor_HFHadronEnergyFraction")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_HFEMEnergy               = MuonTree->fReader->getFloat("JetPFCor_HFEMEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_HFEMEnergyFraction       = MuonTree->fReader->getFloat("JetPFCor_HFEMEnergyFraction")[outputAK5_DEtaSorted.at(0).position_] ;

      NewMuonTree->vbf_maxDeta_j1_ChargedMultiplicity        = MuonTree->fReader->getFloat("JetPFCor_ChargedMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_NeutralMultiplicity        = MuonTree->fReader->getFloat("JetPFCor_NeutralMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_MuonMultiplicity           = MuonTree->fReader->getFloat("JetPFCor_MuonMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_ChargedHadronMultiplicity  = MuonTree->fReader->getFloat("JetPFCor_ChargedHadronMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_NeutralHadronMultiplicity  = MuonTree->fReader->getFloat("JetPFCor_NeutralHadronMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_PhotonMultiplicity         = MuonTree->fReader->getFloat("JetPFCor_PhotonMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ; 
      NewMuonTree->vbf_maxDeta_j1_ElectronMultiplicity       = MuonTree->fReader->getFloat("JetPFCor_ElectronMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_HFHadronMultiplicity       = MuonTree->fReader->getFloat("JetPFCor_HFHadronMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ;
     }
     else if(outputAK5_DEtaSorted.at(0).NameCollection_ == "JetPFCorVBFTag" ){

      NewMuonTree->vbf_maxDeta_j1_QGLikelihood = MuonTree->fReader->getFloat("JetPFCorVBFTag_QGLikelihood")[outputAK5_DEtaSorted.at(0).position_] ;

      NewMuonTree->vbf_maxDeta_j1_isPileUpLoose  = MuonTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetLoose")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_isPileUpMedium = MuonTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetMedium")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_isPileUpTight  = MuonTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetTight")[outputAK5_DEtaSorted.at(0).position_] ;


      NewMuonTree->vbf_maxDeta_j1_bDiscriminatorSSVHE = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorSSVHE")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_bDiscriminatorTCHE  = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorTCHE")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_bDiscriminatorCSV   = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorCSV")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_bDiscriminatorSSVHP = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorSSVHP")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_bDiscriminatorTCHP  = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorTCHP")[outputAK5_DEtaSorted.at(0).position_] ;

      NewMuonTree->vbf_maxDeta_j1_ChargedHadronEnergy      = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_ChargedHadronEnergyFrac  = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronEnergyFrac")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_NeutralHadronEnergy      = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_NeutralHadronEnergyFrac  = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronEnergyFrac")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_ChargedEmEnergy          = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedEmEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_ChargedEmEnergyFrac      = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedEmEnergyFrac")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_ChargedMuEnergy          = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedMuEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_ChargedMuEnergyFrac      = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedMuEnergyFrac")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_NeutralEmEnergy          = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralEmEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_NeutralEmEnergyFrac      = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralEmEnergyFrac")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_PhotonEnergy             = MuonTree->fReader->getFloat("JetPFCorVBFTag_PhotonEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_PhotonEnergyFraction     = MuonTree->fReader->getFloat("JetPFCorVBFTag_PhotonEnergyFraction")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_ElectronEnergy           = MuonTree->fReader->getFloat("JetPFCorVBFTag_ElectronEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_ElectronEnergyFraction   = MuonTree->fReader->getFloat("JetPFCorVBFTag_ElectronEnergyFraction")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_HFHadronEnergy           = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFHadronEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_HFHadronEnergyFraction   = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFHadronEnergyFraction")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_HFEMEnergy               = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFEMEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_HFEMEnergyFraction       = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFEMEnergyFraction")[outputAK5_DEtaSorted.at(0).position_] ;

      NewMuonTree->vbf_maxDeta_j1_ChargedMultiplicity        = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_NeutralMultiplicity        = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_MuonMultiplicity           = MuonTree->fReader->getFloat("JetPFCorVBFTag_MuonMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_ChargedHadronMultiplicity  = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_NeutralHadronMultiplicity  = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_PhotonMultiplicity         = MuonTree->fReader->getFloat("JetPFCorVBFTag_PhotonMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ; 
      NewMuonTree->vbf_maxDeta_j1_ElectronMultiplicity       = MuonTree->fReader->getFloat("JetPFCorVBFTag_ElectronMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_HFHadronMultiplicity       = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFHadronMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ;  

     }
     else { std::cerr<<" problem with High Deta Jet Name Collection "<<std::endl; continue ; }

     if(outputAK5_DEtaSorted.at(1).NameCollection_ == "JetPFCor") {
    
      NewMuonTree->vbf_maxDeta_j2_QGLikelihood = MuonTree->fReader->getFloat("JetPFCor_QGLikelihood")[outputAK5_DEtaSorted.at(1).position_] ;

      NewMuonTree->vbf_maxDeta_j2_isPileUpLoose  = MuonTree->fReader->getBool("JetPFCor_isPileUpJetLoose")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_isPileUpMedium = MuonTree->fReader->getBool("JetPFCor_isPileUpJetMedium")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_isPileUpTight  = MuonTree->fReader->getBool("JetPFCor_isPileUpJetTight")[outputAK5_DEtaSorted.at(1).position_] ;


      NewMuonTree->vbf_maxDeta_j2_bDiscriminatorSSVHE = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorSSVHE")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_bDiscriminatorTCHE  = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorTCHE")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_bDiscriminatorCSV   = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorCSV")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_bDiscriminatorSSVHP = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorSSVHP")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_bDiscriminatorTCHP  = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorTCHP")[outputAK5_DEtaSorted.at(1).position_] ;

      NewMuonTree->vbf_maxDeta_j2_ChargedHadronEnergy      = MuonTree->fReader->getFloat("JetPFCor_ChargedHadronEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_ChargedHadronEnergyFrac  = MuonTree->fReader->getFloat("JetPFCor_ChargedHadronEnergyFrac")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_NeutralHadronEnergy      = MuonTree->fReader->getFloat("JetPFCor_NeutralHadronEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_NeutralHadronEnergyFrac  = MuonTree->fReader->getFloat("JetPFCor_NeutralHadronEnergyFrac")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_ChargedEmEnergy          = MuonTree->fReader->getFloat("JetPFCor_ChargedEmEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_ChargedEmEnergyFrac      = MuonTree->fReader->getFloat("JetPFCor_ChargedEmEnergyFrac")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_ChargedMuEnergy          = MuonTree->fReader->getFloat("JetPFCor_ChargedMuEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_ChargedMuEnergyFrac      = MuonTree->fReader->getFloat("JetPFCor_ChargedMuEnergyFrac")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_NeutralEmEnergy          = MuonTree->fReader->getFloat("JetPFCor_NeutralEmEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_NeutralEmEnergyFrac      = MuonTree->fReader->getFloat("JetPFCor_NeutralEmEnergyFrac")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_PhotonEnergy             = MuonTree->fReader->getFloat("JetPFCor_PhotonEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_PhotonEnergyFraction     = MuonTree->fReader->getFloat("JetPFCor_PhotonEnergyFraction")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_ElectronEnergy           = MuonTree->fReader->getFloat("JetPFCor_ElectronEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_ElectronEnergyFraction   = MuonTree->fReader->getFloat("JetPFCor_ElectronEnergyFraction")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_HFHadronEnergy           = MuonTree->fReader->getFloat("JetPFCor_HFHadronEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_HFHadronEnergyFraction   = MuonTree->fReader->getFloat("JetPFCor_HFHadronEnergyFraction")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_HFEMEnergy               = MuonTree->fReader->getFloat("JetPFCor_HFEMEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_HFEMEnergyFraction       = MuonTree->fReader->getFloat("JetPFCor_HFEMEnergyFraction")[outputAK5_DEtaSorted.at(1).position_] ;

      NewMuonTree->vbf_maxDeta_j2_ChargedMultiplicity        = MuonTree->fReader->getFloat("JetPFCor_ChargedMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_NeutralMultiplicity        = MuonTree->fReader->getFloat("JetPFCor_NeutralMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_MuonMultiplicity           = MuonTree->fReader->getFloat("JetPFCor_MuonMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_ChargedHadronMultiplicity  = MuonTree->fReader->getFloat("JetPFCor_ChargedHadronMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_NeutralHadronMultiplicity  = MuonTree->fReader->getFloat("JetPFCor_NeutralHadronMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_PhotonMultiplicity         = MuonTree->fReader->getFloat("JetPFCor_PhotonMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ; 
      NewMuonTree->vbf_maxDeta_j2_ElectronMultiplicity       = MuonTree->fReader->getFloat("JetPFCor_ElectronMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_HFHadronMultiplicity       = MuonTree->fReader->getFloat("JetPFCor_HFHadronMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ;
     }
     else if(outputAK5_DEtaSorted.at(1).NameCollection_ == "JetPFCorVBFTag" ){

      NewMuonTree->vbf_maxDeta_j2_QGLikelihood = MuonTree->fReader->getFloat("JetPFCorVBFTag_QGLikelihood")[outputAK5_DEtaSorted.at(1).position_] ;

      NewMuonTree->vbf_maxDeta_j2_isPileUpLoose = MuonTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetLoose")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_isPileUpMedium = MuonTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetMedium")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_isPileUpTight = MuonTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetTight")[outputAK5_DEtaSorted.at(1).position_] ;


      NewMuonTree->vbf_maxDeta_j2_bDiscriminatorSSVHE = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorSSVHE")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_bDiscriminatorTCHE  = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorTCHE")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_bDiscriminatorCSV   = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorCSV")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_bDiscriminatorSSVHP = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorSSVHP")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_bDiscriminatorTCHP  = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorTCHP")[outputAK5_DEtaSorted.at(1).position_] ;

      NewMuonTree->vbf_maxDeta_j2_ChargedHadronEnergy      = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_ChargedHadronEnergyFrac  = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronEnergyFrac")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_NeutralHadronEnergy      = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_NeutralHadronEnergyFrac  = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronEnergyFrac")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_ChargedEmEnergy          = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedEmEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_ChargedEmEnergyFrac      = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedEmEnergyFrac")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_ChargedMuEnergy          = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedMuEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_ChargedMuEnergyFrac      = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedMuEnergyFrac")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_NeutralEmEnergy          = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralEmEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_NeutralEmEnergyFrac      = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralEmEnergyFrac")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_PhotonEnergy             = MuonTree->fReader->getFloat("JetPFCorVBFTag_PhotonEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_PhotonEnergyFraction     = MuonTree->fReader->getFloat("JetPFCorVBFTag_PhotonEnergyFraction")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_ElectronEnergy           = MuonTree->fReader->getFloat("JetPFCorVBFTag_ElectronEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_ElectronEnergyFraction   = MuonTree->fReader->getFloat("JetPFCorVBFTag_ElectronEnergyFraction")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_HFHadronEnergy           = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFHadronEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_HFHadronEnergyFraction   = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFHadronEnergyFraction")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_HFEMEnergy               = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFEMEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_HFEMEnergyFraction       = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFEMEnergyFraction")[outputAK5_DEtaSorted.at(1).position_] ;

      NewMuonTree->vbf_maxDeta_j2_ChargedMultiplicity  = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_NeutralMultiplicity  = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_MuonMultiplicity     = MuonTree->fReader->getFloat("JetPFCorVBFTag_MuonMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_ChargedHadronMultiplicity  = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_NeutralHadronMultiplicity  = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_PhotonMultiplicity         = MuonTree->fReader->getFloat("JetPFCorVBFTag_PhotonMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ; 
      NewMuonTree->vbf_maxDeta_j2_ElectronMultiplicity       = MuonTree->fReader->getFloat("JetPFCorVBFTag_ElectronMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_HFHadronMultiplicity       = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFHadronMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ;  

     }
     else { std::cerr<<" problem with High Deta Jet Name Collection "<<std::endl; continue ; }
    
     /////////////////////////////////////////////////////////////////////////////////////////////////
     // Fill Information for Max Mjj Pair of vbf tag jets
     //////////////////////////////////////////////////////////////////////////////////////////////////

     TLorentzVector Total4VMaxMjj = outputAK5_MjjSorted.at(0).Momentum_ + outputAK5_MjjSorted.at(1).Momentum_ ;
    
     NewMuonTree->vbf_maxMjj_jj_e   = Total4VMaxMjj.E(); 
     NewMuonTree->vbf_maxMjj_jj_pt  = Total4VMaxMjj.Pt(); 
     NewMuonTree->vbf_maxMjj_jj_eta = Total4VMaxMjj.Eta(); 
     NewMuonTree->vbf_maxMjj_jj_phi = Total4VMaxMjj.Phi(); 
     NewMuonTree->vbf_maxMjj_jj_m   = Total4VMaxMjj.M(); 
 
     NewMuonTree->vbf_maxMjj_j1_e   = outputAK5_MjjSorted.at(0).Momentum_.E();
     NewMuonTree->vbf_maxMjj_j1_pt  = outputAK5_MjjSorted.at(0).Momentum_.Pt();
     NewMuonTree->vbf_maxMjj_j1_eta = outputAK5_MjjSorted.at(0).Momentum_.Eta();
     NewMuonTree->vbf_maxMjj_j1_phi = outputAK5_MjjSorted.at(0).Momentum_.Phi();
     NewMuonTree->vbf_maxMjj_j1_m   = outputAK5_MjjSorted.at(0).Momentum_.M();

     NewMuonTree->vbf_maxMjj_j2_e   = outputAK5_MjjSorted.at(1).Momentum_.E();
     NewMuonTree->vbf_maxMjj_j2_pt  = outputAK5_MjjSorted.at(1).Momentum_.Pt();
     NewMuonTree->vbf_maxMjj_j2_eta = outputAK5_MjjSorted.at(1).Momentum_.Eta();
     NewMuonTree->vbf_maxMjj_j2_phi = outputAK5_MjjSorted.at(1).Momentum_.Phi();
     NewMuonTree->vbf_maxMjj_j2_m   = outputAK5_MjjSorted.at(1).Momentum_.M();
   
     NewMuonTree->vbf_maxMjj_jj_deta = fabs(outputAK5_MjjSorted.at(0).Momentum_.Eta() - outputAK5_MjjSorted.at(1).Momentum_.Eta()) ;
     if (fabs(outputAK5_MjjSorted.at(0).Momentum_.Phi() - outputAK5_MjjSorted.at(1).Momentum_.Phi())<3.14)
      NewMuonTree->vbf_maxMjj_jj_dphi = fabs(outputAK5_MjjSorted.at(0).Momentum_.Phi() - outputAK5_MjjSorted.at(1).Momentum_.Phi()) ;
     else
      NewMuonTree->vbf_maxMjj_jj_dphi = 6.28 - fabs(outputAK5_MjjSorted.at(0).Momentum_.Phi() - outputAK5_MjjSorted.at(1).Momentum_.Phi()) ;

     if(outputAK5_MjjSorted.at(0).NameCollection_ == "JetPFCor" && outputAK5_MjjSorted.at(1).NameCollection_ == "JetPFCor"){
      
      NewMuonTree->vbf_maxMjj_jj_type = 1 ;

      int nexcj = 0 , nexfj = 0; 
      std::vector<JetAK5>::const_iterator itVec = outputAK5_MjjSorted.begin();
      for( ; itVec != outputAK5_MjjSorted.end() ; itVec++){

	if(itVec->NameCollection_ == "JetPFCor") nexcj ++ ;
   	if(itVec->NameCollection_ == "JetPFCorVBFTag") nexfj ++ ;
   
      }

      NewMuonTree->vbf_maxMjj_n_excj = nexcj ;
      NewMuonTree->vbf_maxMjj_n_exfj = nexfj ;
     }

     else if(outputAK5_MjjSorted.at(0).NameCollection_ == "JetPFCor" && outputAK5_MjjSorted.at(1).NameCollection_ == "JetPFCorVBFTag"){
      
      NewMuonTree->vbf_maxMjj_jj_type = 2 ;

      int nexcj = 0 , nexfj = 0;
      std::vector<JetAK5>::const_iterator itVec = outputAK5_MjjSorted.begin();
      for( ; itVec != outputAK5_MjjSorted.end() ; itVec++){

	if(itVec->NameCollection_ == "JetPFCor") nexcj ++ ;
   	if(itVec->NameCollection_ == "JetPFCorVBFTag") nexfj ++ ;
   
      }

      NewMuonTree->vbf_maxMjj_n_excj = nexcj ;
      NewMuonTree->vbf_maxMjj_n_exfj = nexfj ;
     }

     else if(outputAK5_MjjSorted.at(0).NameCollection_ == "JetPFCorVBFTag" && outputAK5_MjjSorted.at(1).NameCollection_ == "JetPFCor"){
      
      NewMuonTree->vbf_maxMjj_jj_type = 3 ;

      int nexcj = 0 , nexfj = 0;
      std::vector<JetAK5>::const_iterator itVec = outputAK5_MjjSorted.begin();
      for( ; itVec != outputAK5_MjjSorted.end() ; itVec++){

	if(itVec->NameCollection_ == "JetPFCor") nexcj ++ ;
   	if(itVec->NameCollection_ == "JetPFCorVBFTag") nexfj ++ ;
   
      }

      NewMuonTree->vbf_maxMjj_n_excj = nexcj ;
      NewMuonTree->vbf_maxMjj_n_exfj = nexfj ;
     }

     else if(outputAK5_MjjSorted.at(0).NameCollection_ == "JetPFCorVBFTag" && outputAK5_MjjSorted.at(1).NameCollection_ == "JetPFCorVBFTag"){
      
      NewMuonTree->vbf_maxMjj_jj_type = 4 ;

      int nexcj = 0 , nexfj = 0;
      std::vector<JetAK5>::const_iterator itVec = outputAK5_MjjSorted.begin();
      for( ; itVec != outputAK5_MjjSorted.end() ; itVec++){

	if(itVec->NameCollection_ == "JetPFCor") nexcj ++ ;
   	if(itVec->NameCollection_ == "JetPFCorVBFTag") nexfj ++ ;
   
      }

      NewMuonTree->vbf_maxMjj_n_excj = nexcj ;
      NewMuonTree->vbf_maxMjj_n_exfj = nexfj ;
     }
     else{ std::cerr<<" Something Wrong in MaxMjj Jet Categorization "<<std::endl; continue ;}

     if(std::string(SelectionEvents->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEvents->GetXaxis()->SetBinLabel(nStep,"Jet MaxMjj Category");
     if(std::string(SelectionEfficiency->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEfficiency->GetXaxis()->SetBinLabel(nStep,"Jet MaxMjj Category");
     nstepEvents[nStep-1]++;
     nStep = 12 ;

     if(!numberJetBinGen.empty() && !outputGenAK5_MjjSorted.empty()){
      if(numberJetBinGen.at(0)>=2 && outputGenAK5_MjjSorted.size()>=2){

      TLorentzVector Total4VMaxMjjGen = outputGenAK5_MjjSorted.at(0).Momentum_ + outputGenAK5_MjjSorted.at(1).Momentum_ ;
      NewMuonTree->vbf_maxMjj_jj_e_gen   = Total4VMaxMjjGen.E(); 
      NewMuonTree->vbf_maxMjj_jj_pt_gen  = Total4VMaxMjjGen.Pt(); 
      NewMuonTree->vbf_maxMjj_jj_eta_gen = Total4VMaxMjjGen.Eta(); 
      NewMuonTree->vbf_maxMjj_jj_phi_gen = Total4VMaxMjjGen.Phi(); 
      NewMuonTree->vbf_maxMjj_jj_m_gen   = Total4VMaxMjjGen.M(); 
      NewMuonTree->vbf_maxMjj_j1_e_gen   = outputGenAK5_MjjSorted.at(0).Momentum_.E();
      NewMuonTree->vbf_maxMjj_j1_pt_gen  = outputGenAK5_MjjSorted.at(0).Momentum_.Pt();
      NewMuonTree->vbf_maxMjj_j1_eta_gen = outputGenAK5_MjjSorted.at(0).Momentum_.Eta();
      NewMuonTree->vbf_maxMjj_j1_phi_gen = outputGenAK5_MjjSorted.at(0).Momentum_.Phi();
      NewMuonTree->vbf_maxMjj_j1_m_gen   = outputGenAK5_MjjSorted.at(0).Momentum_.M(); 
      NewMuonTree->vbf_maxMjj_j2_e_gen   = outputGenAK5_MjjSorted.at(1).Momentum_.E();
      NewMuonTree->vbf_maxMjj_j2_pt_gen  = outputGenAK5_MjjSorted.at(1).Momentum_.Pt();
      NewMuonTree->vbf_maxMjj_j2_eta_gen = outputGenAK5_MjjSorted.at(1).Momentum_.Eta();
      NewMuonTree->vbf_maxMjj_j2_phi_gen = outputGenAK5_MjjSorted.at(1).Momentum_.Phi();
      NewMuonTree->vbf_maxMjj_j2_m_gen   = outputGenAK5_MjjSorted.at(1).Momentum_.M();

      NewMuonTree->vbf_maxMjj_jj_deta_gen = fabs(outputGenAK5_MjjSorted.at(0).Momentum_.Eta() - outputGenAK5_MjjSorted.at(1).Momentum_.Eta()) ;
      if (fabs(outputGenAK5_MjjSorted.at(0).Momentum_.Phi() - outputGenAK5_MjjSorted.at(1).Momentum_.Phi()) < TMath::Pi())
       NewMuonTree->vbf_maxMjj_jj_dphi_gen = fabs(outputGenAK5_MjjSorted.at(0).Momentum_.Phi() - outputGenAK5_MjjSorted.at(1).Momentum_.Phi()) ;
      else 
       NewMuonTree->vbf_maxMjj_jj_dphi_gen = 2*TMath::Pi() - fabs(outputGenAK5_MjjSorted.at(0).Momentum_.Phi() - outputGenAK5_MjjSorted.at(1).Momentum_.Phi()) ;

      NewMuonTree->vbf_maxMjj_j1_bDiscriminatorSSVHE_gen = MuonTree->fReader->getFloat("JetGen_bDiscriminatorSSVHE")[outputGenAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_bDiscriminatorTCHE_gen  = MuonTree->fReader->getFloat("JetGen_bDiscriminatorTCHE")[outputGenAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_bDiscriminatorCSV_gen   = MuonTree->fReader->getFloat("JetGen_bDiscriminatorCSV")[outputGenAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_bDiscriminatorSSVHP_gen = MuonTree->fReader->getFloat("JetGen_bDiscriminatorSSVHP")[outputGenAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_bDiscriminatorTCHP_gen  = MuonTree->fReader->getFloat("JetGen_bDiscriminatorTCHP")[outputGenAK5_MjjSorted.at(0).position_] ;

      NewMuonTree->vbf_maxMjj_j2_bDiscriminatorSSVHE_gen = MuonTree->fReader->getFloat("JetGen_bDiscriminatorSSVHE")[outputGenAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_bDiscriminatorTCHE_gen  = MuonTree->fReader->getFloat("JetGen_bDiscriminatorTCHE")[outputGenAK5_MjjSorted.at(1).position_] ; 
      NewMuonTree->vbf_maxMjj_j2_bDiscriminatorCSV_gen   = MuonTree->fReader->getFloat("JetGen_bDiscriminatorCSV")[outputGenAK5_MjjSorted.at(1).position_]  ;
      NewMuonTree->vbf_maxMjj_j2_bDiscriminatorSSVHP_gen = MuonTree->fReader->getFloat("JetGen_bDiscriminatorSSVHP")[outputGenAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_bDiscriminatorTCHP_gen  = MuonTree->fReader->getFloat("JetGen_bDiscriminatorTCHP")[outputGenAK5_MjjSorted.at(1).position_] ;
      }
     }

     if(outputAK5_MjjSorted.at(0).NameCollection_ == "JetPFCor") {
    
      NewMuonTree->vbf_maxMjj_j1_QGLikelihood = MuonTree->fReader->getFloat("JetPFCor_QGLikelihood")[outputAK5_MjjSorted.at(0).position_] ;

      NewMuonTree->vbf_maxMjj_j1_isPileUpLoose  = MuonTree->fReader->getBool("JetPFCor_isPileUpJetLoose")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_isPileUpMedium = MuonTree->fReader->getBool("JetPFCor_isPileUpJetMedium")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_isPileUpTight  = MuonTree->fReader->getBool("JetPFCor_isPileUpJetTight")[outputAK5_MjjSorted.at(0).position_] ;


      NewMuonTree->vbf_maxMjj_j1_bDiscriminatorSSVHE = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorSSVHE")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_bDiscriminatorTCHE  = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorTCHE")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_bDiscriminatorCSV   = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorCSV")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_bDiscriminatorSSVHP = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorSSVHP")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_bDiscriminatorTCHP  = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorTCHP")[outputAK5_MjjSorted.at(0).position_] ;

      NewMuonTree->vbf_maxMjj_j1_ChargedHadronEnergy      = MuonTree->fReader->getFloat("JetPFCor_ChargedHadronEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_ChargedHadronEnergyFrac  = MuonTree->fReader->getFloat("JetPFCor_ChargedHadronEnergyFrac")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_NeutralHadronEnergy      = MuonTree->fReader->getFloat("JetPFCor_NeutralHadronEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_NeutralHadronEnergyFrac  = MuonTree->fReader->getFloat("JetPFCor_NeutralHadronEnergyFrac")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_ChargedEmEnergy          = MuonTree->fReader->getFloat("JetPFCor_ChargedEmEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_ChargedEmEnergyFrac      = MuonTree->fReader->getFloat("JetPFCor_ChargedEmEnergyFrac")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_ChargedMuEnergy          = MuonTree->fReader->getFloat("JetPFCor_ChargedMuEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_ChargedMuEnergyFrac      = MuonTree->fReader->getFloat("JetPFCor_ChargedMuEnergyFrac")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_NeutralEmEnergy          = MuonTree->fReader->getFloat("JetPFCor_NeutralEmEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_NeutralEmEnergyFrac      = MuonTree->fReader->getFloat("JetPFCor_NeutralEmEnergyFrac")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_PhotonEnergy             = MuonTree->fReader->getFloat("JetPFCor_PhotonEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_PhotonEnergyFraction     = MuonTree->fReader->getFloat("JetPFCor_PhotonEnergyFraction")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_ElectronEnergy           = MuonTree->fReader->getFloat("JetPFCor_ElectronEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_ElectronEnergyFraction   = MuonTree->fReader->getFloat("JetPFCor_ElectronEnergyFraction")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_HFHadronEnergy           = MuonTree->fReader->getFloat("JetPFCor_HFHadronEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_HFHadronEnergyFraction   = MuonTree->fReader->getFloat("JetPFCor_HFHadronEnergyFraction")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_HFEMEnergy               = MuonTree->fReader->getFloat("JetPFCor_HFEMEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_HFEMEnergyFraction       = MuonTree->fReader->getFloat("JetPFCor_HFEMEnergyFraction")[outputAK5_MjjSorted.at(0).position_] ;

      NewMuonTree->vbf_maxMjj_j1_ChargedMultiplicity        = MuonTree->fReader->getFloat("JetPFCor_ChargedMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_NeutralMultiplicity        = MuonTree->fReader->getFloat("JetPFCor_NeutralMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_MuonMultiplicity           = MuonTree->fReader->getFloat("JetPFCor_MuonMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_ChargedHadronMultiplicity  = MuonTree->fReader->getFloat("JetPFCor_ChargedHadronMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_NeutralHadronMultiplicity  = MuonTree->fReader->getFloat("JetPFCor_NeutralHadronMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_PhotonMultiplicity         = MuonTree->fReader->getFloat("JetPFCor_PhotonMultiplicity")[outputAK5_MjjSorted.at(0).position_] ; 
      NewMuonTree->vbf_maxMjj_j1_ElectronMultiplicity       = MuonTree->fReader->getFloat("JetPFCor_ElectronMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_HFHadronMultiplicity       = MuonTree->fReader->getFloat("JetPFCor_HFHadronMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
     }
     else if(outputAK5_MjjSorted.at(0).NameCollection_ == "JetPFCorVBFTag" ){

      NewMuonTree->vbf_maxMjj_j1_QGLikelihood = MuonTree->fReader->getFloat("JetPFCorVBFTag_QGLikelihood")[outputAK5_MjjSorted.at(0).position_] ;

      NewMuonTree->vbf_maxMjj_j1_isPileUpLoose  = MuonTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetLoose")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_isPileUpMedium = MuonTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetMedium")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_isPileUpTight  = MuonTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetTight")[outputAK5_MjjSorted.at(0).position_] ;


      NewMuonTree->vbf_maxMjj_j1_bDiscriminatorSSVHE = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorSSVHE")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_bDiscriminatorTCHE  = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorTCHE")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_bDiscriminatorCSV   = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorCSV")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_bDiscriminatorSSVHP = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorSSVHP")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_bDiscriminatorTCHP  = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorTCHP")[outputAK5_MjjSorted.at(0).position_] ;

      NewMuonTree->vbf_maxMjj_j1_ChargedHadronEnergy      = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_ChargedHadronEnergyFrac  = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronEnergyFrac")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_NeutralHadronEnergy      = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_NeutralHadronEnergyFrac  = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronEnergyFrac")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_ChargedEmEnergy          = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedEmEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_ChargedEmEnergyFrac      = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedEmEnergyFrac")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_ChargedMuEnergy          = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedMuEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_ChargedMuEnergyFrac      = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedMuEnergyFrac")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_NeutralEmEnergy          = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralEmEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_NeutralEmEnergyFrac  = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralEmEnergyFrac")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_PhotonEnergy         = MuonTree->fReader->getFloat("JetPFCorVBFTag_PhotonEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_PhotonEnergyFraction = MuonTree->fReader->getFloat("JetPFCorVBFTag_PhotonEnergyFraction")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_ElectronEnergy       = MuonTree->fReader->getFloat("JetPFCorVBFTag_ElectronEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_ElectronEnergyFraction   = MuonTree->fReader->getFloat("JetPFCorVBFTag_ElectronEnergyFraction")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_HFHadronEnergy           = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFHadronEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_HFHadronEnergyFraction   = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFHadronEnergyFraction")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_HFEMEnergy               = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFEMEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_HFEMEnergyFraction       = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFEMEnergyFraction")[outputAK5_MjjSorted.at(0).position_] ;

      NewMuonTree->vbf_maxMjj_j1_ChargedMultiplicity        = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_NeutralMultiplicity        = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_MuonMultiplicity           = MuonTree->fReader->getFloat("JetPFCorVBFTag_MuonMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_ChargedHadronMultiplicity  = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_NeutralHadronMultiplicity  = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_PhotonMultiplicity         = MuonTree->fReader->getFloat("JetPFCorVBFTag_PhotonMultiplicity")[outputAK5_MjjSorted.at(0).position_] ; 
      NewMuonTree->vbf_maxMjj_j1_ElectronMultiplicity       = MuonTree->fReader->getFloat("JetPFCorVBFTag_ElectronMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_HFHadronMultiplicity       = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFHadronMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;  

     }
     else { std::cerr<<" problem with High Deta Jet Name Collection "<<std::endl; continue ; }

     if(outputAK5_MjjSorted.at(1).NameCollection_ == "JetPFCor") {
    
      NewMuonTree->vbf_maxMjj_j2_QGLikelihood = MuonTree->fReader->getFloat("JetPFCor_QGLikelihood")[outputAK5_MjjSorted.at(1).position_] ;

      NewMuonTree->vbf_maxMjj_j2_isPileUpLoose = MuonTree->fReader->getBool("JetPFCor_isPileUpJetLoose")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_isPileUpMedium = MuonTree->fReader->getBool("JetPFCor_isPileUpJetMedium")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_isPileUpTight = MuonTree->fReader->getBool("JetPFCor_isPileUpJetTight")[outputAK5_MjjSorted.at(1).position_] ;


      NewMuonTree->vbf_maxMjj_j2_bDiscriminatorSSVHE = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorSSVHE")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_bDiscriminatorTCHE  = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorTCHE")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_bDiscriminatorCSV   = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorCSV")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_bDiscriminatorSSVHP = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorSSVHP")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_bDiscriminatorTCHP  = MuonTree->fReader->getFloat("JetPFCor_bDiscriminatorTCHP")[outputAK5_MjjSorted.at(1).position_] ;

      NewMuonTree->vbf_maxMjj_j2_ChargedHadronEnergy      = MuonTree->fReader->getFloat("JetPFCor_ChargedHadronEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_ChargedHadronEnergyFrac  = MuonTree->fReader->getFloat("JetPFCor_ChargedHadronEnergyFrac")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_NeutralHadronEnergy      = MuonTree->fReader->getFloat("JetPFCor_NeutralHadronEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_NeutralHadronEnergyFrac  = MuonTree->fReader->getFloat("JetPFCor_NeutralHadronEnergyFrac")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_ChargedEmEnergy          = MuonTree->fReader->getFloat("JetPFCor_ChargedEmEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_ChargedEmEnergyFrac      = MuonTree->fReader->getFloat("JetPFCor_ChargedEmEnergyFrac")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_ChargedMuEnergy          = MuonTree->fReader->getFloat("JetPFCor_ChargedMuEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_ChargedMuEnergyFrac      = MuonTree->fReader->getFloat("JetPFCor_ChargedMuEnergyFrac")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_NeutralEmEnergy          = MuonTree->fReader->getFloat("JetPFCor_NeutralEmEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_NeutralEmEnergyFrac  = MuonTree->fReader->getFloat("JetPFCor_NeutralEmEnergyFrac")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_PhotonEnergy         = MuonTree->fReader->getFloat("JetPFCor_PhotonEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_PhotonEnergyFraction = MuonTree->fReader->getFloat("JetPFCor_PhotonEnergyFraction")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_ElectronEnergy       = MuonTree->fReader->getFloat("JetPFCor_ElectronEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_ElectronEnergyFraction   = MuonTree->fReader->getFloat("JetPFCor_ElectronEnergyFraction")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_HFHadronEnergy           = MuonTree->fReader->getFloat("JetPFCor_HFHadronEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_HFHadronEnergyFraction   = MuonTree->fReader->getFloat("JetPFCor_HFHadronEnergyFraction")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_HFEMEnergy               = MuonTree->fReader->getFloat("JetPFCor_HFEMEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_HFEMEnergyFraction       = MuonTree->fReader->getFloat("JetPFCor_HFEMEnergyFraction")[outputAK5_MjjSorted.at(1).position_] ;

      NewMuonTree->vbf_maxMjj_j2_ChargedMultiplicity        = MuonTree->fReader->getFloat("JetPFCor_ChargedMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_NeutralMultiplicity        = MuonTree->fReader->getFloat("JetPFCor_NeutralMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_MuonMultiplicity           = MuonTree->fReader->getFloat("JetPFCor_MuonMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_ChargedHadronMultiplicity  = MuonTree->fReader->getFloat("JetPFCor_ChargedHadronMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_NeutralHadronMultiplicity  = MuonTree->fReader->getFloat("JetPFCor_NeutralHadronMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_PhotonMultiplicity         = MuonTree->fReader->getFloat("JetPFCor_PhotonMultiplicity")[outputAK5_MjjSorted.at(1).position_] ; 
      NewMuonTree->vbf_maxMjj_j2_ElectronMultiplicity       = MuonTree->fReader->getFloat("JetPFCor_ElectronMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_HFHadronMultiplicity       = MuonTree->fReader->getFloat("JetPFCor_HFHadronMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
     }
     else if(outputAK5_MjjSorted.at(1).NameCollection_ == "JetPFCorVBFTag" ){

      NewMuonTree->vbf_maxMjj_j2_QGLikelihood = MuonTree->fReader->getFloat("JetPFCorVBFTag_QGLikelihood")[outputAK5_MjjSorted.at(1).position_] ;

      NewMuonTree->vbf_maxMjj_j2_isPileUpLoose  = MuonTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetLoose")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_isPileUpMedium = MuonTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetMedium")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_isPileUpTight  = MuonTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetTight")[outputAK5_MjjSorted.at(1).position_] ;


      NewMuonTree->vbf_maxMjj_j2_bDiscriminatorSSVHE = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorSSVHE")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_bDiscriminatorTCHE  = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorTCHE")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_bDiscriminatorCSV   = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorCSV")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_bDiscriminatorSSVHP = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorSSVHP")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_bDiscriminatorTCHP  = MuonTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorTCHP")[outputAK5_MjjSorted.at(1).position_] ;

      NewMuonTree->vbf_maxMjj_j2_ChargedHadronEnergy      = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_ChargedHadronEnergyFrac  = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronEnergyFrac")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_NeutralHadronEnergy      = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_NeutralHadronEnergyFrac  = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronEnergyFrac")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_ChargedEmEnergy          = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedEmEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_ChargedEmEnergyFrac      = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedEmEnergyFrac")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_ChargedMuEnergy          = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedMuEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_ChargedMuEnergyFrac      = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedMuEnergyFrac")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_NeutralEmEnergy          = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralEmEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_NeutralEmEnergyFrac  = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralEmEnergyFrac")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_PhotonEnergy         = MuonTree->fReader->getFloat("JetPFCorVBFTag_PhotonEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_PhotonEnergyFraction = MuonTree->fReader->getFloat("JetPFCorVBFTag_PhotonEnergyFraction")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_ElectronEnergy       = MuonTree->fReader->getFloat("JetPFCorVBFTag_ElectronEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_ElectronEnergyFraction   = MuonTree->fReader->getFloat("JetPFCorVBFTag_ElectronEnergyFraction")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_HFHadronEnergy           = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFHadronEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_HFHadronEnergyFraction   = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFHadronEnergyFraction")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_HFEMEnergy               = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFEMEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_HFEMEnergyFraction       = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFEMEnergyFraction")[outputAK5_MjjSorted.at(1).position_] ;

      NewMuonTree->vbf_maxMjj_j2_ChargedMultiplicity  = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_NeutralMultiplicity  = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_MuonMultiplicity     = MuonTree->fReader->getFloat("JetPFCorVBFTag_MuonMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_ChargedHadronMultiplicity  = MuonTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_NeutralHadronMultiplicity  = MuonTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_PhotonMultiplicity         = MuonTree->fReader->getFloat("JetPFCorVBFTag_PhotonMultiplicity")[outputAK5_MjjSorted.at(1).position_] ; 
      NewMuonTree->vbf_maxMjj_j2_ElectronMultiplicity       = MuonTree->fReader->getFloat("JetPFCorVBFTag_ElectronMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_HFHadronMultiplicity       = MuonTree->fReader->getFloat("JetPFCorVBFTag_HFHadronMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;  
    
     }
     else { std::cerr<<" problem with High Deta Jet Name Collection "<<std::endl; continue ; }
          
    }
    NewMuonTree->fTree->Fill(); // Fill the events 
    
   } // End of Loop on the event

   // Save Results in the output
   NewMuonTree->fTree->Write(TreeName.c_str());

   int ibin = 0 ;
   nStep = 0 ;
   for( ; nstepEvents[ibin]!=0 ; ibin ++){ nStep ++ ; SelectionEvents->SetBinContent(ibin+1,nstepEvents[ibin]); }

   for(ibin =1 ; ibin<SelectionEvents->GetNbinsX() ; ibin ++) {
    if(SelectionEvents->GetBinContent(ibin+1)!=0) SelectionEfficiency->SetBinContent(ibin+1,SelectionEvents->GetBinContent(ibin+1)/SelectionEvents->GetBinContent(ibin)); 
   }

   SelectionEvents->GetXaxis()->SetRangeUser(0,nStep-1);
   SelectionEvents->Write();

   SelectionEfficiency->GetXaxis()->SetRangeUser(1,nStep-1);
   SelectionEfficiency->Write();

   std::cout << " Finish :: " << outputFile->GetName() << "    "<<  MuonTree->fTree->GetEntries ()  << std::endl;
   outputFile->Close();
  } // End of Muon Analysis
  
  // Electron Sample Processing 
  if(LeptonType == "Electron"){

   std::cout<<" Enter in the Electron Category "<<std::endl;
   std::cout<<"                            "<<std::endl;

   // create and open the output file 

   TFile *outputFile = new TFile((OutputRootDirectory+"/"+OutputRootFile).c_str(),"RECREATE");
   outputFile->cd();

   std::cout<<" Open Input File : "<<inputFile->GetName()<<" TreeName  "<<TreeName<<std::endl;

   VBFElectronClass* ElectronTree = new VBFElectronClass (inputFile,TreeName);
   
   ElectronTree->SetReader(ElectronTree->fTree);

   std::cout<<"                                   "<<std::endl;
   std::cout<<" Clone Tree  "<<std::endl;
   
   TTree *newtree = ElectronTree->fTree->CloneTree(0);

   // Add new Branches 
   VBFElectronClass* NewElectronTree = new VBFElectronClass(newtree);
   
   NewElectronTree->SetNewBranches(NewElectronTree->fTree);
  
   // Loop on the events 

   std::cout<<"                                   "<<std::endl;
   std::cout << "Input Tree Number of Entries : " <<  ElectronTree->fTree->GetEntries ()  << std::endl ;
   std::cout<<"                                   "<<std::endl;

   for(int iEntry = 0 ; iEntry <  ElectronTree->fTree->GetEntries ()  ; iEntry++){

    nStep = 1 ;
    if(std::string(SelectionEvents->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEvents->GetXaxis()->SetBinLabel(nStep,"All Events");
    if(std::string(SelectionEfficiency->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEfficiency->GetXaxis()->SetBinLabel(nStep,"All Events");
    nstepEvents[nStep-1]++;
    nStep = 2;
    
    ElectronTree->fTree->GetEntry(iEntry); 

    NewElectronTree->InitializateVariables(); 

    if (iEntry % 10000 == 0) std::cout << "reading event " << iEntry << std::endl ;

    // Basic Selections for boosted region 
    if(!(ElectronTree->fReader->getInt("isgengdboostedWevt")[0]) || (ElectronTree->fReader->getFloat("GroomedJet_CA8_deltaR_lca8jet")[0]) < TMath::Pi()/ 2.0  || 
        (ElectronTree->fReader->getFloat("GroomedJet_CA8_pt")[0])< JetPtWboostedMin || (ElectronTree->fReader->getFloat("W_pt")[0]) < JetPtWboostedMin ) continue ;
    
    if(std::string(SelectionEvents->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEvents->GetXaxis()->SetBinLabel(nStep,"Base Boosted W");    
    if(std::string(SelectionEfficiency->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEfficiency->GetXaxis()->SetBinLabel(nStep,"Base Boosted W");
    nstepEvents[nStep-1]++;
    nStep = 3;

    // Join forward and central PFJetCor Collection    
    std::vector<std::vector<TLorentzVector> > GroomedJet_CA8_Collection (JetPtCutMin.size());

    std::vector< std::vector <JetAK5> > JetPFCor_AK5_Collection (JetPtCutMin.size()) ;
    std::vector< std::vector <JetAK5> > CleanedJetPFCor_AK5_Collection(JetPtCutMin.size()) ;
    std::vector< std::vector <JetAK5> > HadronicW_AK5_Collection(JetPtCutMin.size()) ;

    std::vector<std::vector<TLorentzVector> > GenGroomedJet_CA8_Collection(JetPtCutMin.size()); 
    std::vector< std::vector <JetAK5> > GenJetPFCor_AK5_Collection(JetPtCutMin.size()); 
    std::vector< std::vector <JetAK5> > GenCleanedJetPFCor_AK5_Collection(JetPtCutMin.size());
    std::vector< std::vector <JetAK5> > GenHadronicW_AK5_Collection(JetPtCutMin.size());
   
    for(int iJetPtCutMin = 0 ; iJetPtCutMin < int(JetPtCutMin.size()); iJetPtCutMin++){

     for(int iJet = 0 ; iJet < JetCollectionDimension ; iJet++) { // run over the whole CA8 jet collection and took the 4vetcor 
      TLorentzVector JetTemp , GenJetTemp ;
      std::string nameCollection ; 
      JetTemp.SetPtEtaPhiE(ElectronTree->fReader->getFloat("GroomedJet_CA8_pt")[iJet],ElectronTree->fReader->getFloat("GroomedJet_CA8_eta")[iJet], 
                           ElectronTree->fReader->getFloat("GroomedJet_CA8_phi")[iJet],ElectronTree->fReader->getFloat("GroomedJet_CA8_e")[iJet]);

      if(NewElectronTree->fTree->FindBranch("GenGroomedJet_CA8_pt") && NewElectronTree->fTree->FindBranch("JetGen_Pt")){
       GenJetTemp.SetPtEtaPhiE(ElectronTree->fReader->getFloat("GenGroomedJet_CA8_pt")[iJet],ElectronTree->fReader->getFloat("GenGroomedJet_CA8_eta")[iJet], 
                              ElectronTree->fReader->getFloat("GenGroomedJet_CA8_phi")[iJet],ElectronTree->fReader->getFloat("GenGroomedJet_CA8_e")[iJet]);
      }

     // Selection on CA8 Jets -> pt cut on each jet over the threshold and acceptance
      if(fabs(JetTemp.Eta())<JetEtaCutMax && JetTemp.Pt()>JetPtCutMin.at(iJetPtCutMin))
		GroomedJet_CA8_Collection.at(iJetPtCutMin).push_back(JetTemp);
               

      if(NewElectronTree->fTree->FindBranch("GenGroomedJet_CA8_pt") && NewElectronTree->fTree->FindBranch("JetGen_Pt")){
       if(fabs(GenJetTemp.Eta())<JetEtaCutMax && GenJetTemp.Pt()>JetPtCutMin.at(iJetPtCutMin))
		GenGroomedJet_CA8_Collection.at(iJetPtCutMin).push_back(GenJetTemp);
      }
      // take the central AK5 jets 
      JetTemp.SetPtEtaPhiE(ElectronTree->fReader->getFloat("JetPFCor_Pt")[iJet],ElectronTree->fReader->getFloat("JetPFCor_Eta")[iJet],
 			  ElectronTree->fReader->getFloat("JetPFCor_Phi")[iJet],ElectronTree->fReader->getFloat("JetPFCor_E")[iJet]);

      if(NewElectronTree->fTree->FindBranch("GenGroomedJet_CA8_pt") && NewElectronTree->fTree->FindBranch("JetGen_Pt"))
       GenJetTemp.SetPtEtaPhiE(ElectronTree->fReader->getFloat("JetGen_Pt")[iJet],ElectronTree->fReader->getFloat("JetGen_Eta")[iJet],
   	   		       ElectronTree->fReader->getFloat("JetGen_Phi")[iJet],ElectronTree->fReader->getFloat("JetGen_E")[iJet]);
      
     // Selection on PF Cor Central Jets --> AK5
     if(fabs(JetTemp.Eta())<JetEtaCutMax && JetTemp.Pt()>JetPtCutMin.at(iJetPtCutMin)){
       JetAK5 tempJetAK5 (iJet,"JetPFCor",JetTemp);
       JetPFCor_AK5_Collection.at(iJetPtCutMin).push_back(tempJetAK5);
     }

     if(NewElectronTree->fTree->FindBranch("GenGroomedJet_CA8_pt") && NewElectronTree->fTree->FindBranch("JetGen_Pt")){
      if(fabs(GenJetTemp.Eta())<JetEtaCutMax && GenJetTemp.Pt()>JetPtCutMin.at(iJetPtCutMin)){
       JetAK5 tempJetAK5 (iJet,"GenJet",GenJetTemp);
       GenJetPFCor_AK5_Collection.at(iJetPtCutMin).push_back(tempJetAK5);
      }
      }
    }
   
    for(int iJet = 0 ; iJet < JetCollectionDimension ; iJet++) { //only AK5 forward jet over the pt threshold 
     TLorentzVector JetTemp ; 
     JetTemp.SetPtEtaPhiE(ElectronTree->fReader->getFloat("JetPFCorVBFTag_Pt")[iJet],ElectronTree->fReader->getFloat("JetPFCorVBFTag_Eta")[iJet], 
                          ElectronTree->fReader->getFloat("JetPFCorVBFTag_Phi")[iJet],ElectronTree->fReader->getFloat("JetPFCorVBFTag_E")[iJet]);

     // Selection on PF Cor Forward Jets --> AK5
     if(fabs(JetTemp.Eta())<JetEtaCutMax && JetTemp.Pt()>JetPtCutMin.at(iJetPtCutMin)){
       JetAK5 tempJetAK5 (iJet,"JetPFCorVBFTag",JetTemp);
       JetPFCor_AK5_Collection.at(iJetPtCutMin).push_back(tempJetAK5);
     }
    }
   
    //choose the jet corresponding to the hadronic W and fill new branches with its variables --> apply the pT cut on the whole CA8 jet and select the one with mass closer
    // to the W mass just to have another solution 
    float difference = 1000.;
    int iWHadronic = 0;

    for(size_t i = 0; i < GroomedJet_CA8_Collection.size() ; i ++){
      if ( ElectronTree->fReader->getFloat("GroomedJet_CA8_pt")[i]>JetPtWboostedMin ) {
	if ( fabs (ElectronTree->fReader->getFloat("GroomedJet_CA8_mass_pr")[i] - Wmass) < difference ) {
	  difference = fabs (ElectronTree->fReader->getFloat("GroomedJet_CA8_mass_pr")[i] - Wmass);
	  iWHadronic = i;
	}
	}
    }
    if(iJetPtCutMin == 0){

     NewElectronTree -> WHadposition = iWHadronic;   //position of the hadronic W in the CA8Jet collection with pt over threshold + mass closer to the W one

     NewElectronTree->Hadronic_W_Jet_mass_uncorr    = ElectronTree->fReader->getFloat("GroomedJet_CA8_mass_uncorr")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_mass_tr_uncorr = ElectronTree->fReader->getFloat("GroomedJet_CA8_mass_tr_uncorr")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_mass_ft_uncorr = ElectronTree->fReader->getFloat("GroomedJet_CA8_mass_ft_uncorr")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_mass_pr_uncorr = ElectronTree->fReader->getFloat("GroomedJet_CA8_mass_pr_uncorr")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_massdrop_pr_uncorr = ElectronTree->fReader->getFloat("GroomedJet_CA8_massdrop_pr_uncorr")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_tau2tau1           = ElectronTree->fReader->getFloat("GroomedJet_CA8_tau2tau1")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_tau1 = ElectronTree->fReader->getFloat("GroomedJet_CA8_tau1")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_tau2 = ElectronTree->fReader->getFloat("GroomedJet_CA8_tau2")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_tau3 = ElectronTree->fReader->getFloat("GroomedJet_CA8_tau3")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_tau4 = ElectronTree->fReader->getFloat("GroomedJet_CA8_tau4")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_pt   = ElectronTree->fReader->getFloat("GroomedJet_CA8_pt")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_eta  = ElectronTree->fReader->getFloat("GroomedJet_CA8_eta")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_phi  = ElectronTree->fReader->getFloat("GroomedJet_CA8_phi")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_e    = ElectronTree->fReader->getFloat("GroomedJet_CA8_e")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_pt_tr_uncorr = ElectronTree->fReader->getFloat("GroomedJet_CA8_pt_tr_uncorr")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_pt_tr  = ElectronTree->fReader->getFloat("GroomedJet_CA8_pt_tr")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_eta_tr = ElectronTree->fReader->getFloat("GroomedJet_CA8_eta_tr")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_phi_tr = ElectronTree->fReader->getFloat("GroomedJet_CA8_phi_tr")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_e_tr   = ElectronTree->fReader->getFloat("GroomedJet_CA8_e_tr")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_pt_ft_uncorr = ElectronTree->fReader->getFloat("GroomedJet_CA8_pt_ft_uncorr")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_pt_ft  = ElectronTree->fReader->getFloat("GroomedJet_CA8_pt_ft")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_eta_ft = ElectronTree->fReader->getFloat("GroomedJet_CA8_eta_ft")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_phi_ft = ElectronTree->fReader->getFloat("GroomedJet_CA8_phi_ft")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_e_ft   = ElectronTree->fReader->getFloat("GroomedJet_CA8_e_ft")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_pt_pr_uncorr = ElectronTree->fReader->getFloat("GroomedJet_CA8_pt_pr_uncorr")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_pt_pr  = ElectronTree->fReader->getFloat("GroomedJet_CA8_pt_tr")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_eta_pr = ElectronTree->fReader->getFloat("GroomedJet_CA8_eta_pr")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_phi_pr = ElectronTree->fReader->getFloat("GroomedJet_CA8_phi_pr")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_e_pr   = ElectronTree->fReader->getFloat("GroomedJet_CA8_e_pr")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_prsubjet1_px = ElectronTree->fReader->getFloat("GroomedJet_CA8_prsubjet1_px")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_prsubjet1_py = ElectronTree->fReader->getFloat("GroomedJet_CA8_prsubjet1_py")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_prsubjet1_pz = ElectronTree->fReader->getFloat("GroomedJet_CA8_prsubjet1_pz")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_prsubjet1_e  = ElectronTree->fReader->getFloat("GroomedJet_CA8_prsubjet1_e")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_prsubjet2_px = ElectronTree->fReader->getFloat("GroomedJet_CA8_prsubjet2_px")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_prsubjet2_py = ElectronTree->fReader->getFloat("GroomedJet_CA8_prsubjet2_py")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_prsubjet2_pz = ElectronTree->fReader->getFloat("GroomedJet_CA8_prsubjet2_pz")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_prsubjet2_e  = ElectronTree->fReader->getFloat("GroomedJet_CA8_prsubjet2_e")[iWHadronic];  
     NewElectronTree->Hadronic_W_Jet_mass    = ElectronTree->fReader->getFloat("GroomedJet_CA8_mass")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_mass_tr = ElectronTree->fReader->getFloat("GroomedJet_CA8_mass_tr")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_mass_ft = ElectronTree->fReader->getFloat("GroomedJet_CA8_mass_ft")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_mass_pr = ElectronTree->fReader->getFloat("GroomedJet_CA8_mass_pr")[iWHadronic];  
     NewElectronTree->Hadronic_W_Jet_massdrop = ElectronTree->fReader->getFloat("GroomedJet_CA8_massdrop_pr")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_area    = ElectronTree->fReader->getFloat("GroomedJet_CA8_area")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_area_tr = ElectronTree->fReader->getFloat("GroomedJet_CA8_area_tr")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_area_ft = ElectronTree->fReader->getFloat("GroomedJet_CA8_area_ft")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_area_pr = ElectronTree->fReader->getFloat("GroomedJet_CA8_area_pr")[iWHadronic]; 
     NewElectronTree->Hadronic_W_Jet_jetconsituents = ElectronTree->fReader->getFloat("GroomedJet_CA8_jetconstituents")[iWHadronic]; 
     NewElectronTree->Hadronic_W_Jet_jetcharge = ElectronTree->fReader->getFloat("GroomedJet_CA8_jetcharge")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_rcores = ElectronTree->fReader->getFloat("GroomedJet_CA8_rcores")[iWHadronic];  
     NewElectronTree->Hadronic_W_Jet_ptcores = ElectronTree->fReader->getFloat("GroomedJet_CA8_ptcores")[iWHadronic];  
     NewElectronTree->Hadronic_W_Jet_planarflow = ElectronTree->fReader->getFloat("GroomedJet_CA8_planarflow")[iWHadronic];  
     NewElectronTree->Hadronic_W_Jet_qjetmass = ElectronTree->fReader->getFloat("GroomedJet_CA8_qjetmass")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_qjetmassdrop = ElectronTree->fReader->getFloat("GroomedJet_CA8_qjetmassdrop")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_deltaR_ljet = ElectronTree->fReader->getFloat("GroomedJet_CA8_deltaR_lca8jet")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_deltaphi_METjet = ElectronTree->fReader->getFloat("GroomedJet_CA8_deltaphi_METca8jet")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_deltaphi_Vca8jet = ElectronTree->fReader->getFloat("GroomedJet_CA8_deltaphi_Vca8jet")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_rcores01 = ElectronTree->fReader->getFloat("GroomedJet_CA8_rcores01")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_rcores02 = ElectronTree->fReader->getFloat("GroomedJet_CA8_rcores02")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_rcores03 = ElectronTree->fReader->getFloat("GroomedJet_CA8_rcores03")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_rcores04 = ElectronTree->fReader->getFloat("GroomedJet_CA8_rcores04")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_rcores05 = ElectronTree->fReader->getFloat("GroomedJet_CA8_rcores05")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_rcores06 = ElectronTree->fReader->getFloat("GroomedJet_CA8_rcores06")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_rcores07 = ElectronTree->fReader->getFloat("GroomedJet_CA8_rcores07")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_rcores08 = ElectronTree->fReader->getFloat("GroomedJet_CA8_rcores08")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_rcores09 = ElectronTree->fReader->getFloat("GroomedJet_CA8_rcores09")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_rcores10 = ElectronTree->fReader->getFloat("GroomedJet_CA8_rcores10")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_rcores11 = ElectronTree->fReader->getFloat("GroomedJet_CA8_rcores11")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_ptcores01 = ElectronTree->fReader->getFloat("GroomedJet_CA8_ptcores01")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_ptcores02 = ElectronTree->fReader->getFloat("GroomedJet_CA8_ptcores02")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_ptcores03 = ElectronTree->fReader->getFloat("GroomedJet_CA8_ptcores03")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_ptcores04 = ElectronTree->fReader->getFloat("GroomedJet_CA8_ptcores04")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_ptcores05 = ElectronTree->fReader->getFloat("GroomedJet_CA8_ptcores05")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_ptcores06 = ElectronTree->fReader->getFloat("GroomedJet_CA8_ptcores06")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_ptcores07 = ElectronTree->fReader->getFloat("GroomedJet_CA8_ptcores07")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_ptcores08 = ElectronTree->fReader->getFloat("GroomedJet_CA8_ptcores08")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_ptcores09 = ElectronTree->fReader->getFloat("GroomedJet_CA8_ptcores09")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_ptcores10 = ElectronTree->fReader->getFloat("GroomedJet_CA8_ptcores10")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_ptcores11 = ElectronTree->fReader->getFloat("GroomedJet_CA8_ptcores11")[iWHadronic]; 
     NewElectronTree->Hadronic_W_Jet_planarflow01 = ElectronTree->fReader->getFloat("GroomedJet_CA8_planarflow01")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_planarflow02 = ElectronTree->fReader->getFloat("GroomedJet_CA8_planarflow02")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_planarflow03 = ElectronTree->fReader->getFloat("GroomedJet_CA8_planarflow03")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_planarflow04 = ElectronTree->fReader->getFloat("GroomedJet_CA8_planarflow04")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_planarflow05 = ElectronTree->fReader->getFloat("GroomedJet_CA8_planarflow05")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_planarflow06 = ElectronTree->fReader->getFloat("GroomedJet_CA8_planarflow06")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_planarflow07 = ElectronTree->fReader->getFloat("GroomedJet_CA8_planarflow07")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_planarflow08 = ElectronTree->fReader->getFloat("GroomedJet_CA8_planarflow08")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_planarflow09 = ElectronTree->fReader->getFloat("GroomedJet_CA8_planarflow09")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_planarflow10 = ElectronTree->fReader->getFloat("GroomedJet_CA8_planarflow10")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_planarflow11 = ElectronTree->fReader->getFloat("GroomedJet_CA8_planarflow11")[iWHadronic];
     NewElectronTree->Hadronic_W_Jet_mass_sensi_tr = ElectronTree->fReader->getFloat("GroomedJet_CA8_mass_sensi_tr")[iWHadronic]; 
     NewElectronTree->Hadronic_W_Jet_mass_sensi_ft = ElectronTree->fReader->getFloat("GroomedJet_CA8_mass_sensi_ft")[iWHadronic];     
     NewElectronTree->Hadronic_W_Jet_mass_sensi_pr = ElectronTree->fReader->getFloat("GroomedJet_CA8_mass_sensi_pr")[iWHadronic]; 
     NewElectronTree->Hadronic_W_Jet_qjetmassvolatility = ElectronTree->fReader->getFloat("GroomedJet_CA8_qjetmassvolatility")[iWHadronic]; 
     NewElectronTree->Hadronic_W_Jet_prsubjet1ptoverjetpt = ElectronTree->fReader->getFloat("GroomedJet_CA8_prsubjet1ptoverjetpt")[iWHadronic]; 
     NewElectronTree->Hadronic_W_Jet_prsubjet2ptoverjetpt = ElectronTree->fReader->getFloat("GroomedJet_CA8_prsubjet2ptoverjetpt")[iWHadronic]; 
     NewElectronTree->Hadronic_W_Jet_prsubjet1subjet2_deltaR = ElectronTree->fReader->getFloat("GroomedJet_CA8_prsubjet1subjet2_deltaR")[iWHadronic]; 

     }
   }
    // Calculate Neutrino Pz using all the possible choices : type0 -> if real roots, pick the one nearest to the lepton Pz except when the Pz so chosen 
    //                                                                 is greater than 300 GeV in which case pick the most central root.               
    //                                                        type1 -> type = 1: if real roots, choose the one closest to the lepton Pz                                                                                                                                    if complex roots, use only the real part.     
    //                                                        type = 2: if real roots, choose the most central solution.                                                                                                                                          if complex roots, use only the real part.                                                                                                                                       type = 3: if real roots, pick the largest value of the cosine*                         

    TLorentzVector W_electron, W_Met;
   
    W_electron.SetPxPyPzE(ElectronTree->fReader->getFloat("W_electron_px")[0],ElectronTree->fReader->getFloat("W_electron_py")[0],
                          ElectronTree->fReader->getFloat("W_electron_pz")[0],ElectronTree->fReader->getFloat("W_electron_e")[0]);
    W_Met.SetPxPyPzE(ElectronTree->fReader->getFloat("event_met_pfmet")[0] * TMath::Cos(ElectronTree->fReader->getFloat("event_met_pfmetPhi")[0]),
                     ElectronTree->fReader->getFloat("event_met_pfmet")[0] * TMath::Sin(ElectronTree->fReader->getFloat("event_met_pfmetPhi")[0]),0.,
                     fabs(ElectronTree->fReader->getFloat("event_met_pfmet")[0]));

    if(W_electron.Pt()<=0 || W_Met.Pt() <= 0 ){ std::cerr<<" Negative Lepton - Neutrino Pt "<<std::endl; continue ; }

    if(std::string(SelectionEvents->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEvents->GetXaxis()->SetBinLabel(nStep,"Leptonic 4V");
    if(std::string(SelectionEfficiency->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEfficiency->GetXaxis()->SetBinLabel(nStep,"Leptonic 4V");
    nstepEvents[nStep-1]++;
    nStep = 4;

    // type0 calculation of neutrino pZ
    METzCalculator<TLorentzVector> NeutrinoPz_type0;
    NeutrinoPz_type0.SetMET(W_Met);
    NeutrinoPz_type0.SetLepton(W_electron);
    NeutrinoPz_type0.SetLeptonType("muon");

    double pz1_type0 = NeutrinoPz_type0.Calculate(); // Default one -> according to type0
    double pz2_type0 = NeutrinoPz_type0.getOther(); // Default one

    // don't touch the neutrino pT
    TLorentzVector W_neutrino_type0_met; W_neutrino_type0_met.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type0,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type0*pz1_type0));
    NewElectronTree->W_mass_type0_met   = (W_neutrino_type0_met+W_electron).M();  
    NewElectronTree->W_pz_type0_met     = (W_neutrino_type0_met+W_electron).Pz();   
    NewElectronTree->W_nu1_pz_type0_met = pz1_type0; 
    NewElectronTree->W_nu2_pz_type0_met = pz2_type0;

    // chenge the neutrino pT in case of complex solution in order to make it real
    TLorentzVector W_neutrino_type0;  W_neutrino_type0.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type0,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type0*pz1_type0));

    if (NeutrinoPz_type0.IsComplex()) {// if this is a complix, change MET                                                                                                                            
     double nu_pt1 = NeutrinoPz_type0.getPtneutrino(1);
     double nu_pt2 = NeutrinoPz_type0.getPtneutrino(2);
   
     TLorentzVector W_neutrino_1;
     W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(ElectronTree->fReader->getFloat("event_met_pfmetPhi")[0]), 
                             nu_pt1 * TMath::Sin(ElectronTree->fReader->getFloat("event_met_pfmetPhi")[0]), pz1_type0, sqrt(nu_pt1*nu_pt1 + pz1_type0*pz1_type0) );
     TLorentzVector W_neutrino_2;
     W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(ElectronTree->fReader->getFloat("event_met_pfmetPhi")[0]),
                             nu_pt2 * TMath::Sin(ElectronTree->fReader->getFloat("event_met_pfmetPhi")[0]), pz1_type0, sqrt(nu_pt2*nu_pt2 + pz1_type0*pz1_type0) );

     if ( fabs((W_electron+W_neutrino_1).M()-Wmass) < fabs((W_electron+W_neutrino_2).M()-Wmass) )  W_neutrino_type0 = W_neutrino_1;
     else W_neutrino_type0 = W_neutrino_2;

    }

    NewElectronTree->W_mass_type0 = (W_electron+W_neutrino_type0).M();  
    NewElectronTree->W_pz_type0   = (W_electron+W_neutrino_type0).Pz();  
    NewElectronTree->W_nu1_pz_type0 = pz1_type0;  
    NewElectronTree->W_nu2_pz_type0 = pz2_type0;

  

    // type2 calculation of neutrino pZ
    METzCalculator<TLorentzVector> NeutrinoPz_type2;
    NeutrinoPz_type2.SetMET(W_Met);
    NeutrinoPz_type2.SetLepton(W_electron);
    NeutrinoPz_type2.SetLeptonType("muon");
    double pz1_type2 = NeutrinoPz_type2.Calculate(2); // Default one -> according to type2
    double pz2_type2 = NeutrinoPz_type2.getOther(); // Default one

    // don't touch the neutrino pT
    TLorentzVector W_neutrino_type2_met; W_neutrino_type2_met.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type2,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type2*pz1_type2));
    NewElectronTree->W_mass_type2_met   = (W_neutrino_type2_met+W_electron).M();  
    NewElectronTree->W_pz_type2_met     = (W_neutrino_type2_met+W_electron).Pz();   
    NewElectronTree->W_nu1_pz_type2_met = pz1_type2; 
    NewElectronTree->W_nu2_pz_type2_met = pz2_type2;

    // chenge the neutrino pT in case of complex solution in order to make it real
    TLorentzVector W_neutrino_type2;  W_neutrino_type2.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type2,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type2*pz1_type2));

    if (NeutrinoPz_type2.IsComplex()) {// if this is a complix, change MET                                                                                                                            
     double nu_pt1 = NeutrinoPz_type2.getPtneutrino(1);
     double nu_pt2 = NeutrinoPz_type2.getPtneutrino(2);
   
     TLorentzVector W_neutrino_1;
     W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(ElectronTree->fReader->getFloat("event_met_pfmetPhi")[0]), 
                             nu_pt1 * TMath::Sin(ElectronTree->fReader->getFloat("event_met_pfmetPhi")[0]), pz1_type2, sqrt(nu_pt1*nu_pt1 + pz1_type2*pz1_type2) );
     TLorentzVector W_neutrino_2;
     W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(ElectronTree->fReader->getFloat("event_met_pfmetPhi")[0]),
                             nu_pt2 * TMath::Sin(ElectronTree->fReader->getFloat("event_met_pfmetPhi")[0]), pz1_type2, sqrt(nu_pt2*nu_pt2 + pz1_type2*pz1_type2) );

     if ( fabs((W_electron+W_neutrino_1).M()-Wmass) < fabs((W_electron+W_neutrino_2).M()-Wmass) )  W_neutrino_type2 = W_neutrino_1;
     else W_neutrino_type2 = W_neutrino_2;

    }

    NewElectronTree->W_mass_type2 = (W_electron+W_neutrino_type2).M();  
    NewElectronTree->W_pz_type2   = (W_electron+W_neutrino_type2).Pz();  
    NewElectronTree->W_nu1_pz_type2 = pz1_type2;  
    NewElectronTree->W_nu2_pz_type2 = pz2_type2;


    //////////////////////////////
    
    TLorentzVector W_subjet1, W_subjet2 ;  // take the two subjet of the hardest CA8 after pruning
   
    W_subjet1.SetPxPyPzE(ElectronTree->fReader->getFloat("GroomedJet_CA8_prsubjet1_px")[0],ElectronTree->fReader->getFloat("GroomedJet_CA8_prsubjet1_py")[0],
         		 ElectronTree->fReader->getFloat("GroomedJet_CA8_prsubjet1_pz")[0],ElectronTree->fReader->getFloat("GroomedJet_CA8_prsubjet1_e")[0] );
    W_subjet2.SetPxPyPzE(ElectronTree->fReader->getFloat("GroomedJet_CA8_prsubjet2_px")[0],ElectronTree->fReader->getFloat("GroomedJet_CA8_prsubjet2_py")[0],
  		  	 ElectronTree->fReader->getFloat("GroomedJet_CA8_prsubjet2_pz")[0],ElectronTree->fReader->getFloat("GroomedJet_CA8_prsubjet2_e")[0] );

    if(W_subjet1.Pt() <= 0 || W_subjet2.Pt() <= 0){ std::cerr<<" Problem with subjets "<<std::endl; continue ; }

    if(std::string(SelectionEvents->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEvents->GetXaxis()->SetBinLabel(nStep,"Subjets 4V");
    if(std::string(SelectionEfficiency->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEfficiency->GetXaxis()->SetBinLabel(nStep,"Subjets 4V");
    nstepEvents[nStep-1]++;
    nStep = 5;


    TLorentzVector W_GroomedJet_CA8; 
    W_GroomedJet_CA8.SetPtEtaPhiE(ElectronTree->fReader->getFloat("GroomedJet_CA8_pt")[0], ElectronTree->fReader->getFloat("GroomedJet_CA8_eta")[0],
                                     ElectronTree->fReader->getFloat("GroomedJet_CA8_phi")[0], ElectronTree->fReader->getFloat("GroomedJet_CA8_e")[0]);

    if(W_GroomedJet_CA8.Pt() <=0){ std::cerr<<" Problem with pruned CA8 "<<std::endl; continue ;}

    if(std::string(SelectionEvents->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEvents->GetXaxis()->SetBinLabel(nStep,"Pruned CA8 4V");
    if(std::string(SelectionEfficiency->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEfficiency->GetXaxis()->SetBinLabel(nStep,"Pruned CA8 4V");
    nstepEvents[nStep-1]++;
    nStep = 6;

   
    // Kinematic Fit  for each neutrino type --> just type 0  and type 2                                                                                                                
    TLorentzVector fit_electron_type0(0,0,0,0), fit_neutrino_type0(0,0,0,0), fit_W_subjet1_type0(0,0,0,0), fit_W_subjet2_type0(0,0,0,0) ;
    TLorentzVector fit_electron_type2(0,0,0,0), fit_neutrino_type2(0,0,0,0), fit_W_subjet1_type2(0,0,0,0), fit_W_subjet2_type2(0,0,0,0) ;
    TLorentzVector fit_electron_type0_met(0,0,0,0), fit_neutrino_type0_met(0,0,0,0), fit_W_subjet1_type0_met(0,0,0,0), fit_W_subjet2_type0_met(0,0,0,0) ;
    TLorentzVector fit_electron_type2_met(0,0,0,0), fit_neutrino_type2_met(0,0,0,0), fit_W_subjet1_type2_met(0,0,0,0), fit_W_subjet2_type2_met(0,0,0,0) ;

    doKinematicFit(1, W_electron, W_neutrino_type0, W_subjet1, W_subjet2,  fit_electron_type0, fit_neutrino_type0, fit_W_subjet1_type0, fit_W_subjet2_type0, NewElectronTree->fit_chi2_type0, 
                   NewElectronTree->fit_NDF_type0, NewElectronTree->fit_status_type0, LeptonType);
    doKinematicFit(1, W_electron, W_neutrino_type2, W_subjet1, W_subjet2,  fit_electron_type2, fit_neutrino_type2, fit_W_subjet1_type2, fit_W_subjet2_type2, NewElectronTree->fit_chi2_type2, 
                   NewElectronTree->fit_NDF_type2, NewElectronTree->fit_status_type2, LeptonType);
    doKinematicFit(1, W_electron, W_neutrino_type0_met, W_subjet1, W_subjet2,  fit_electron_type0_met, fit_neutrino_type0_met, fit_W_subjet1_type0_met, fit_W_subjet2_type0_met, 
                   NewElectronTree->fit_chi2_type0_met, NewElectronTree->fit_NDF_type0_met, NewElectronTree->fit_status_type0_met, LeptonType);
    doKinematicFit(1, W_electron, W_neutrino_type2_met, W_subjet1, W_subjet2,  fit_electron_type2_met, fit_neutrino_type2_met, fit_W_subjet1_type2_met, fit_W_subjet2_type2_met,
                   NewElectronTree->fit_chi2_type2_met, NewElectronTree->fit_NDF_type2_met, NewElectronTree->fit_status_type2_met, LeptonType);
    
    if(fit_electron_type0.Pt() >0 && fit_neutrino_type0.Pt()>0 && fit_W_subjet1_type0.Pt()>0 && fit_W_subjet2_type0.Pt()>0){

     NewElectronTree->fit_el_px_type0 = fit_electron_type0.Px();
     NewElectronTree->fit_el_py_type0 = fit_electron_type0.Py(); 
     NewElectronTree->fit_el_pz_type0 = fit_electron_type0.Pz(); 
     NewElectronTree->fit_el_e_type0  = fit_electron_type0.E();
     NewElectronTree->fit_nv_px_type0 = fit_neutrino_type0.Px(); 
     NewElectronTree->fit_nv_py_type0 = fit_neutrino_type0.Py(); 
     NewElectronTree->fit_nv_pz_type0 = fit_neutrino_type0.Pz(); 
     NewElectronTree->fit_nv_e_type0  = fit_neutrino_type0.E();

     NewElectronTree->fit_subjet1_px_type0 = fit_W_subjet1_type0.Px();  NewElectronTree->fit_subjet2_px_type0 = fit_W_subjet2_type0.Px();
     NewElectronTree->fit_subjet1_py_type0 = fit_W_subjet1_type0.Py();  NewElectronTree->fit_subjet2_py_type0 = fit_W_subjet2_type0.Py();  
     NewElectronTree->fit_subjet1_pz_type0 = fit_W_subjet1_type0.Pz();  NewElectronTree->fit_subjet2_pz_type0 = fit_W_subjet2_type0.Pz();
     NewElectronTree->fit_subjet1_e_type0  = fit_W_subjet1_type0.E();   NewElectronTree->fit_subjet2_e_type0  = fit_W_subjet2_type0.E();
     NewElectronTree->fit_subjet1_m_type0  = fit_W_subjet1_type0.M();   NewElectronTree->fit_subjet2_m_type0  = fit_W_subjet2_type0.M();

     NewElectronTree->fit_lvj_m_type0   = (fit_electron_type0+fit_neutrino_type0+fit_W_subjet1_type0+fit_W_subjet2_type0).M();
     NewElectronTree->fit_lv_m_type0    = (fit_electron_type0+fit_neutrino_type0).M();
     NewElectronTree->fit_j_m_type0     = (fit_W_subjet1_type0+fit_W_subjet2_type0).M();
     NewElectronTree->fit_lvj_pt_type0  = (fit_electron_type0+fit_neutrino_type0+fit_W_subjet1_type0+fit_W_subjet2_type0).M();
     NewElectronTree->fit_lvj_eta_type0 = (fit_electron_type0+fit_neutrino_type0+fit_W_subjet1_type0+fit_W_subjet2_type0).Eta();
     NewElectronTree->fit_lvj_phi_type0 = (fit_electron_type0+fit_neutrino_type0+fit_W_subjet1_type0+fit_W_subjet2_type0).Phi();
     NewElectronTree->fit_lvj_e_type0   = (fit_electron_type0+fit_neutrino_type0+fit_W_subjet1_type0+fit_W_subjet2_type0).E();
    }

    if(fit_electron_type2.Pt() >0 && fit_neutrino_type2.Pt()>0 && fit_W_subjet1_type2.Pt()>0 && fit_W_subjet2_type2.Pt()>0){

     NewElectronTree->fit_el_px_type2 = fit_electron_type2.Px();
     NewElectronTree->fit_el_py_type2 = fit_electron_type2.Py(); 
     NewElectronTree->fit_el_pz_type2 = fit_electron_type2.Pz(); 
     NewElectronTree->fit_el_e_type2  = fit_electron_type2.E();
     NewElectronTree->fit_nv_px_type2 = fit_neutrino_type2.Px(); 
     NewElectronTree->fit_nv_py_type2 = fit_neutrino_type2.Py(); 
     NewElectronTree->fit_nv_pz_type2 = fit_neutrino_type2.Pz(); 
     NewElectronTree->fit_nv_e_type2  = fit_neutrino_type2.E();

     NewElectronTree->fit_subjet1_px_type2 = fit_W_subjet1_type2.Px();  NewElectronTree->fit_subjet2_px_type2 = fit_W_subjet2_type2.Px();
     NewElectronTree->fit_subjet1_py_type2 = fit_W_subjet1_type2.Py();  NewElectronTree->fit_subjet2_py_type2 = fit_W_subjet2_type2.Py();  
     NewElectronTree->fit_subjet1_pz_type2 = fit_W_subjet1_type2.Pz();  NewElectronTree->fit_subjet2_pz_type2 = fit_W_subjet2_type2.Pz();
     NewElectronTree->fit_subjet1_e_type2  = fit_W_subjet1_type2.E();   NewElectronTree->fit_subjet2_e_type2  = fit_W_subjet2_type2.E();
     NewElectronTree->fit_subjet1_m_type2  = fit_W_subjet1_type2.M();   NewElectronTree->fit_subjet2_m_type2  = fit_W_subjet2_type2.M();

     NewElectronTree->fit_lvj_m_type2   = (fit_electron_type2+fit_neutrino_type2+fit_W_subjet1_type2+fit_W_subjet2_type2).M();
     NewElectronTree->fit_lv_m_type2    = (fit_electron_type2+fit_neutrino_type2).M();
     NewElectronTree->fit_j_m_type2     = (fit_W_subjet1_type2+fit_W_subjet2_type2).M();
     NewElectronTree->fit_lvj_pt_type2  = (fit_electron_type2+fit_neutrino_type2+fit_W_subjet1_type2+fit_W_subjet2_type2).M();
     NewElectronTree->fit_lvj_eta_type2 = (fit_electron_type2+fit_neutrino_type2+fit_W_subjet1_type2+fit_W_subjet2_type2).Eta();
     NewElectronTree->fit_lvj_phi_type2 = (fit_electron_type2+fit_neutrino_type2+fit_W_subjet1_type2+fit_W_subjet2_type2).Phi();
     NewElectronTree->fit_lvj_e_type2   = (fit_electron_type2+fit_neutrino_type2+fit_W_subjet1_type2+fit_W_subjet2_type2).E();
    }

    if(fit_electron_type0_met.Pt() >0 && fit_neutrino_type0_met.Pt()>0 && fit_W_subjet1_type0_met.Pt()>0 && fit_W_subjet2_type0_met.Pt()>0){

     NewElectronTree->fit_el_px_type0_met = fit_electron_type0_met.Px();
     NewElectronTree->fit_el_py_type0_met = fit_electron_type0_met.Py(); 
     NewElectronTree->fit_el_pz_type0_met = fit_electron_type0_met.Pz(); 
     NewElectronTree->fit_el_e_type0_met  = fit_electron_type0_met.E();
     NewElectronTree->fit_nv_px_type0_met = fit_neutrino_type0_met.Px(); 
     NewElectronTree->fit_nv_py_type0_met = fit_neutrino_type0_met.Py(); 
     NewElectronTree->fit_nv_pz_type0_met = fit_neutrino_type0_met.Pz(); 
     NewElectronTree->fit_nv_e_type0_met  = fit_neutrino_type0_met.E();

     NewElectronTree->fit_subjet1_px_type0_met = fit_W_subjet1_type0_met.Px();  NewElectronTree->fit_subjet2_px_type0_met = fit_W_subjet2_type0_met.Px();
     NewElectronTree->fit_subjet1_py_type0_met = fit_W_subjet1_type0_met.Py();  NewElectronTree->fit_subjet2_py_type0_met = fit_W_subjet2_type0_met.Py();  
     NewElectronTree->fit_subjet1_pz_type0_met = fit_W_subjet1_type0_met.Pz();  NewElectronTree->fit_subjet2_pz_type0_met = fit_W_subjet2_type0_met.Pz();
     NewElectronTree->fit_subjet1_e_type0_met  = fit_W_subjet1_type0_met.E();   NewElectronTree->fit_subjet2_e_type0_met  = fit_W_subjet2_type0_met.E();
     NewElectronTree->fit_subjet1_m_type0_met  = fit_W_subjet1_type0_met.M();   NewElectronTree->fit_subjet2_m_type0_met  = fit_W_subjet2_type0_met.M();

     NewElectronTree->fit_lvj_m_type0_met   = (fit_electron_type0_met+fit_neutrino_type0_met+fit_W_subjet1_type0_met+fit_W_subjet2_type0_met).M();
     NewElectronTree->fit_lv_m_type0_met    = (fit_electron_type0_met+fit_neutrino_type0_met).M();
     NewElectronTree->fit_j_m_type0_met     = (fit_W_subjet1_type0_met+fit_W_subjet2_type0_met).M();
     NewElectronTree->fit_lvj_pt_type0_met  = (fit_electron_type0_met+fit_neutrino_type0_met+fit_W_subjet1_type0_met+fit_W_subjet2_type0_met).M();
     NewElectronTree->fit_lvj_eta_type0_met = (fit_electron_type0_met+fit_neutrino_type0_met+fit_W_subjet1_type0_met+fit_W_subjet2_type0_met).Eta();
     NewElectronTree->fit_lvj_phi_type0_met = (fit_electron_type0_met+fit_neutrino_type0_met+fit_W_subjet1_type0_met+fit_W_subjet2_type0_met).Phi();
     NewElectronTree->fit_lvj_e_type0_met   = (fit_electron_type0_met+fit_neutrino_type0_met+fit_W_subjet1_type0_met+fit_W_subjet2_type0_met).E();
    }


    if(fit_electron_type2_met.Pt() >0 && fit_neutrino_type2_met.Pt()>0 && fit_W_subjet1_type2_met.Pt()>0 && fit_W_subjet2_type2_met.Pt()>0){

     NewElectronTree->fit_el_px_type2_met = fit_electron_type2_met.Px();
     NewElectronTree->fit_el_py_type2_met = fit_electron_type2_met.Py(); 
     NewElectronTree->fit_el_pz_type2_met = fit_electron_type2_met.Pz(); 
     NewElectronTree->fit_el_e_type2_met  = fit_electron_type2_met.E();
     NewElectronTree->fit_nv_px_type2_met = fit_neutrino_type2_met.Px(); 
     NewElectronTree->fit_nv_py_type2_met = fit_neutrino_type2_met.Py(); 
     NewElectronTree->fit_nv_pz_type2_met = fit_neutrino_type2_met.Pz(); 
     NewElectronTree->fit_nv_e_type2_met  = fit_neutrino_type2_met.E();

     NewElectronTree->fit_subjet1_px_type2_met = fit_W_subjet1_type2_met.Px();  NewElectronTree->fit_subjet2_px_type2_met = fit_W_subjet2_type2_met.Px();
     NewElectronTree->fit_subjet1_py_type2_met = fit_W_subjet1_type2_met.Py();  NewElectronTree->fit_subjet2_py_type2_met = fit_W_subjet2_type2_met.Py();  
     NewElectronTree->fit_subjet1_pz_type2_met = fit_W_subjet1_type2_met.Pz();  NewElectronTree->fit_subjet2_pz_type2_met = fit_W_subjet2_type2_met.Pz();
     NewElectronTree->fit_subjet1_e_type2_met  = fit_W_subjet1_type2_met.E();   NewElectronTree->fit_subjet2_e_type2_met  = fit_W_subjet2_type2_met.E();
     NewElectronTree->fit_subjet1_m_type2_met  = fit_W_subjet1_type2_met.M();   NewElectronTree->fit_subjet2_m_type2_met  = fit_W_subjet2_type2_met.M();

     NewElectronTree->fit_lvj_m_type2_met   = (fit_electron_type2_met+fit_neutrino_type2_met+fit_W_subjet1_type2_met+fit_W_subjet2_type2_met).M();
     NewElectronTree->fit_lv_m_type2_met    = (fit_electron_type2_met+fit_neutrino_type2_met).M();
     NewElectronTree->fit_j_m_type2_met     = (fit_W_subjet1_type2_met+fit_W_subjet2_type2_met).M();
     NewElectronTree->fit_lvj_pt_type2_met  = (fit_electron_type2_met+fit_neutrino_type2_met+fit_W_subjet1_type2_met+fit_W_subjet2_type2_met).M();
     NewElectronTree->fit_lvj_eta_type2_met = (fit_electron_type2_met+fit_neutrino_type2_met+fit_W_subjet1_type2_met+fit_W_subjet2_type2_met).Eta();
     NewElectronTree->fit_lvj_phi_type2_met = (fit_electron_type2_met+fit_neutrino_type2_met+fit_W_subjet1_type2_met+fit_W_subjet2_type2_met).Phi();
     NewElectronTree->fit_lvj_e_type2_met   = (fit_electron_type2_met+fit_neutrino_type2_met+fit_W_subjet1_type2_met+fit_W_subjet2_type2_met).E();
    }

    
    NewElectronTree->boosted_lvj_m_type0   = (W_electron+W_neutrino_type0+W_subjet1+W_subjet2).M();
    NewElectronTree->boosted_lv_m_type0    = (W_electron+W_neutrino_type0).M();
    NewElectronTree->boosted_j_m_type0     = (W_subjet1+W_subjet2).M();
    NewElectronTree->boosted_lvj_pt_type0  = (W_electron+W_neutrino_type0+W_subjet1+W_subjet2).Pt();
    NewElectronTree->boosted_lvj_eta_type0 = (W_electron+W_neutrino_type0+W_subjet1+W_subjet2).Eta();
    NewElectronTree->boosted_lvj_phi_type0 = (W_electron+W_neutrino_type0+W_subjet1+W_subjet2).Phi();
    NewElectronTree->boosted_lvj_e_type0   = (W_electron+W_neutrino_type0+W_subjet1+W_subjet2).E();
 
    NewElectronTree->boostedW_lvj_m_type0   = (W_electron+W_neutrino_type0+W_GroomedJet_CA8).M();
    NewElectronTree->boostedW_lv_m_type0    = (W_electron+W_neutrino_type0).M();
    NewElectronTree->boostedW_j_m_type0     = ElectronTree->fReader->getFloat("GroomedJet_CA8_mass_pr")[0];
    NewElectronTree->boostedW_lvj_pt_type0  = (W_electron+W_neutrino_type0+W_GroomedJet_CA8).Pt();
    NewElectronTree->boostedW_lvj_eta_type0 = (W_electron+W_neutrino_type0+W_GroomedJet_CA8).Eta();
    NewElectronTree->boostedW_lvj_phi_type0 = (W_electron+W_neutrino_type0+W_GroomedJet_CA8).Phi();
    NewElectronTree->boostedW_lvj_e_type0   = (W_electron+W_neutrino_type0+W_GroomedJet_CA8).E();

    NewElectronTree->boosted_lvj_m_type2   = (W_electron+W_neutrino_type2+W_subjet1+W_subjet2).M();
    NewElectronTree->boosted_lv_m_type2    = (W_electron+W_neutrino_type2).M();
    NewElectronTree->boosted_j_m_type2     = (W_subjet1+W_subjet2).M();
    NewElectronTree->boosted_lvj_pt_type2  = (W_electron+W_neutrino_type2+W_subjet1+W_subjet2).Pt();
    NewElectronTree->boosted_lvj_eta_type2 = (W_electron+W_neutrino_type2+W_subjet1+W_subjet2).Eta();
    NewElectronTree->boosted_lvj_phi_type2 = (W_electron+W_neutrino_type2+W_subjet1+W_subjet2).Phi();
    NewElectronTree->boosted_lvj_e_type2   = (W_electron+W_neutrino_type2+W_subjet1+W_subjet2).E();
 
    NewElectronTree->boostedW_lvj_m_type2   = (W_electron+W_neutrino_type2+W_GroomedJet_CA8).M();
    NewElectronTree->boostedW_lv_m_type2    = (W_electron+W_neutrino_type2).M();
    NewElectronTree->boostedW_j_m_type2     = ElectronTree->fReader->getFloat("GroomedJet_CA8_mass_pr")[0];
    NewElectronTree->boostedW_lvj_pt_type2  = (W_electron+W_neutrino_type2+W_GroomedJet_CA8).Pt();
    NewElectronTree->boostedW_lvj_eta_type2 = (W_electron+W_neutrino_type2+W_GroomedJet_CA8).Eta();
    NewElectronTree->boostedW_lvj_phi_type2 = (W_electron+W_neutrino_type2+W_GroomedJet_CA8).Phi();
    NewElectronTree->boostedW_lvj_e_type2   = (W_electron+W_neutrino_type2+W_GroomedJet_CA8).E();

    NewElectronTree->boosted_lvj_m_type0_met   = (W_electron+W_neutrino_type0_met+W_subjet1+W_subjet2).M();
    NewElectronTree->boosted_lv_m_type0_met    = (W_electron+W_neutrino_type0_met).M();
    NewElectronTree->boosted_j_m_type0_met     = (W_subjet1+W_subjet2).M();
    NewElectronTree->boosted_lvj_pt_type0_met  = (W_electron+W_neutrino_type0_met+W_subjet1+W_subjet2).Pt();
    NewElectronTree->boosted_lvj_eta_type0_met = (W_electron+W_neutrino_type0_met+W_subjet1+W_subjet2).Eta();
    NewElectronTree->boosted_lvj_phi_type0_met = (W_electron+W_neutrino_type0_met+W_subjet1+W_subjet2).Phi();
    NewElectronTree->boosted_lvj_e_type0_met   = (W_electron+W_neutrino_type0_met+W_subjet1+W_subjet2).E();
 
    NewElectronTree->boostedW_lvj_m_type0_met   = (W_electron+W_neutrino_type0_met+W_GroomedJet_CA8).M();
    NewElectronTree->boostedW_lv_m_type0_met    = (W_electron+W_neutrino_type0_met).M();
    NewElectronTree->boostedW_j_m_type0_met     = ElectronTree->fReader->getFloat("GroomedJet_CA8_mass_pr")[0];
    NewElectronTree->boostedW_lvj_pt_type0_met  = (W_electron+W_neutrino_type0_met+W_GroomedJet_CA8).Pt();
    NewElectronTree->boostedW_lvj_eta_type0_met = (W_electron+W_neutrino_type0_met+W_GroomedJet_CA8).Eta();
    NewElectronTree->boostedW_lvj_phi_type0_met = (W_electron+W_neutrino_type0_met+W_GroomedJet_CA8).Phi();
    NewElectronTree->boostedW_lvj_e_type0_met   = (W_electron+W_neutrino_type0_met+W_GroomedJet_CA8).E();

    NewElectronTree->boosted_lvj_m_type2_met   = (W_electron+W_neutrino_type2_met+W_subjet1+W_subjet2).M();
    NewElectronTree->boosted_lv_m_type2_met    = (W_electron+W_neutrino_type2_met).M();
    NewElectronTree->boosted_j_m_type2_met     = (W_subjet1+W_subjet2).M();
    NewElectronTree->boosted_lvj_pt_type2_met  = (W_electron+W_neutrino_type2_met+W_subjet1+W_subjet2).Pt();
    NewElectronTree->boosted_lvj_eta_type2_met = (W_electron+W_neutrino_type2_met+W_subjet1+W_subjet2).Eta();
    NewElectronTree->boosted_lvj_phi_type2_met = (W_electron+W_neutrino_type2_met+W_subjet1+W_subjet2).Phi();
    NewElectronTree->boosted_lvj_e_type2_met   = (W_electron+W_neutrino_type2_met+W_subjet1+W_subjet2).E();
 
    NewElectronTree->boostedW_lvj_m_type2_met   = (W_electron+W_neutrino_type2_met+W_GroomedJet_CA8).M();
    NewElectronTree->boostedW_lv_m_type2_met    = (W_electron+W_neutrino_type2_met).M();
    NewElectronTree->boostedW_j_m_type2_met     = ElectronTree->fReader->getFloat("GroomedJet_CA8_mass_pr")[0];
    NewElectronTree->boostedW_lvj_pt_type2_met  = (W_electron+W_neutrino_type2_met+W_GroomedJet_CA8).Pt();
    NewElectronTree->boostedW_lvj_eta_type2_met = (W_electron+W_neutrino_type2_met+W_GroomedJet_CA8).Eta();
    NewElectronTree->boostedW_lvj_phi_type2_met = (W_electron+W_neutrino_type2_met+W_GroomedJet_CA8).Phi();
    
    
    // Angles for the central Higgs Kinematics
    double costheta1, costheta2, phi, costhetastar, phistar1, phistar2;

    //Use the Subjet in the Boosted W Analyisis                                                                                                                                             
    if (ElectronTree->fReader->getFloat("W_electron_charge")[0] < 0) calculateAngles(W_electron, W_neutrino_type0,W_subjet1,W_subjet2,costheta1,costheta2,phi,costhetastar,phistar1,phistar2);
    else calculateAngles(W_neutrino_type0, W_electron, W_subjet1, W_subjet2, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);
   
    NewElectronTree->boosted_wjj_ang_ha_type0   = costheta1;
    NewElectronTree->boosted_wjj_ang_hb_type0   = fabs(costheta2); 
    NewElectronTree->boosted_wjj_ang_hs_type0   = costhetastar;
    NewElectronTree->boosted_wjj_ang_phi_type0  = phi;
    NewElectronTree->boosted_wjj_ang_phia_type0 = phistar1;																
    NewElectronTree->boosted_wjj_ang_phib_type0 = phistar2;

    if (ElectronTree->fReader->getFloat("W_electron_charge")[0] < 0) calculateAngles(W_electron, W_neutrino_type2,W_subjet1,W_subjet2,costheta1,costheta2,phi,costhetastar,phistar1,phistar2);
    else calculateAngles(W_neutrino_type2, W_electron, W_subjet1, W_subjet2, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);
   
    NewElectronTree->boosted_wjj_ang_ha_type2   = costheta1;
    NewElectronTree->boosted_wjj_ang_hb_type2   = fabs(costheta2); 
    NewElectronTree->boosted_wjj_ang_hs_type2   = costhetastar;
    NewElectronTree->boosted_wjj_ang_phi_type2  = phi;
    NewElectronTree->boosted_wjj_ang_phia_type2 = phistar1;															     
    NewElectronTree->boosted_wjj_ang_phib_type2 = phistar2;

    if (ElectronTree->fReader->getFloat("W_electron_charge")[0] < 0) calculateAngles(W_electron, W_neutrino_type0_met,W_subjet1,W_subjet2,costheta1,costheta2,phi,costhetastar,phistar1,phistar2);
    else calculateAngles(W_neutrino_type0_met, W_electron, W_subjet1, W_subjet2, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);
   
    NewElectronTree->boosted_wjj_ang_ha_type0_met   = costheta1;
    NewElectronTree->boosted_wjj_ang_hb_type0_met   = fabs(costheta2); 
    NewElectronTree->boosted_wjj_ang_hs_type0_met   = costhetastar;
    NewElectronTree->boosted_wjj_ang_phi_type0_met  = phi;
    NewElectronTree->boosted_wjj_ang_phia_type0_met = phistar1;															
    NewElectronTree->boosted_wjj_ang_phib_type0_met = phistar2;

    if (ElectronTree->fReader->getFloat("W_electron_charge")[0] < 0) calculateAngles(W_electron, W_neutrino_type2_met,W_subjet1,W_subjet2,costheta1,costheta2,phi,costhetastar,phistar1,phistar2);
    else calculateAngles(W_neutrino_type2_met, W_electron, W_subjet1, W_subjet2, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);
   
    NewElectronTree->boosted_wjj_ang_ha_type2_met   = costheta1;
    NewElectronTree->boosted_wjj_ang_hb_type2_met   = fabs(costheta2); 
    NewElectronTree->boosted_wjj_ang_hs_type2_met   = costhetastar;
    NewElectronTree->boosted_wjj_ang_phi_type2_met  = phi;
    NewElectronTree->boosted_wjj_ang_phia_type2_met = phistar1;														
    NewElectronTree->boosted_wjj_ang_phib_type2_met = phistar2;

    // Clean AK5 Jet Collection from the hadronic W and sotre the jet binning
    std::vector<int> numberJetBin ;
    for( size_t iJetPtCutMin = 0; iJetPtCutMin < JetPtCutMin.size(); iJetPtCutMin++){
      for(size_t iJet = 0; iJet < JetPFCor_AK5_Collection.at(iJetPtCutMin).size() ; iJet ++){       
	if(deltaR(JetPFCor_AK5_Collection.at(iJetPtCutMin).at(iJet).Momentum_.Phi(),GroomedJet_CA8_Collection.at(iJetPtCutMin).at(0).Phi(),
		  JetPFCor_AK5_Collection.at(iJetPtCutMin).at(iJet).Momentum_.Eta(),GroomedJet_CA8_Collection.at(iJetPtCutMin).at(0).Eta()) < CleaningTreshold ){
	  HadronicW_AK5_Collection.at(iJetPtCutMin).push_back(JetPFCor_AK5_Collection.at(iJetPtCutMin).at(iJet)); continue ;}

	CleanedJetPFCor_AK5_Collection.at(iJetPtCutMin).push_back(JetPFCor_AK5_Collection.at(iJetPtCutMin).at(iJet));

      }
 
      if(!CleanedJetPFCor_AK5_Collection.at(iJetPtCutMin).empty())
       numberJetBin.push_back(CleanedJetPFCor_AK5_Collection.at(iJetPtCutMin).size());
      else
       numberJetBin.push_back(0);

    }

    std::vector<int> numberJetBinGen ;
    if(NewElectronTree->fTree->FindBranch("GenGroomedJet_CA8_pt") && NewElectronTree->fTree->FindBranch("JetGen_Pt")){
     for( size_t iJetPtCutMin = 0; iJetPtCutMin < JetPtCutMin.size(); iJetPtCutMin++){
      for(size_t iJet = 0; iJet < GenJetPFCor_AK5_Collection.at(iJetPtCutMin).size() ; iJet ++){
	if(GenGroomedJet_CA8_Collection.at(iJetPtCutMin).empty()) break;
	if(deltaR(GenJetPFCor_AK5_Collection.at(iJetPtCutMin).at(iJet).Momentum_.Phi(),GenGroomedJet_CA8_Collection.at(iJetPtCutMin).at(0).Phi(),
		  GenJetPFCor_AK5_Collection.at(iJetPtCutMin).at(iJet).Momentum_.Eta(),GenGroomedJet_CA8_Collection.at(iJetPtCutMin).at(0).Eta()) < CleaningTreshold ){
	  GenHadronicW_AK5_Collection.at(iJetPtCutMin).push_back(GenJetPFCor_AK5_Collection.at(iJetPtCutMin).at(iJet)); continue ;}

	GenCleanedJetPFCor_AK5_Collection.at(iJetPtCutMin).push_back(GenJetPFCor_AK5_Collection.at(iJetPtCutMin).at(iJet));
        
      }
      
      if(!GenCleanedJetPFCor_AK5_Collection.at(iJetPtCutMin).empty())
	numberJetBinGen.push_back(GenCleanedJetPFCor_AK5_Collection.at(iJetPtCutMin).size());
      else
       numberJetBinGen.push_back(0);
     
     }
    }

    NewElectronTree -> numberJetBin = numberJetBin;
    NewElectronTree -> numberJetBinGen = numberJetBinGen;

    if (NewElectronTree -> numberJetBin.at(0) == 1 && CleanedJetPFCor_AK5_Collection.at(0).size() == 1){

      std::sort(CleanedJetPFCor_AK5_Collection.at(0).begin(),CleanedJetPFCor_AK5_Collection.at(0).end(),TLVP_PtSort());
      std::sort(GenCleanedJetPFCor_AK5_Collection.at(0).begin(),GenCleanedJetPFCor_AK5_Collection.at(0).end(),TLVP_PtSort());

      NewElectronTree->vbf_maxpt_j1_e   = CleanedJetPFCor_AK5_Collection.at(0).at(0).Momentum_.E();
      NewElectronTree->vbf_maxpt_j1_pt  = CleanedJetPFCor_AK5_Collection.at(0).at(0).Momentum_.Pt();
      NewElectronTree->vbf_maxpt_j1_eta = CleanedJetPFCor_AK5_Collection.at(0).at(0).Momentum_.Eta();
      NewElectronTree->vbf_maxpt_j1_phi = CleanedJetPFCor_AK5_Collection.at(0).at(0).Momentum_.Phi();
      NewElectronTree->vbf_maxpt_j1_m   = CleanedJetPFCor_AK5_Collection.at(0).at(0).Momentum_.M(); 
      
      if(!numberJetBinGen.empty()){
	if(numberJetBinGen.at(0)==1 && GenCleanedJetPFCor_AK5_Collection.at(0).size()==1){
        NewElectronTree->vbf_maxpt_j1_e_gen   = GenCleanedJetPFCor_AK5_Collection.at(0).at(0).Momentum_.E();
        NewElectronTree->vbf_maxpt_j1_pt_gen  = GenCleanedJetPFCor_AK5_Collection.at(0).at(0).Momentum_.Pt();
        NewElectronTree->vbf_maxpt_j1_eta_gen = GenCleanedJetPFCor_AK5_Collection.at(0).at(0).Momentum_.Eta();
        NewElectronTree->vbf_maxpt_j1_phi_gen = GenCleanedJetPFCor_AK5_Collection.at(0).at(0).Momentum_.Phi();
        NewElectronTree->vbf_maxpt_j1_m_gen   = GenCleanedJetPFCor_AK5_Collection.at(0).at(0).Momentum_.M(); 
        NewElectronTree->vbf_maxpt_j1_bDiscriminatorSSVHE_gen = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorSSVHE")[GenCleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
        NewElectronTree->vbf_maxpt_j1_bDiscriminatorTCHE_gen  = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorTCHE")[GenCleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
        NewElectronTree->vbf_maxpt_j1_bDiscriminatorCSV_gen   = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorCSV")[GenCleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
        NewElectronTree->vbf_maxpt_j1_bDiscriminatorSSVHP_gen = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorSSVHP")[GenCleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
        NewElectronTree->vbf_maxpt_j1_bDiscriminatorTCHP_gen  = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorTCHP")[GenCleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
       }
      }
      
      NewElectronTree->vbf_maxpt_j1_QGLikelihood = ElectronTree->fReader->getFloat("JetPFCor_QGLikelihood")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;

      NewElectronTree->vbf_maxpt_j1_isPileUpLoose  = ElectronTree->fReader->getBool("JetPFCor_isPileUpJetLoose")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_isPileUpMedium = ElectronTree->fReader->getBool("JetPFCor_isPileUpJetMedium")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_isPileUpTight  = ElectronTree->fReader->getBool("JetPFCor_isPileUpJetTight")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;

      NewElectronTree->vbf_maxpt_j1_bDiscriminatorSSVHE = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorSSVHE")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_bDiscriminatorTCHE  = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorTCHE")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_bDiscriminatorCSV   = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorCSV")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_bDiscriminatorSSVHP = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorSSVHP")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_bDiscriminatorTCHP  = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorTCHP")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;

      NewElectronTree->vbf_maxpt_j1_ChargedHadronEnergy      = ElectronTree->fReader->getFloat("JetPFCor_ChargedHadronEnergy")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_ChargedHadronEnergyFrac  = ElectronTree->fReader->getFloat("JetPFCor_ChargedHadronEnergyFrac")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_NeutralHadronEnergy      = ElectronTree->fReader->getFloat("JetPFCor_NeutralHadronEnergy")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_NeutralHadronEnergyFrac  = ElectronTree->fReader->getFloat("JetPFCor_NeutralHadronEnergyFrac")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_ChargedEmEnergy          = ElectronTree->fReader->getFloat("JetPFCor_ChargedEmEnergy")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_ChargedEmEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCor_ChargedEmEnergyFrac")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_ChargedMuEnergy          = ElectronTree->fReader->getFloat("JetPFCor_ChargedMuEnergy")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_ChargedMuEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCor_ChargedMuEnergyFrac")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_NeutralEmEnergy          = ElectronTree->fReader->getFloat("JetPFCor_NeutralEmEnergy")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_NeutralEmEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCor_NeutralEmEnergyFrac")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_PhotonEnergy             = ElectronTree->fReader->getFloat("JetPFCor_PhotonEnergy")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_PhotonEnergyFraction     = ElectronTree->fReader->getFloat("JetPFCor_PhotonEnergyFraction")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_ElectronEnergy           = ElectronTree->fReader->getFloat("JetPFCor_ElectronEnergy")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_ElectronEnergyFraction   = ElectronTree->fReader->getFloat("JetPFCor_ElectronEnergyFraction")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_HFHadronEnergy           = ElectronTree->fReader->getFloat("JetPFCor_HFHadronEnergy")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_HFHadronEnergyFraction   = ElectronTree->fReader->getFloat("JetPFCor_HFHadronEnergyFraction")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_HFEMEnergy               = ElectronTree->fReader->getFloat("JetPFCor_HFEMEnergy")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_HFEMEnergyFraction       = ElectronTree->fReader->getFloat("JetPFCor_HFEMEnergyFraction")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;

      NewElectronTree->vbf_maxpt_j1_ChargedMultiplicity        = ElectronTree->fReader->getFloat("JetPFCor_ChargedMultiplicity")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_NeutralMultiplicity        = ElectronTree->fReader->getFloat("JetPFCor_NeutralMultiplicity")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_MuonMultiplicity           = ElectronTree->fReader->getFloat("JetPFCor_MuonMultiplicity")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_ChargedHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCor_ChargedHadronMultiplicity")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_NeutralHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCor_NeutralHadronMultiplicity")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_PhotonMultiplicity         = ElectronTree->fReader->getFloat("JetPFCor_PhotonMultiplicity")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ; 
      NewElectronTree->vbf_maxpt_j1_MuonMultiplicity           = ElectronTree->fReader->getFloat("JetPFCor_MuonMultiplicity")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_HFHadronMultiplicity       = ElectronTree->fReader->getFloat("JetPFCor_HFHadronMultiplicity")[CleanedJetPFCor_AK5_Collection.at(0).at(0).position_] ;


    }
    /// store info only for the VBF case 
    else if (NewElectronTree -> numberJetBin.at(0) >= 2){
    
     // vbf Tag Jet Selection
    
     std::vector<JetAK5> outputAK5_PtSorted;
     std::vector<JetAK5> outputAK5_DEtaSorted;
     std::vector<JetAK5> outputAK5_MjjSorted;

     std::vector<JetAK5> outputGenAK5_PtSorted;
     std::vector<JetAK5> outputGenAK5_DEtaSorted;
     std::vector<JetAK5> outputGenAK5_MjjSorted;

     // Sorting of AK5 Cleaned Collection in Pt

     std::sort(CleanedJetPFCor_AK5_Collection.at(0).begin(),CleanedJetPFCor_AK5_Collection.at(0).end(),TLVP_PtSort());
     outputAK5_PtSorted.push_back(CleanedJetPFCor_AK5_Collection.at(0).at(0));
     outputAK5_PtSorted.push_back(CleanedJetPFCor_AK5_Collection.at(0).at(1));
     if(outputAK5_PtSorted.size() < 2) continue ;
     
     if(!numberJetBinGen.empty() && !GenCleanedJetPFCor_AK5_Collection.at(0).empty()){
      if(numberJetBinGen.at(0)>=2 && GenCleanedJetPFCor_AK5_Collection.at(0).size()>=2){
       std::sort(GenCleanedJetPFCor_AK5_Collection.at(0).begin(),GenCleanedJetPFCor_AK5_Collection.at(0).end(),TLVP_PtSort());
       outputGenAK5_PtSorted.push_back(GenCleanedJetPFCor_AK5_Collection.at(0).at(0));
       outputGenAK5_PtSorted.push_back(GenCleanedJetPFCor_AK5_Collection.at(0).at(1));
      }
     }
         
     // Sorting of AK5 Cleaned Collection in DeltaEta

     std::sort(CleanedJetPFCor_AK5_Collection.at(0).begin(),CleanedJetPFCor_AK5_Collection.at(0).end(),TLVP_EtaSort());
     outputAK5_DEtaSorted.push_back(CleanedJetPFCor_AK5_Collection.at(0).front());
     outputAK5_DEtaSorted.push_back(CleanedJetPFCor_AK5_Collection.at(0).back());
     if(outputAK5_DEtaSorted.size() < 2) continue ;
     
     if(!numberJetBinGen.empty() && !GenCleanedJetPFCor_AK5_Collection.at(0).empty()){
      if(numberJetBinGen.at(0)>=2 && GenCleanedJetPFCor_AK5_Collection.at(0).size()>=2){
       std::sort(GenCleanedJetPFCor_AK5_Collection.at(0).begin(),GenCleanedJetPFCor_AK5_Collection.at(0).end(),TLVP_EtaSort());
       outputGenAK5_DEtaSorted.push_back(GenCleanedJetPFCor_AK5_Collection.at(0).front());
       outputGenAK5_DEtaSorted.push_back(GenCleanedJetPFCor_AK5_Collection.at(0).back());
      }
     }          
     // Sorting of AK5 Cleaned Collection in Mjj
     float maxMjj = 0. ;
     int iJ1 = 0 ;
     int iJ2 = 0 ;

     for (size_t iJet = 0 ; iJet < CleanedJetPFCor_AK5_Collection.at(0).size()-1 ; ++iJet){
       for (size_t jJet = iJet + 1 ; jJet < CleanedJetPFCor_AK5_Collection.at(0).size() ; ++jJet){

        TLorentzVector SumMomentum = CleanedJetPFCor_AK5_Collection.at(0).at(iJet).Momentum_ + CleanedJetPFCor_AK5_Collection.at(0).at(jJet).Momentum_ ;
        float Mjj = SumMomentum.M();
        if(Mjj > maxMjj){
         
          iJ1 = iJet ;
          iJ2 = jJet ;
          maxMjj = Mjj ;
        }
      }
     }

     outputAK5_MjjSorted.push_back (CleanedJetPFCor_AK5_Collection.at(0).at (iJ1)) ;
     outputAK5_MjjSorted.push_back (CleanedJetPFCor_AK5_Collection.at(0).at (iJ2)) ;
     if(outputAK5_MjjSorted.size() < 2) continue ;
    
     maxMjj = 0. ;
     iJ1 = 0 ; iJ2 = 0;

     if(!numberJetBinGen.empty() && !GenCleanedJetPFCor_AK5_Collection.at(0).empty()){
       if(numberJetBinGen.at(0)>=2 && GenCleanedJetPFCor_AK5_Collection.at(0).size()>=2){
       for (size_t iJet = 0 ; iJet < GenCleanedJetPFCor_AK5_Collection.at(0).size()-1 ; ++iJet){
        for (size_t jJet = iJet + 1 ; jJet < GenCleanedJetPFCor_AK5_Collection.at(0).size() ; ++jJet){

         TLorentzVector SumMomentum = GenCleanedJetPFCor_AK5_Collection.at(0).at(iJet).Momentum_ + GenCleanedJetPFCor_AK5_Collection.at(0).at(jJet).Momentum_ ;
         float Mjj = SumMomentum.M();
         if(Mjj > maxMjj){         
          iJ1 = iJet ;
          iJ2 = jJet ;
          maxMjj = Mjj ;
        }
       }
      }
      outputGenAK5_MjjSorted.push_back (GenCleanedJetPFCor_AK5_Collection.at(0).at (iJ1)) ;
      outputGenAK5_MjjSorted.push_back (GenCleanedJetPFCor_AK5_Collection.at(0).at (iJ2)) ;
      }     
     }
      
     //////////////////////////////////////////////////////////////////////////////////////////////////
     // Fill Information for Max Pt Pair of vbf tag jets
     //////////////////////////////////////////////////////////////////////////////////////////////////

     TLorentzVector Total4VMaxPt = outputAK5_PtSorted.at(0).Momentum_ + outputAK5_PtSorted.at(1).Momentum_ ;
    
     NewElectronTree->vbf_maxpt_jj_e   = Total4VMaxPt.E(); 
     NewElectronTree->vbf_maxpt_jj_pt  = Total4VMaxPt.Pt(); 
     NewElectronTree->vbf_maxpt_jj_eta = Total4VMaxPt.Eta(); 
     NewElectronTree->vbf_maxpt_jj_phi = Total4VMaxPt.Phi(); 
     NewElectronTree->vbf_maxpt_jj_m   = Total4VMaxPt.M(); 

     NewElectronTree->vbf_maxpt_j1_e   = outputAK5_PtSorted.at(0).Momentum_.E();
     NewElectronTree->vbf_maxpt_j1_pt  = outputAK5_PtSorted.at(0).Momentum_.Pt();
     NewElectronTree->vbf_maxpt_j1_eta = outputAK5_PtSorted.at(0).Momentum_.Eta();
     NewElectronTree->vbf_maxpt_j1_phi = outputAK5_PtSorted.at(0).Momentum_.Phi();
     NewElectronTree->vbf_maxpt_j1_m   = outputAK5_PtSorted.at(0).Momentum_.M(); 
 
     NewElectronTree->vbf_maxpt_j2_e   = outputAK5_PtSorted.at(1).Momentum_.E();
     NewElectronTree->vbf_maxpt_j2_pt  = outputAK5_PtSorted.at(1).Momentum_.Pt();
     NewElectronTree->vbf_maxpt_j2_eta = outputAK5_PtSorted.at(1).Momentum_.Eta();
     NewElectronTree->vbf_maxpt_j2_phi = outputAK5_PtSorted.at(1).Momentum_.Phi();
     NewElectronTree->vbf_maxpt_j2_m   = outputAK5_PtSorted.at(1).Momentum_.M();

   
     NewElectronTree->vbf_maxpt_jj_deta = fabs(outputAK5_PtSorted.at(0).Momentum_.Eta() - outputAK5_PtSorted.at(1).Momentum_.Eta()) ;
     if (fabs(outputAK5_PtSorted.at(0).Momentum_.Phi() - outputAK5_PtSorted.at(1).Momentum_.Phi()) < TMath::Pi())
       NewElectronTree->vbf_maxpt_jj_dphi = fabs(outputAK5_PtSorted.at(0).Momentum_.Phi() - outputAK5_PtSorted.at(1).Momentum_.Phi()) ;
     else 
       NewElectronTree->vbf_maxpt_jj_dphi = 2*TMath::Pi() - fabs(outputAK5_PtSorted.at(0).Momentum_.Phi() - outputAK5_PtSorted.at(1).Momentum_.Phi()) ;


     if(outputAK5_PtSorted.at(0).NameCollection_ == "JetPFCor" && outputAK5_PtSorted.at(1).NameCollection_ == "JetPFCor"){
       NewElectronTree->vbf_maxpt_jj_type = 1 ; /// both central 

      int nexcj = 0 , nexfj = 0; 
      std::vector<JetAK5>::const_iterator itVec = outputAK5_PtSorted.begin();
      for( ; itVec != outputAK5_PtSorted.end() ; itVec++){

	if(itVec->NameCollection_ == "JetPFCor") nexcj ++ ;
   	if(itVec->NameCollection_ == "JetPFCorVBFTag") nexfj ++ ;
   
      }

      NewElectronTree->vbf_maxpt_n_excj = nexcj ; // number of central jets 
      NewElectronTree->vbf_maxpt_n_exfj = nexfj ; // number of forward jets
     }

     else if(outputAK5_PtSorted.at(0).NameCollection_ == "JetPFCor" && outputAK5_PtSorted.at(1).NameCollection_ == "JetPFCorVBFTag"){
      
      NewElectronTree->vbf_maxpt_jj_type = 2 ;
      int nexcj = 0 , nexfj = 0;
      std::vector<JetAK5>::const_iterator itVec = outputAK5_PtSorted.begin();
      for( ; itVec != outputAK5_PtSorted.end() ; itVec++){

	if(itVec->NameCollection_ == "JetPFCor") nexcj ++ ;
   	if(itVec->NameCollection_ == "JetPFCorVBFTag") nexfj ++ ;
   
      }

      NewElectronTree->vbf_maxpt_n_excj = nexcj ;
      NewElectronTree->vbf_maxpt_n_exfj = nexfj ;
     }

     else if(outputAK5_PtSorted.at(0).NameCollection_ == "JetPFCorVBFTag" && outputAK5_PtSorted.at(1).NameCollection_ == "JetPFCor"){
      
      NewElectronTree->vbf_maxpt_jj_type = 3 ;
      int nexcj = 0 , nexfj = 0;
      std::vector<JetAK5>::const_iterator itVec = outputAK5_PtSorted.begin();
      for( ; itVec != outputAK5_PtSorted.end() ; itVec++){

	if(itVec->NameCollection_ == "JetPFCor") nexcj ++ ;
   	if(itVec->NameCollection_ == "JetPFCorVBFTag") nexfj ++ ;
   
      }

      NewElectronTree->vbf_maxpt_n_excj = nexcj ;
      NewElectronTree->vbf_maxpt_n_exfj = nexfj ;
     }

     else if(outputAK5_PtSorted.at(0).NameCollection_ == "JetPFCorVBFTag" && outputAK5_PtSorted.at(1).NameCollection_ == "JetPFCorVBFTag"){
      
      NewElectronTree->vbf_maxpt_jj_type = 4 ;
      int nexcj = 0 , nexfj = 0;
      std::vector<JetAK5>::const_iterator itVec = outputAK5_PtSorted.begin();
      for( ; itVec != outputAK5_PtSorted.end() ; itVec++){

	if(itVec->NameCollection_ == "JetPFCor") nexcj ++ ;
   	if(itVec->NameCollection_ == "JetPFCorVBFTag") nexfj ++ ;
   
      }

      NewElectronTree->vbf_maxpt_n_excj = nexcj ;
      NewElectronTree->vbf_maxpt_n_exfj = nexfj ;
     }
     else{ std::cerr<<" Something Wrong in MaxPt Jet Categorization "<<std::endl; continue ;}

     if(std::string(SelectionEvents->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEvents->GetXaxis()->SetBinLabel(nStep,"Jet MaxPt Category");
     if(std::string(SelectionEfficiency->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEfficiency->GetXaxis()->SetBinLabel(nStep,"Jet MaxPt Category");
     nstepEvents[nStep-1]++;
     nStep = 7;

     if( NewElectronTree->vbf_maxpt_jj_type < 0 || NewElectronTree->vbf_maxpt_n_excj < 0 || NewElectronTree->vbf_maxpt_n_exfj < 0 ) continue ;

     
     if(!numberJetBinGen.empty() && !outputGenAK5_PtSorted.empty()){
      if(numberJetBinGen.at(0)>=2 && outputGenAK5_PtSorted.size()>=2){

      TLorentzVector Total4VMaxPtGen = outputGenAK5_PtSorted.at(0).Momentum_ + outputGenAK5_PtSorted.at(1).Momentum_ ;
      NewElectronTree->vbf_maxpt_jj_e_gen   = Total4VMaxPtGen.E(); 
      NewElectronTree->vbf_maxpt_jj_pt_gen  = Total4VMaxPtGen.Pt(); 
      NewElectronTree->vbf_maxpt_jj_eta_gen = Total4VMaxPtGen.Eta(); 
      NewElectronTree->vbf_maxpt_jj_phi_gen = Total4VMaxPtGen.Phi(); 
      NewElectronTree->vbf_maxpt_jj_m_gen   = Total4VMaxPtGen.M(); 
      NewElectronTree->vbf_maxpt_j1_e_gen   = outputGenAK5_PtSorted.at(0).Momentum_.E();
      NewElectronTree->vbf_maxpt_j1_pt_gen  = outputGenAK5_PtSorted.at(0).Momentum_.Pt();
      NewElectronTree->vbf_maxpt_j1_eta_gen = outputGenAK5_PtSorted.at(0).Momentum_.Eta();
      NewElectronTree->vbf_maxpt_j1_phi_gen = outputGenAK5_PtSorted.at(0).Momentum_.Phi();
      NewElectronTree->vbf_maxpt_j1_m_gen   = outputGenAK5_PtSorted.at(0).Momentum_.M(); 
      NewElectronTree->vbf_maxpt_j2_e_gen   = outputGenAK5_PtSorted.at(1).Momentum_.E();
      NewElectronTree->vbf_maxpt_j2_pt_gen  = outputGenAK5_PtSorted.at(1).Momentum_.Pt();
      NewElectronTree->vbf_maxpt_j2_eta_gen = outputGenAK5_PtSorted.at(1).Momentum_.Eta();
      NewElectronTree->vbf_maxpt_j2_phi_gen = outputGenAK5_PtSorted.at(1).Momentum_.Phi();
      NewElectronTree->vbf_maxpt_j2_m_gen   = outputGenAK5_PtSorted.at(1).Momentum_.M();

      NewElectronTree->vbf_maxpt_jj_deta_gen = fabs(outputGenAK5_PtSorted.at(0).Momentum_.Eta() - outputGenAK5_PtSorted.at(1).Momentum_.Eta()) ;
      if (fabs(outputGenAK5_PtSorted.at(0).Momentum_.Phi() - outputGenAK5_PtSorted.at(1).Momentum_.Phi()) < TMath::Pi())
       NewElectronTree->vbf_maxpt_jj_dphi_gen = fabs(outputGenAK5_PtSorted.at(0).Momentum_.Phi() - outputGenAK5_PtSorted.at(1).Momentum_.Phi()) ;
      else 
       NewElectronTree->vbf_maxpt_jj_dphi_gen = 2*TMath::Pi() - fabs(outputGenAK5_PtSorted.at(0).Momentum_.Phi() - outputGenAK5_PtSorted.at(1).Momentum_.Phi()) ;

      NewElectronTree->vbf_maxpt_j1_bDiscriminatorSSVHE_gen = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorSSVHE")[outputGenAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_bDiscriminatorTCHE_gen  = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorTCHE")[outputGenAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_bDiscriminatorCSV_gen   = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorCSV")[outputGenAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_bDiscriminatorSSVHP_gen = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorSSVHP")[outputGenAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_bDiscriminatorTCHP_gen  = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorTCHP")[outputGenAK5_PtSorted.at(0).position_] ;

      NewElectronTree->vbf_maxpt_j2_bDiscriminatorSSVHE_gen = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorSSVHE")[outputGenAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_bDiscriminatorTCHE_gen  = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorTCHE")[outputGenAK5_PtSorted.at(1).position_] ; 
      NewElectronTree->vbf_maxpt_j2_bDiscriminatorCSV_gen   = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorCSV")[outputGenAK5_PtSorted.at(1).position_]  ;
      NewElectronTree->vbf_maxpt_j2_bDiscriminatorSSVHP_gen = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorSSVHP")[outputGenAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_bDiscriminatorTCHP_gen  = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorTCHP")[outputGenAK5_PtSorted.at(1).position_] ;
      }
      }    
     if(outputAK5_PtSorted.at(0).NameCollection_ == "JetPFCor") {
    
      NewElectronTree->vbf_maxpt_j1_QGLikelihood = ElectronTree->fReader->getFloat("JetPFCor_QGLikelihood")[outputAK5_PtSorted.at(0).position_] ;

      NewElectronTree->vbf_maxpt_j1_isPileUpLoose  = ElectronTree->fReader->getBool("JetPFCor_isPileUpJetLoose")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_isPileUpMedium = ElectronTree->fReader->getBool("JetPFCor_isPileUpJetMedium")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_isPileUpTight  = ElectronTree->fReader->getBool("JetPFCor_isPileUpJetTight")[outputAK5_PtSorted.at(0).position_] ;


      NewElectronTree->vbf_maxpt_j1_bDiscriminatorSSVHE = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorSSVHE")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_bDiscriminatorTCHE  = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorTCHE")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_bDiscriminatorCSV   = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorCSV")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_bDiscriminatorSSVHP = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorSSVHP")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_bDiscriminatorTCHP  = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorTCHP")[outputAK5_PtSorted.at(0).position_] ;

      NewElectronTree->vbf_maxpt_j1_ChargedHadronEnergy      = ElectronTree->fReader->getFloat("JetPFCor_ChargedHadronEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_ChargedHadronEnergyFrac  = ElectronTree->fReader->getFloat("JetPFCor_ChargedHadronEnergyFrac")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_NeutralHadronEnergy      = ElectronTree->fReader->getFloat("JetPFCor_NeutralHadronEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_NeutralHadronEnergyFrac  = ElectronTree->fReader->getFloat("JetPFCor_NeutralHadronEnergyFrac")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_ChargedEmEnergy          = ElectronTree->fReader->getFloat("JetPFCor_ChargedEmEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_ChargedEmEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCor_ChargedEmEnergyFrac")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_ChargedMuEnergy          = ElectronTree->fReader->getFloat("JetPFCor_ChargedMuEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_ChargedMuEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCor_ChargedMuEnergyFrac")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_NeutralEmEnergy          = ElectronTree->fReader->getFloat("JetPFCor_NeutralEmEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_NeutralEmEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCor_NeutralEmEnergyFrac")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_PhotonEnergy             = ElectronTree->fReader->getFloat("JetPFCor_PhotonEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_PhotonEnergyFraction     = ElectronTree->fReader->getFloat("JetPFCor_PhotonEnergyFraction")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_ElectronEnergy           = ElectronTree->fReader->getFloat("JetPFCor_ElectronEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_ElectronEnergyFraction   = ElectronTree->fReader->getFloat("JetPFCor_ElectronEnergyFraction")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_HFHadronEnergy           = ElectronTree->fReader->getFloat("JetPFCor_HFHadronEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_HFHadronEnergyFraction   = ElectronTree->fReader->getFloat("JetPFCor_HFHadronEnergyFraction")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_HFEMEnergy               = ElectronTree->fReader->getFloat("JetPFCor_HFEMEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_HFEMEnergyFraction       = ElectronTree->fReader->getFloat("JetPFCor_HFEMEnergyFraction")[outputAK5_PtSorted.at(0).position_] ;

      NewElectronTree->vbf_maxpt_j1_ChargedMultiplicity        = ElectronTree->fReader->getFloat("JetPFCor_ChargedMultiplicity")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_NeutralMultiplicity        = ElectronTree->fReader->getFloat("JetPFCor_NeutralMultiplicity")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_MuonMultiplicity           = ElectronTree->fReader->getFloat("JetPFCor_MuonMultiplicity")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_ChargedHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCor_ChargedHadronMultiplicity")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_NeutralHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCor_NeutralHadronMultiplicity")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_PhotonMultiplicity         = ElectronTree->fReader->getFloat("JetPFCor_PhotonMultiplicity")[outputAK5_PtSorted.at(0).position_] ; 
      NewElectronTree->vbf_maxpt_j1_ElectronMultiplicity       = ElectronTree->fReader->getFloat("JetPFCor_ElectronMultiplicity")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_HFHadronMultiplicity       = ElectronTree->fReader->getFloat("JetPFCor_HFHadronMultiplicity")[outputAK5_PtSorted.at(0).position_] ;
     }
     else if(outputAK5_PtSorted.at(0).NameCollection_ == "JetPFCorVBFTag" ){

      NewElectronTree->vbf_maxpt_j1_QGLikelihood = ElectronTree->fReader->getFloat("JetPFCorVBFTag_QGLikelihood")[outputAK5_PtSorted.at(0).position_] ;

      NewElectronTree->vbf_maxpt_j1_isPileUpLoose  = ElectronTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetLoose")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_isPileUpMedium = ElectronTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetMedium")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_isPileUpTight  = ElectronTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetTight")[outputAK5_PtSorted.at(0).position_] ;


      NewElectronTree->vbf_maxpt_j1_bDiscriminatorSSVHE = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorSSVHE")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_bDiscriminatorTCHE  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorTCHE")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_bDiscriminatorCSV   = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorCSV")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_bDiscriminatorSSVHP = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorSSVHP")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_bDiscriminatorTCHP  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorTCHP")[outputAK5_PtSorted.at(0).position_] ;

      NewElectronTree->vbf_maxpt_j1_ChargedHadronEnergy      = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_ChargedHadronEnergyFrac  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronEnergyFrac")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_NeutralHadronEnergy      = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_NeutralHadronEnergyFrac  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronEnergyFrac")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_ChargedEmEnergy          = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedEmEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_ChargedEmEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedEmEnergyFrac")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_ChargedMuEnergy          = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedMuEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_ChargedMuEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedMuEnergyFrac")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_NeutralEmEnergy          = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralEmEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_NeutralEmEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralEmEnergyFrac")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_PhotonEnergy             = ElectronTree->fReader->getFloat("JetPFCorVBFTag_PhotonEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_PhotonEnergyFraction     = ElectronTree->fReader->getFloat("JetPFCorVBFTag_PhotonEnergyFraction")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_ElectronEnergy           = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ElectronEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_ElectronEnergyFraction   = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ElectronEnergyFraction")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_HFHadronEnergy           = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFHadronEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_HFHadronEnergyFraction   = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFHadronEnergyFraction")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_HFEMEnergy               = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFEMEnergy")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_HFEMEnergyFraction       = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFEMEnergyFraction")[outputAK5_PtSorted.at(0).position_] ;

      NewElectronTree->vbf_maxpt_j1_ChargedMultiplicity        = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedMultiplicity")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_NeutralMultiplicity        = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralMultiplicity")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_MuonMultiplicity           = ElectronTree->fReader->getFloat("JetPFCorVBFTag_MuonMultiplicity")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_ChargedHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronMultiplicity")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_NeutralHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronMultiplicity")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_PhotonMultiplicity         = ElectronTree->fReader->getFloat("JetPFCorVBFTag_PhotonMultiplicity")[outputAK5_PtSorted.at(0).position_] ; 
      NewElectronTree->vbf_maxpt_j1_ElectronMultiplicity       = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ElectronMultiplicity")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_HFHadronMultiplicity       = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFHadronMultiplicity")[outputAK5_PtSorted.at(0).position_] ;  

     }
     else { std::cerr<<" problem with High pT Jet Name Collection "<<std::endl; continue ; }

     if(outputAK5_PtSorted.at(1).NameCollection_ == "JetPFCor") {
    
      NewElectronTree->vbf_maxpt_j2_QGLikelihood = ElectronTree->fReader->getFloat("JetPFCor_QGLikelihood")[outputAK5_PtSorted.at(1).position_] ;

      NewElectronTree->vbf_maxpt_j2_isPileUpLoose  = ElectronTree->fReader->getBool("JetPFCor_isPileUpJetLoose")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_isPileUpMedium = ElectronTree->fReader->getBool("JetPFCor_isPileUpJetMedium")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_isPileUpTight  = ElectronTree->fReader->getBool("JetPFCor_isPileUpJetTight")[outputAK5_PtSorted.at(1).position_] ;


      NewElectronTree->vbf_maxpt_j2_bDiscriminatorSSVHE = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorSSVHE")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_bDiscriminatorTCHE  = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorTCHE")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_bDiscriminatorCSV   = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorCSV")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_bDiscriminatorSSVHP = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorSSVHP")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_bDiscriminatorTCHP  = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorTCHP")[outputAK5_PtSorted.at(1).position_] ;

      NewElectronTree->vbf_maxpt_j2_ChargedHadronEnergy      = ElectronTree->fReader->getFloat("JetPFCor_ChargedHadronEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_ChargedHadronEnergyFrac  = ElectronTree->fReader->getFloat("JetPFCor_ChargedHadronEnergyFrac")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_NeutralHadronEnergy      = ElectronTree->fReader->getFloat("JetPFCor_NeutralHadronEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_NeutralHadronEnergyFrac  = ElectronTree->fReader->getFloat("JetPFCor_NeutralHadronEnergyFrac")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_ChargedEmEnergy          = ElectronTree->fReader->getFloat("JetPFCor_ChargedEmEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_ChargedEmEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCor_ChargedEmEnergyFrac")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_ChargedMuEnergy          = ElectronTree->fReader->getFloat("JetPFCor_ChargedMuEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_ChargedMuEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCor_ChargedMuEnergyFrac")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_NeutralEmEnergy          = ElectronTree->fReader->getFloat("JetPFCor_NeutralEmEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_NeutralEmEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCor_NeutralEmEnergyFrac")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_PhotonEnergy             = ElectronTree->fReader->getFloat("JetPFCor_PhotonEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_PhotonEnergyFraction     = ElectronTree->fReader->getFloat("JetPFCor_PhotonEnergyFraction")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_ElectronEnergy           = ElectronTree->fReader->getFloat("JetPFCor_ElectronEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_ElectronEnergyFraction   = ElectronTree->fReader->getFloat("JetPFCor_ElectronEnergyFraction")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_HFHadronEnergy           = ElectronTree->fReader->getFloat("JetPFCor_HFHadronEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_HFHadronEnergyFraction   = ElectronTree->fReader->getFloat("JetPFCor_HFHadronEnergyFraction")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_HFEMEnergy               = ElectronTree->fReader->getFloat("JetPFCor_HFEMEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_HFEMEnergyFraction       = ElectronTree->fReader->getFloat("JetPFCor_HFEMEnergyFraction")[outputAK5_PtSorted.at(1).position_] ;

      NewElectronTree->vbf_maxpt_j2_ChargedMultiplicity        = ElectronTree->fReader->getFloat("JetPFCor_ChargedMultiplicity")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_NeutralMultiplicity        = ElectronTree->fReader->getFloat("JetPFCor_NeutralMultiplicity")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_MuonMultiplicity           = ElectronTree->fReader->getFloat("JetPFCor_MuonMultiplicity")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_ChargedHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCor_ChargedHadronMultiplicity")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_NeutralHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCor_NeutralHadronMultiplicity")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_PhotonMultiplicity         = ElectronTree->fReader->getFloat("JetPFCor_PhotonMultiplicity")[outputAK5_PtSorted.at(1).position_] ; 
      NewElectronTree->vbf_maxpt_j2_ElectronMultiplicity       = ElectronTree->fReader->getFloat("JetPFCor_ElectronMultiplicity")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_HFHadronMultiplicity       = ElectronTree->fReader->getFloat("JetPFCor_HFHadronMultiplicity")[outputAK5_PtSorted.at(1).position_] ;
     }
     else if(outputAK5_PtSorted.at(1).NameCollection_ == "JetPFCorVBFTag" ){

      NewElectronTree->vbf_maxpt_j2_QGLikelihood = ElectronTree->fReader->getFloat("JetPFCorVBFTag_QGLikelihood")[outputAK5_PtSorted.at(1).position_] ;

      NewElectronTree->vbf_maxpt_j2_isPileUpLoose  = ElectronTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetLoose")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_isPileUpMedium = ElectronTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetMedium")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_isPileUpTight  = ElectronTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetTight")[outputAK5_PtSorted.at(1).position_] ;


      NewElectronTree->vbf_maxpt_j2_bDiscriminatorSSVHE = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorSSVHE")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_bDiscriminatorTCHE  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorTCHE")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_bDiscriminatorCSV   = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorCSV")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_bDiscriminatorSSVHP = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorSSVHP")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_bDiscriminatorTCHP  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorTCHP")[outputAK5_PtSorted.at(1).position_] ;

      NewElectronTree->vbf_maxpt_j2_ChargedHadronEnergy      = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_ChargedHadronEnergyFrac  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronEnergyFrac")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_NeutralHadronEnergy      = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_NeutralHadronEnergyFrac  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronEnergyFrac")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_ChargedEmEnergy          = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedEmEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_ChargedEmEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedEmEnergyFrac")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_ChargedMuEnergy          = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedMuEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_ChargedMuEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedMuEnergyFrac")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_NeutralEmEnergy          = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralEmEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_NeutralEmEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralEmEnergyFrac")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_PhotonEnergy             = ElectronTree->fReader->getFloat("JetPFCorVBFTag_PhotonEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_PhotonEnergyFraction     = ElectronTree->fReader->getFloat("JetPFCorVBFTag_PhotonEnergyFraction")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_ElectronEnergy           = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ElectronEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_ElectronEnergyFraction   = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ElectronEnergyFraction")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_HFHadronEnergy           = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFHadronEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_HFHadronEnergyFraction   = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFHadronEnergyFraction")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_HFEMEnergy               = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFEMEnergy")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_HFEMEnergyFraction       = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFEMEnergyFraction")[outputAK5_PtSorted.at(1).position_] ;

      NewElectronTree->vbf_maxpt_j2_ChargedMultiplicity        = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedMultiplicity")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_NeutralMultiplicity        = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralMultiplicity")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_MuonMultiplicity           = ElectronTree->fReader->getFloat("JetPFCorVBFTag_MuonMultiplicity")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_ChargedHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronMultiplicity")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_NeutralHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronMultiplicity")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_PhotonMultiplicity         = ElectronTree->fReader->getFloat("JetPFCorVBFTag_PhotonMultiplicity")[outputAK5_PtSorted.at(1).position_] ; 
      NewElectronTree->vbf_maxpt_j2_ElectronMultiplicity       = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ElectronMultiplicity")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_HFHadronMultiplicity       = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFHadronMultiplicity")[outputAK5_PtSorted.at(1).position_] ;  

     }
     else { std::cerr<<" problem with High pT Jet Name Collection "<<std::endl; continue ; }

     /////////////////////////////////////////////////////////////////////////////////////////////////
     // Fill Information for Max Deta Pair of vbf tag jets
     //////////////////////////////////////////////////////////////////////////////////////////////////
   
     TLorentzVector Total4VMaxDeta = outputAK5_DEtaSorted.at(0).Momentum_ + outputAK5_DEtaSorted.at(1).Momentum_ ;
    
     NewElectronTree->vbf_maxDeta_jj_e   = Total4VMaxDeta.E(); 
     NewElectronTree->vbf_maxDeta_jj_pt  = Total4VMaxDeta.Pt(); 
     NewElectronTree->vbf_maxDeta_jj_eta = Total4VMaxDeta.Eta(); 
     NewElectronTree->vbf_maxDeta_jj_phi = Total4VMaxDeta.Phi(); 
     NewElectronTree->vbf_maxDeta_jj_m   = Total4VMaxDeta.M(); 
 
     NewElectronTree->vbf_maxDeta_j1_e   = outputAK5_DEtaSorted.at(0).Momentum_.E();
     NewElectronTree->vbf_maxDeta_j1_pt  = outputAK5_DEtaSorted.at(0).Momentum_.Pt();
     NewElectronTree->vbf_maxDeta_j1_eta = outputAK5_DEtaSorted.at(0).Momentum_.Eta();
     NewElectronTree->vbf_maxDeta_j1_phi = outputAK5_DEtaSorted.at(0).Momentum_.Phi();
     NewElectronTree->vbf_maxDeta_j1_m   = outputAK5_DEtaSorted.at(0).Momentum_.M();

     NewElectronTree->vbf_maxDeta_j2_e   = outputAK5_DEtaSorted.at(1).Momentum_.E();
     NewElectronTree->vbf_maxDeta_j2_pt  = outputAK5_DEtaSorted.at(1).Momentum_.Pt();
     NewElectronTree->vbf_maxDeta_j2_eta = outputAK5_DEtaSorted.at(1).Momentum_.Eta();
     NewElectronTree->vbf_maxDeta_j2_phi = outputAK5_DEtaSorted.at(1).Momentum_.Phi();
     NewElectronTree->vbf_maxDeta_j2_m   = outputAK5_DEtaSorted.at(1).Momentum_.M();
   
     NewElectronTree->vbf_maxDeta_jj_deta = fabs(outputAK5_DEtaSorted.at(0).Momentum_.Eta() - outputAK5_DEtaSorted.at(1).Momentum_.Eta()) ;
     if (fabs(outputAK5_DEtaSorted.at(0).Momentum_.Phi() - outputAK5_DEtaSorted.at(1).Momentum_.Phi()) < 3.14)
      NewElectronTree->vbf_maxDeta_jj_dphi = fabs(outputAK5_DEtaSorted.at(0).Momentum_.Phi() - outputAK5_DEtaSorted.at(1).Momentum_.Phi()) ;
     else  
      NewElectronTree->vbf_maxDeta_jj_dphi = 6.28 - fabs(outputAK5_DEtaSorted.at(0).Momentum_.Phi() - outputAK5_DEtaSorted.at(1).Momentum_.Phi()) ;

     if(outputAK5_DEtaSorted.at(0).NameCollection_ == "JetPFCor" && outputAK5_DEtaSorted.at(1).NameCollection_ == "JetPFCor"){
      
      NewElectronTree->vbf_maxDeta_jj_type = 1 ;

      int nexcj = 0 , nexfj = 0; 
      std::vector<JetAK5>::const_iterator itVec = outputAK5_DEtaSorted.begin();
      for( ; itVec != outputAK5_DEtaSorted.end() ; itVec++){

	if(itVec->NameCollection_ == "JetPFCor") nexcj ++ ;
   	if(itVec->NameCollection_ == "JetPFCorVBFTag") nexfj ++ ;
   
      }

      NewElectronTree->vbf_maxDeta_n_excj = nexcj ;
      NewElectronTree->vbf_maxDeta_n_exfj = nexfj ;
     }

     else if(outputAK5_DEtaSorted.at(0).NameCollection_ == "JetPFCor" && outputAK5_DEtaSorted.at(1).NameCollection_ == "JetPFCorVBFTag"){
      
      NewElectronTree->vbf_maxDeta_jj_type = 2 ;

      int nexcj = 0 , nexfj = 0;
      std::vector<JetAK5>::const_iterator itVec = outputAK5_DEtaSorted.begin();
      for( ; itVec != outputAK5_DEtaSorted.end() ; itVec++){

	if(itVec->NameCollection_ == "JetPFCor") nexcj ++ ;
   	if(itVec->NameCollection_ == "JetPFCorVBFTag") nexfj ++ ;
   
      }

      NewElectronTree->vbf_maxDeta_n_excj = nexcj ;
      NewElectronTree->vbf_maxDeta_n_exfj = nexfj ;
     }

     else if(outputAK5_DEtaSorted.at(0).NameCollection_ == "JetPFCorVBFTag" && outputAK5_DEtaSorted.at(1).NameCollection_ == "JetPFCor"){
      
      NewElectronTree->vbf_maxDeta_jj_type = 3 ;

      int nexcj = 0 , nexfj = 0;
      std::vector<JetAK5>::const_iterator itVec = outputAK5_DEtaSorted.begin();
      for( ; itVec != outputAK5_DEtaSorted.end() ; itVec++){

	if(itVec->NameCollection_ == "JetPFCor") nexcj ++ ;
   	if(itVec->NameCollection_ == "JetPFCorVBFTag") nexfj ++ ;
   
      }

      NewElectronTree->vbf_maxDeta_n_excj = nexcj ;
      NewElectronTree->vbf_maxDeta_n_exfj = nexfj ;
     }

     else if(outputAK5_DEtaSorted.at(0).NameCollection_ == "JetPFCorVBFTag" && outputAK5_DEtaSorted.at(1).NameCollection_ == "JetPFCorVBFTag"){
      
      NewElectronTree->vbf_maxDeta_jj_type = 4 ;

      int nexcj = 0 , nexfj = 0;
      std::vector<JetAK5>::const_iterator itVec = outputAK5_DEtaSorted.begin();
      for( ; itVec != outputAK5_DEtaSorted.end() ; itVec++){

	if(itVec->NameCollection_ == "JetPFCor") nexcj ++ ;
   	if(itVec->NameCollection_ == "JetPFCorVBFTag") nexfj ++ ;
   
      }

      NewElectronTree->vbf_maxDeta_n_excj = nexcj ;
      NewElectronTree->vbf_maxDeta_n_exfj = nexfj ;
     }
     else{ std::cerr<<" Something Wrong in MaxDeta Jet Categorization "<<std::endl; continue ;}

     if(std::string(SelectionEvents->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEvents->GetXaxis()->SetBinLabel(nStep,"Jet MaxDeta Category");
     if(std::string(SelectionEfficiency->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEfficiency->GetXaxis()->SetBinLabel(nStep,"Jet MaxDeta Category");
     nstepEvents[nStep-1]++;
     nStep = 8;

     
     if(!numberJetBinGen.empty() && !outputGenAK5_DEtaSorted.empty()){
      if(numberJetBinGen.at(0)>=2 && outputGenAK5_DEtaSorted.size()>=2){

      TLorentzVector Total4VMaxDetaGen = outputGenAK5_DEtaSorted.at(0).Momentum_ + outputGenAK5_DEtaSorted.at(1).Momentum_ ;
      NewElectronTree->vbf_maxDeta_jj_e_gen   = Total4VMaxDetaGen.E(); 
      NewElectronTree->vbf_maxDeta_jj_pt_gen  = Total4VMaxDetaGen.Pt(); 
      NewElectronTree->vbf_maxDeta_jj_eta_gen = Total4VMaxDetaGen.Eta(); 
      NewElectronTree->vbf_maxDeta_jj_phi_gen = Total4VMaxDetaGen.Phi(); 
      NewElectronTree->vbf_maxDeta_jj_m_gen   = Total4VMaxDetaGen.M(); 
      NewElectronTree->vbf_maxDeta_j1_e_gen   = outputGenAK5_DEtaSorted.at(0).Momentum_.E();
      NewElectronTree->vbf_maxDeta_j1_pt_gen  = outputGenAK5_DEtaSorted.at(0).Momentum_.Pt();
      NewElectronTree->vbf_maxDeta_j1_eta_gen = outputGenAK5_DEtaSorted.at(0).Momentum_.Eta();
      NewElectronTree->vbf_maxDeta_j1_phi_gen = outputGenAK5_DEtaSorted.at(0).Momentum_.Phi();
      NewElectronTree->vbf_maxDeta_j1_m_gen   = outputGenAK5_DEtaSorted.at(0).Momentum_.M(); 
      NewElectronTree->vbf_maxDeta_j2_e_gen   = outputGenAK5_DEtaSorted.at(1).Momentum_.E();
      NewElectronTree->vbf_maxDeta_j2_pt_gen  = outputGenAK5_DEtaSorted.at(1).Momentum_.Pt();
      NewElectronTree->vbf_maxDeta_j2_eta_gen = outputGenAK5_DEtaSorted.at(1).Momentum_.Eta();
      NewElectronTree->vbf_maxDeta_j2_phi_gen = outputGenAK5_DEtaSorted.at(1).Momentum_.Phi();
      NewElectronTree->vbf_maxDeta_j2_m_gen   = outputGenAK5_DEtaSorted.at(1).Momentum_.M();

      NewElectronTree->vbf_maxDeta_jj_deta_gen = fabs(outputGenAK5_DEtaSorted.at(0).Momentum_.Eta() - outputGenAK5_DEtaSorted.at(1).Momentum_.Eta()) ;
      if (fabs(outputGenAK5_DEtaSorted.at(0).Momentum_.Phi() - outputGenAK5_DEtaSorted.at(1).Momentum_.Phi()) < TMath::Pi())
       NewElectronTree->vbf_maxDeta_jj_dphi_gen = fabs(outputGenAK5_DEtaSorted.at(0).Momentum_.Phi() - outputGenAK5_DEtaSorted.at(1).Momentum_.Phi()) ;
      else 
       NewElectronTree->vbf_maxDeta_jj_dphi_gen = 2*TMath::Pi() - fabs(outputGenAK5_DEtaSorted.at(0).Momentum_.Phi() - outputGenAK5_DEtaSorted.at(1).Momentum_.Phi()) ;

      NewElectronTree->vbf_maxDeta_j1_bDiscriminatorSSVHE_gen = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorSSVHE")[outputGenAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_bDiscriminatorTCHE_gen  = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorTCHE")[outputGenAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_bDiscriminatorCSV_gen   = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorCSV")[outputGenAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_bDiscriminatorSSVHP_gen = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorSSVHP")[outputGenAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_bDiscriminatorTCHP_gen  = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorTCHP")[outputGenAK5_DEtaSorted.at(0).position_] ;

      NewElectronTree->vbf_maxDeta_j2_bDiscriminatorSSVHE_gen = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorSSVHE")[outputGenAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_bDiscriminatorTCHE_gen  = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorTCHE")[outputGenAK5_DEtaSorted.at(1).position_] ; 
      NewElectronTree->vbf_maxDeta_j2_bDiscriminatorCSV_gen   = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorCSV")[outputGenAK5_DEtaSorted.at(1).position_]  ;
      NewElectronTree->vbf_maxDeta_j2_bDiscriminatorSSVHP_gen = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorSSVHP")[outputGenAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_bDiscriminatorTCHP_gen  = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorTCHP")[outputGenAK5_DEtaSorted.at(1).position_] ;
     }
    }

    if(outputAK5_DEtaSorted.at(0).NameCollection_ == "JetPFCor") {
    
      NewElectronTree->vbf_maxDeta_j1_QGLikelihood = ElectronTree->fReader->getFloat("JetPFCor_QGLikelihood")[outputAK5_DEtaSorted.at(0).position_] ;

      NewElectronTree->vbf_maxDeta_j1_isPileUpLoose  = ElectronTree->fReader->getBool("JetPFCor_isPileUpJetLoose")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_isPileUpMedium = ElectronTree->fReader->getBool("JetPFCor_isPileUpJetMedium")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_isPileUpTight  = ElectronTree->fReader->getBool("JetPFCor_isPileUpJetTight")[outputAK5_DEtaSorted.at(0).position_] ;


      NewElectronTree->vbf_maxDeta_j1_bDiscriminatorSSVHE = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorSSVHE")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_bDiscriminatorTCHE  = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorTCHE")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_bDiscriminatorCSV   = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorCSV")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_bDiscriminatorSSVHP = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorSSVHP")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_bDiscriminatorTCHP  = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorTCHP")[outputAK5_DEtaSorted.at(0).position_] ;

      NewElectronTree->vbf_maxDeta_j1_ChargedHadronEnergy      = ElectronTree->fReader->getFloat("JetPFCor_ChargedHadronEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_ChargedHadronEnergyFrac  = ElectronTree->fReader->getFloat("JetPFCor_ChargedHadronEnergyFrac")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_NeutralHadronEnergy      = ElectronTree->fReader->getFloat("JetPFCor_NeutralHadronEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_NeutralHadronEnergyFrac  = ElectronTree->fReader->getFloat("JetPFCor_NeutralHadronEnergyFrac")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_ChargedEmEnergy          = ElectronTree->fReader->getFloat("JetPFCor_ChargedEmEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_ChargedEmEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCor_ChargedEmEnergyFrac")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_ChargedMuEnergy          = ElectronTree->fReader->getFloat("JetPFCor_ChargedMuEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_ChargedMuEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCor_ChargedMuEnergyFrac")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_NeutralEmEnergy          = ElectronTree->fReader->getFloat("JetPFCor_NeutralEmEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_NeutralEmEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCor_NeutralEmEnergyFrac")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_PhotonEnergy             = ElectronTree->fReader->getFloat("JetPFCor_PhotonEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_PhotonEnergyFraction     = ElectronTree->fReader->getFloat("JetPFCor_PhotonEnergyFraction")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_ElectronEnergy           = ElectronTree->fReader->getFloat("JetPFCor_ElectronEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_ElectronEnergyFraction   = ElectronTree->fReader->getFloat("JetPFCor_ElectronEnergyFraction")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_HFHadronEnergy           = ElectronTree->fReader->getFloat("JetPFCor_HFHadronEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_HFHadronEnergyFraction   = ElectronTree->fReader->getFloat("JetPFCor_HFHadronEnergyFraction")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_HFEMEnergy               = ElectronTree->fReader->getFloat("JetPFCor_HFEMEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_HFEMEnergyFraction       = ElectronTree->fReader->getFloat("JetPFCor_HFEMEnergyFraction")[outputAK5_DEtaSorted.at(0).position_] ;

      NewElectronTree->vbf_maxDeta_j1_ChargedMultiplicity        = ElectronTree->fReader->getFloat("JetPFCor_ChargedMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_NeutralMultiplicity        = ElectronTree->fReader->getFloat("JetPFCor_NeutralMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_MuonMultiplicity           = ElectronTree->fReader->getFloat("JetPFCor_MuonMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_ChargedHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCor_ChargedHadronMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_NeutralHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCor_NeutralHadronMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_PhotonMultiplicity         = ElectronTree->fReader->getFloat("JetPFCor_PhotonMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ; 
      NewElectronTree->vbf_maxDeta_j1_ElectronMultiplicity       = ElectronTree->fReader->getFloat("JetPFCor_ElectronMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_HFHadronMultiplicity       = ElectronTree->fReader->getFloat("JetPFCor_HFHadronMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ;
     }
     else if(outputAK5_DEtaSorted.at(0).NameCollection_ == "JetPFCorVBFTag" ){

      NewElectronTree->vbf_maxDeta_j1_QGLikelihood = ElectronTree->fReader->getFloat("JetPFCorVBFTag_QGLikelihood")[outputAK5_DEtaSorted.at(0).position_] ;

      NewElectronTree->vbf_maxDeta_j1_isPileUpLoose  = ElectronTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetLoose")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_isPileUpMedium = ElectronTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetMedium")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_isPileUpTight  = ElectronTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetTight")[outputAK5_DEtaSorted.at(0).position_] ;


      NewElectronTree->vbf_maxDeta_j1_bDiscriminatorSSVHE = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorSSVHE")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_bDiscriminatorTCHE  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorTCHE")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_bDiscriminatorCSV   = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorCSV")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_bDiscriminatorSSVHP = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorSSVHP")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_bDiscriminatorTCHP  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorTCHP")[outputAK5_DEtaSorted.at(0).position_] ;

      NewElectronTree->vbf_maxDeta_j1_ChargedHadronEnergy      = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_ChargedHadronEnergyFrac  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronEnergyFrac")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_NeutralHadronEnergy      = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_NeutralHadronEnergyFrac  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronEnergyFrac")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_ChargedEmEnergy          = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedEmEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_ChargedEmEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedEmEnergyFrac")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_ChargedMuEnergy          = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedMuEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_ChargedMuEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedMuEnergyFrac")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_NeutralEmEnergy          = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralEmEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_NeutralEmEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralEmEnergyFrac")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_PhotonEnergy             = ElectronTree->fReader->getFloat("JetPFCorVBFTag_PhotonEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_PhotonEnergyFraction     = ElectronTree->fReader->getFloat("JetPFCorVBFTag_PhotonEnergyFraction")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_ElectronEnergy           = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ElectronEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_ElectronEnergyFraction   = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ElectronEnergyFraction")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_HFHadronEnergy           = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFHadronEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_HFHadronEnergyFraction   = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFHadronEnergyFraction")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_HFEMEnergy               = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFEMEnergy")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_HFEMEnergyFraction       = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFEMEnergyFraction")[outputAK5_DEtaSorted.at(0).position_] ;

      NewElectronTree->vbf_maxDeta_j1_ChargedMultiplicity        = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_NeutralMultiplicity        = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_MuonMultiplicity           = ElectronTree->fReader->getFloat("JetPFCorVBFTag_MuonMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_ChargedHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_NeutralHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_PhotonMultiplicity         = ElectronTree->fReader->getFloat("JetPFCorVBFTag_PhotonMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ; 
      NewElectronTree->vbf_maxDeta_j1_ElectronMultiplicity       = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ElectronMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_HFHadronMultiplicity       = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFHadronMultiplicity")[outputAK5_DEtaSorted.at(0).position_] ;  

     }
     else { std::cerr<<" problem with High Deta Jet Name Collection "<<std::endl; continue ; }

     if(outputAK5_DEtaSorted.at(1).NameCollection_ == "JetPFCor") {
    
      NewElectronTree->vbf_maxDeta_j2_QGLikelihood = ElectronTree->fReader->getFloat("JetPFCor_QGLikelihood")[outputAK5_DEtaSorted.at(1).position_] ;

      NewElectronTree->vbf_maxDeta_j2_isPileUpLoose  = ElectronTree->fReader->getBool("JetPFCor_isPileUpJetLoose")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_isPileUpMedium = ElectronTree->fReader->getBool("JetPFCor_isPileUpJetMedium")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_isPileUpTight  = ElectronTree->fReader->getBool("JetPFCor_isPileUpJetTight")[outputAK5_DEtaSorted.at(1).position_] ;


      NewElectronTree->vbf_maxDeta_j2_bDiscriminatorSSVHE = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorSSVHE")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_bDiscriminatorTCHE  = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorTCHE")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_bDiscriminatorCSV   = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorCSV")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_bDiscriminatorSSVHP = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorSSVHP")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_bDiscriminatorTCHP  = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorTCHP")[outputAK5_DEtaSorted.at(1).position_] ;

      NewElectronTree->vbf_maxDeta_j2_ChargedHadronEnergy      = ElectronTree->fReader->getFloat("JetPFCor_ChargedHadronEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_ChargedHadronEnergyFrac  = ElectronTree->fReader->getFloat("JetPFCor_ChargedHadronEnergyFrac")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_NeutralHadronEnergy      = ElectronTree->fReader->getFloat("JetPFCor_NeutralHadronEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_NeutralHadronEnergyFrac  = ElectronTree->fReader->getFloat("JetPFCor_NeutralHadronEnergyFrac")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_ChargedEmEnergy          = ElectronTree->fReader->getFloat("JetPFCor_ChargedEmEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_ChargedEmEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCor_ChargedEmEnergyFrac")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_ChargedMuEnergy          = ElectronTree->fReader->getFloat("JetPFCor_ChargedMuEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_ChargedMuEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCor_ChargedMuEnergyFrac")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_NeutralEmEnergy          = ElectronTree->fReader->getFloat("JetPFCor_NeutralEmEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_NeutralEmEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCor_NeutralEmEnergyFrac")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_PhotonEnergy             = ElectronTree->fReader->getFloat("JetPFCor_PhotonEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_PhotonEnergyFraction     = ElectronTree->fReader->getFloat("JetPFCor_PhotonEnergyFraction")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_ElectronEnergy           = ElectronTree->fReader->getFloat("JetPFCor_ElectronEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_ElectronEnergyFraction   = ElectronTree->fReader->getFloat("JetPFCor_ElectronEnergyFraction")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_HFHadronEnergy           = ElectronTree->fReader->getFloat("JetPFCor_HFHadronEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_HFHadronEnergyFraction   = ElectronTree->fReader->getFloat("JetPFCor_HFHadronEnergyFraction")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_HFEMEnergy               = ElectronTree->fReader->getFloat("JetPFCor_HFEMEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_HFEMEnergyFraction       = ElectronTree->fReader->getFloat("JetPFCor_HFEMEnergyFraction")[outputAK5_DEtaSorted.at(1).position_] ;

      NewElectronTree->vbf_maxDeta_j2_ChargedMultiplicity        = ElectronTree->fReader->getFloat("JetPFCor_ChargedMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_NeutralMultiplicity        = ElectronTree->fReader->getFloat("JetPFCor_NeutralMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_MuonMultiplicity           = ElectronTree->fReader->getFloat("JetPFCor_MuonMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_ChargedHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCor_ChargedHadronMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_NeutralHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCor_NeutralHadronMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_PhotonMultiplicity         = ElectronTree->fReader->getFloat("JetPFCor_PhotonMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ; 
      NewElectronTree->vbf_maxDeta_j2_ElectronMultiplicity       = ElectronTree->fReader->getFloat("JetPFCor_ElectronMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_HFHadronMultiplicity       = ElectronTree->fReader->getFloat("JetPFCor_HFHadronMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ;
     }
     else if(outputAK5_DEtaSorted.at(1).NameCollection_ == "JetPFCorVBFTag" ){

      NewElectronTree->vbf_maxDeta_j2_QGLikelihood = ElectronTree->fReader->getFloat("JetPFCorVBFTag_QGLikelihood")[outputAK5_DEtaSorted.at(1).position_] ;

      NewElectronTree->vbf_maxDeta_j2_isPileUpLoose = ElectronTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetLoose")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_isPileUpMedium = ElectronTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetMedium")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_isPileUpTight = ElectronTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetTight")[outputAK5_DEtaSorted.at(1).position_] ;


      NewElectronTree->vbf_maxDeta_j2_bDiscriminatorSSVHE = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorSSVHE")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_bDiscriminatorTCHE  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorTCHE")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_bDiscriminatorCSV   = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorCSV")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_bDiscriminatorSSVHP = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorSSVHP")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_bDiscriminatorTCHP  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorTCHP")[outputAK5_DEtaSorted.at(1).position_] ;

      NewElectronTree->vbf_maxDeta_j2_ChargedHadronEnergy      = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_ChargedHadronEnergyFrac  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronEnergyFrac")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_NeutralHadronEnergy      = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_NeutralHadronEnergyFrac  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronEnergyFrac")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_ChargedEmEnergy          = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedEmEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_ChargedEmEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedEmEnergyFrac")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_ChargedMuEnergy          = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedMuEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_ChargedMuEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedMuEnergyFrac")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_NeutralEmEnergy          = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralEmEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_NeutralEmEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralEmEnergyFrac")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_PhotonEnergy             = ElectronTree->fReader->getFloat("JetPFCorVBFTag_PhotonEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_PhotonEnergyFraction     = ElectronTree->fReader->getFloat("JetPFCorVBFTag_PhotonEnergyFraction")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_ElectronEnergy           = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ElectronEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_ElectronEnergyFraction   = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ElectronEnergyFraction")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_HFHadronEnergy           = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFHadronEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_HFHadronEnergyFraction   = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFHadronEnergyFraction")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_HFEMEnergy               = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFEMEnergy")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_HFEMEnergyFraction       = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFEMEnergyFraction")[outputAK5_DEtaSorted.at(1).position_] ;

      NewElectronTree->vbf_maxDeta_j2_ChargedMultiplicity  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_NeutralMultiplicity  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_MuonMultiplicity     = ElectronTree->fReader->getFloat("JetPFCorVBFTag_MuonMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_ChargedHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_NeutralHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_PhotonMultiplicity         = ElectronTree->fReader->getFloat("JetPFCorVBFTag_PhotonMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ; 
      NewElectronTree->vbf_maxDeta_j2_ElectronMultiplicity       = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ElectronMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_HFHadronMultiplicity       = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFHadronMultiplicity")[outputAK5_DEtaSorted.at(1).position_] ;  

     }
     else { std::cerr<<" problem with High Deta Jet Name Collection "<<std::endl; continue ; }
    
     /////////////////////////////////////////////////////////////////////////////////////////////////
     // Fill Information for Max Mjj Pair of vbf tag jets
     //////////////////////////////////////////////////////////////////////////////////////////////////

     TLorentzVector Total4VMaxMjj = outputAK5_MjjSorted.at(0).Momentum_ + outputAK5_MjjSorted.at(1).Momentum_ ;
    
     NewElectronTree->vbf_maxMjj_jj_e   = Total4VMaxMjj.E(); 
     NewElectronTree->vbf_maxMjj_jj_pt  = Total4VMaxMjj.Pt(); 
     NewElectronTree->vbf_maxMjj_jj_eta = Total4VMaxMjj.Eta(); 
     NewElectronTree->vbf_maxMjj_jj_phi = Total4VMaxMjj.Phi(); 
     NewElectronTree->vbf_maxMjj_jj_m   = Total4VMaxMjj.M(); 
 
     NewElectronTree->vbf_maxMjj_j1_e   = outputAK5_MjjSorted.at(0).Momentum_.E();
     NewElectronTree->vbf_maxMjj_j1_pt  = outputAK5_MjjSorted.at(0).Momentum_.Pt();
     NewElectronTree->vbf_maxMjj_j1_eta = outputAK5_MjjSorted.at(0).Momentum_.Eta();
     NewElectronTree->vbf_maxMjj_j1_phi = outputAK5_MjjSorted.at(0).Momentum_.Phi();
     NewElectronTree->vbf_maxMjj_j1_m   = outputAK5_MjjSorted.at(0).Momentum_.M();

     NewElectronTree->vbf_maxMjj_j2_e   = outputAK5_MjjSorted.at(1).Momentum_.E();
     NewElectronTree->vbf_maxMjj_j2_pt  = outputAK5_MjjSorted.at(1).Momentum_.Pt();
     NewElectronTree->vbf_maxMjj_j2_eta = outputAK5_MjjSorted.at(1).Momentum_.Eta();
     NewElectronTree->vbf_maxMjj_j2_phi = outputAK5_MjjSorted.at(1).Momentum_.Phi();
     NewElectronTree->vbf_maxMjj_j2_m   = outputAK5_MjjSorted.at(1).Momentum_.M();
   
     NewElectronTree->vbf_maxMjj_jj_deta = fabs(outputAK5_MjjSorted.at(0).Momentum_.Eta() - outputAK5_MjjSorted.at(1).Momentum_.Eta()) ;
     if (fabs(outputAK5_MjjSorted.at(0).Momentum_.Phi() - outputAK5_MjjSorted.at(1).Momentum_.Phi())<3.14)
      NewElectronTree->vbf_maxMjj_jj_dphi = fabs(outputAK5_MjjSorted.at(0).Momentum_.Phi() - outputAK5_MjjSorted.at(1).Momentum_.Phi()) ;
     else
      NewElectronTree->vbf_maxMjj_jj_dphi = 6.28 - fabs(outputAK5_MjjSorted.at(0).Momentum_.Phi() - outputAK5_MjjSorted.at(1).Momentum_.Phi()) ;

     if(outputAK5_MjjSorted.at(0).NameCollection_ == "JetPFCor" && outputAK5_MjjSorted.at(1).NameCollection_ == "JetPFCor"){
      
      NewElectronTree->vbf_maxMjj_jj_type = 1 ;

      int nexcj = 0 , nexfj = 0; 
      std::vector<JetAK5>::const_iterator itVec = outputAK5_MjjSorted.begin();
      for( ; itVec != outputAK5_MjjSorted.end() ; itVec++){

	if(itVec->NameCollection_ == "JetPFCor") nexcj ++ ;
   	if(itVec->NameCollection_ == "JetPFCorVBFTag") nexfj ++ ;
   
      }

      NewElectronTree->vbf_maxMjj_n_excj = nexcj ;
      NewElectronTree->vbf_maxMjj_n_exfj = nexfj ;
     }

     else if(outputAK5_MjjSorted.at(0).NameCollection_ == "JetPFCor" && outputAK5_MjjSorted.at(1).NameCollection_ == "JetPFCorVBFTag"){
      
      NewElectronTree->vbf_maxMjj_jj_type = 2 ;

      int nexcj = 0 , nexfj = 0;
      std::vector<JetAK5>::const_iterator itVec = outputAK5_MjjSorted.begin();
      for( ; itVec != outputAK5_MjjSorted.end() ; itVec++){

	if(itVec->NameCollection_ == "JetPFCor") nexcj ++ ;
   	if(itVec->NameCollection_ == "JetPFCorVBFTag") nexfj ++ ;
   
      }

      NewElectronTree->vbf_maxMjj_n_excj = nexcj ;
      NewElectronTree->vbf_maxMjj_n_exfj = nexfj ;
     }

     else if(outputAK5_MjjSorted.at(0).NameCollection_ == "JetPFCorVBFTag" && outputAK5_MjjSorted.at(1).NameCollection_ == "JetPFCor"){
      
      NewElectronTree->vbf_maxMjj_jj_type = 3 ;

      int nexcj = 0 , nexfj = 0;
      std::vector<JetAK5>::const_iterator itVec = outputAK5_MjjSorted.begin();
      for( ; itVec != outputAK5_MjjSorted.end() ; itVec++){

	if(itVec->NameCollection_ == "JetPFCor") nexcj ++ ;
   	if(itVec->NameCollection_ == "JetPFCorVBFTag") nexfj ++ ;
   
      }

      NewElectronTree->vbf_maxMjj_n_excj = nexcj ;
      NewElectronTree->vbf_maxMjj_n_exfj = nexfj ;
     }

     else if(outputAK5_MjjSorted.at(0).NameCollection_ == "JetPFCorVBFTag" && outputAK5_MjjSorted.at(1).NameCollection_ == "JetPFCorVBFTag"){
      
      NewElectronTree->vbf_maxMjj_jj_type = 4 ;

      int nexcj = 0 , nexfj = 0;
      std::vector<JetAK5>::const_iterator itVec = outputAK5_MjjSorted.begin();
      for( ; itVec != outputAK5_MjjSorted.end() ; itVec++){

	if(itVec->NameCollection_ == "JetPFCor") nexcj ++ ;
   	if(itVec->NameCollection_ == "JetPFCorVBFTag") nexfj ++ ;
   
      }

      NewElectronTree->vbf_maxMjj_n_excj = nexcj ;
      NewElectronTree->vbf_maxMjj_n_exfj = nexfj ;
     }
     else{ std::cerr<<" Something Wrong in MaxMjj Jet Categorization "<<std::endl; continue ;}

     if(std::string(SelectionEvents->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEvents->GetXaxis()->SetBinLabel(nStep,"Jet MaxMjj Category");
     if(std::string(SelectionEfficiency->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEfficiency->GetXaxis()->SetBinLabel(nStep,"Jet MaxMjj Category");
     nstepEvents[nStep-1]++;
     nStep = 12 ;

     if(!numberJetBinGen.empty() && !outputGenAK5_MjjSorted.empty()){
      if(numberJetBinGen.at(0)>=2 && outputGenAK5_MjjSorted.size()>=2){

      TLorentzVector Total4VMaxMjjGen = outputGenAK5_MjjSorted.at(0).Momentum_ + outputGenAK5_MjjSorted.at(1).Momentum_ ;
      NewElectronTree->vbf_maxMjj_jj_e_gen   = Total4VMaxMjjGen.E(); 
      NewElectronTree->vbf_maxMjj_jj_pt_gen  = Total4VMaxMjjGen.Pt(); 
      NewElectronTree->vbf_maxMjj_jj_eta_gen = Total4VMaxMjjGen.Eta(); 
      NewElectronTree->vbf_maxMjj_jj_phi_gen = Total4VMaxMjjGen.Phi(); 
      NewElectronTree->vbf_maxMjj_jj_m_gen   = Total4VMaxMjjGen.M(); 
      NewElectronTree->vbf_maxMjj_j1_e_gen   = outputGenAK5_MjjSorted.at(0).Momentum_.E();
      NewElectronTree->vbf_maxMjj_j1_pt_gen  = outputGenAK5_MjjSorted.at(0).Momentum_.Pt();
      NewElectronTree->vbf_maxMjj_j1_eta_gen = outputGenAK5_MjjSorted.at(0).Momentum_.Eta();
      NewElectronTree->vbf_maxMjj_j1_phi_gen = outputGenAK5_MjjSorted.at(0).Momentum_.Phi();
      NewElectronTree->vbf_maxMjj_j1_m_gen   = outputGenAK5_MjjSorted.at(0).Momentum_.M(); 
      NewElectronTree->vbf_maxMjj_j2_e_gen   = outputGenAK5_MjjSorted.at(1).Momentum_.E();
      NewElectronTree->vbf_maxMjj_j2_pt_gen  = outputGenAK5_MjjSorted.at(1).Momentum_.Pt();
      NewElectronTree->vbf_maxMjj_j2_eta_gen = outputGenAK5_MjjSorted.at(1).Momentum_.Eta();
      NewElectronTree->vbf_maxMjj_j2_phi_gen = outputGenAK5_MjjSorted.at(1).Momentum_.Phi();
      NewElectronTree->vbf_maxMjj_j2_m_gen   = outputGenAK5_MjjSorted.at(1).Momentum_.M();

      NewElectronTree->vbf_maxMjj_jj_deta_gen = fabs(outputGenAK5_MjjSorted.at(0).Momentum_.Eta() - outputGenAK5_MjjSorted.at(1).Momentum_.Eta()) ;
      if (fabs(outputGenAK5_MjjSorted.at(0).Momentum_.Phi() - outputGenAK5_MjjSorted.at(1).Momentum_.Phi()) < TMath::Pi())
       NewElectronTree->vbf_maxMjj_jj_dphi_gen = fabs(outputGenAK5_MjjSorted.at(0).Momentum_.Phi() - outputGenAK5_MjjSorted.at(1).Momentum_.Phi()) ;
      else 
       NewElectronTree->vbf_maxMjj_jj_dphi_gen = 2*TMath::Pi() - fabs(outputGenAK5_MjjSorted.at(0).Momentum_.Phi() - outputGenAK5_MjjSorted.at(1).Momentum_.Phi()) ;

      NewElectronTree->vbf_maxMjj_j1_bDiscriminatorSSVHE_gen = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorSSVHE")[outputGenAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_bDiscriminatorTCHE_gen  = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorTCHE")[outputGenAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_bDiscriminatorCSV_gen   = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorCSV")[outputGenAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_bDiscriminatorSSVHP_gen = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorSSVHP")[outputGenAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_bDiscriminatorTCHP_gen  = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorTCHP")[outputGenAK5_MjjSorted.at(0).position_] ;

      NewElectronTree->vbf_maxMjj_j2_bDiscriminatorSSVHE_gen = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorSSVHE")[outputGenAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_bDiscriminatorTCHE_gen  = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorTCHE")[outputGenAK5_MjjSorted.at(1).position_] ; 
      NewElectronTree->vbf_maxMjj_j2_bDiscriminatorCSV_gen   = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorCSV")[outputGenAK5_MjjSorted.at(1).position_]  ;
      NewElectronTree->vbf_maxMjj_j2_bDiscriminatorSSVHP_gen = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorSSVHP")[outputGenAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_bDiscriminatorTCHP_gen  = ElectronTree->fReader->getFloat("JetGen_bDiscriminatorTCHP")[outputGenAK5_MjjSorted.at(1).position_] ;
      }
     }

     if(outputAK5_MjjSorted.at(0).NameCollection_ == "JetPFCor") {
    
      NewElectronTree->vbf_maxMjj_j1_QGLikelihood = ElectronTree->fReader->getFloat("JetPFCor_QGLikelihood")[outputAK5_MjjSorted.at(0).position_] ;

      NewElectronTree->vbf_maxMjj_j1_isPileUpLoose  = ElectronTree->fReader->getBool("JetPFCor_isPileUpJetLoose")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_isPileUpMedium = ElectronTree->fReader->getBool("JetPFCor_isPileUpJetMedium")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_isPileUpTight  = ElectronTree->fReader->getBool("JetPFCor_isPileUpJetTight")[outputAK5_MjjSorted.at(0).position_] ;


      NewElectronTree->vbf_maxMjj_j1_bDiscriminatorSSVHE = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorSSVHE")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_bDiscriminatorTCHE  = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorTCHE")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_bDiscriminatorCSV   = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorCSV")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_bDiscriminatorSSVHP = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorSSVHP")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_bDiscriminatorTCHP  = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorTCHP")[outputAK5_MjjSorted.at(0).position_] ;

      NewElectronTree->vbf_maxMjj_j1_ChargedHadronEnergy      = ElectronTree->fReader->getFloat("JetPFCor_ChargedHadronEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_ChargedHadronEnergyFrac  = ElectronTree->fReader->getFloat("JetPFCor_ChargedHadronEnergyFrac")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_NeutralHadronEnergy      = ElectronTree->fReader->getFloat("JetPFCor_NeutralHadronEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_NeutralHadronEnergyFrac  = ElectronTree->fReader->getFloat("JetPFCor_NeutralHadronEnergyFrac")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_ChargedEmEnergy          = ElectronTree->fReader->getFloat("JetPFCor_ChargedEmEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_ChargedEmEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCor_ChargedEmEnergyFrac")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_ChargedMuEnergy          = ElectronTree->fReader->getFloat("JetPFCor_ChargedMuEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_ChargedMuEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCor_ChargedMuEnergyFrac")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_NeutralEmEnergy          = ElectronTree->fReader->getFloat("JetPFCor_NeutralEmEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_NeutralEmEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCor_NeutralEmEnergyFrac")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_PhotonEnergy             = ElectronTree->fReader->getFloat("JetPFCor_PhotonEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_PhotonEnergyFraction     = ElectronTree->fReader->getFloat("JetPFCor_PhotonEnergyFraction")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_ElectronEnergy           = ElectronTree->fReader->getFloat("JetPFCor_ElectronEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_ElectronEnergyFraction   = ElectronTree->fReader->getFloat("JetPFCor_ElectronEnergyFraction")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_HFHadronEnergy           = ElectronTree->fReader->getFloat("JetPFCor_HFHadronEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_HFHadronEnergyFraction   = ElectronTree->fReader->getFloat("JetPFCor_HFHadronEnergyFraction")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_HFEMEnergy               = ElectronTree->fReader->getFloat("JetPFCor_HFEMEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_HFEMEnergyFraction       = ElectronTree->fReader->getFloat("JetPFCor_HFEMEnergyFraction")[outputAK5_MjjSorted.at(0).position_] ;

      NewElectronTree->vbf_maxMjj_j1_ChargedMultiplicity        = ElectronTree->fReader->getFloat("JetPFCor_ChargedMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_NeutralMultiplicity        = ElectronTree->fReader->getFloat("JetPFCor_NeutralMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_MuonMultiplicity           = ElectronTree->fReader->getFloat("JetPFCor_MuonMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_ChargedHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCor_ChargedHadronMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_NeutralHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCor_NeutralHadronMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_PhotonMultiplicity         = ElectronTree->fReader->getFloat("JetPFCor_PhotonMultiplicity")[outputAK5_MjjSorted.at(0).position_] ; 
      NewElectronTree->vbf_maxMjj_j1_ElectronMultiplicity       = ElectronTree->fReader->getFloat("JetPFCor_ElectronMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_HFHadronMultiplicity       = ElectronTree->fReader->getFloat("JetPFCor_HFHadronMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
     }
     else if(outputAK5_MjjSorted.at(0).NameCollection_ == "JetPFCorVBFTag" ){

      NewElectronTree->vbf_maxMjj_j1_QGLikelihood = ElectronTree->fReader->getFloat("JetPFCorVBFTag_QGLikelihood")[outputAK5_MjjSorted.at(0).position_] ;

      NewElectronTree->vbf_maxMjj_j1_isPileUpLoose  = ElectronTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetLoose")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_isPileUpMedium = ElectronTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetMedium")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_isPileUpTight  = ElectronTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetTight")[outputAK5_MjjSorted.at(0).position_] ;


      NewElectronTree->vbf_maxMjj_j1_bDiscriminatorSSVHE = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorSSVHE")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_bDiscriminatorTCHE  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorTCHE")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_bDiscriminatorCSV   = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorCSV")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_bDiscriminatorSSVHP = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorSSVHP")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_bDiscriminatorTCHP  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorTCHP")[outputAK5_MjjSorted.at(0).position_] ;

      NewElectronTree->vbf_maxMjj_j1_ChargedHadronEnergy      = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_ChargedHadronEnergyFrac  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronEnergyFrac")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_NeutralHadronEnergy      = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_NeutralHadronEnergyFrac  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronEnergyFrac")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_ChargedEmEnergy          = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedEmEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_ChargedEmEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedEmEnergyFrac")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_ChargedMuEnergy          = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedMuEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_ChargedMuEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedMuEnergyFrac")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_NeutralEmEnergy          = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralEmEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_NeutralEmEnergyFrac  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralEmEnergyFrac")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_PhotonEnergy         = ElectronTree->fReader->getFloat("JetPFCorVBFTag_PhotonEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_PhotonEnergyFraction = ElectronTree->fReader->getFloat("JetPFCorVBFTag_PhotonEnergyFraction")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_ElectronEnergy       = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ElectronEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_ElectronEnergyFraction   = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ElectronEnergyFraction")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_HFHadronEnergy           = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFHadronEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_HFHadronEnergyFraction   = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFHadronEnergyFraction")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_HFEMEnergy               = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFEMEnergy")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_HFEMEnergyFraction       = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFEMEnergyFraction")[outputAK5_MjjSorted.at(0).position_] ;

      NewElectronTree->vbf_maxMjj_j1_ChargedMultiplicity        = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_NeutralMultiplicity        = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_MuonMultiplicity           = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ElectronMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_ChargedHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_NeutralHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_PhotonMultiplicity         = ElectronTree->fReader->getFloat("JetPFCorVBFTag_PhotonMultiplicity")[outputAK5_MjjSorted.at(0).position_] ; 
      NewElectronTree->vbf_maxMjj_j1_ElectronMultiplicity       = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ElectronMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_HFHadronMultiplicity       = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFHadronMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;  

     }
     else { std::cerr<<" problem with High Deta Jet Name Collection "<<std::endl; continue ; }

     if(outputAK5_MjjSorted.at(1).NameCollection_ == "JetPFCor") {
    
      NewElectronTree->vbf_maxMjj_j2_QGLikelihood = ElectronTree->fReader->getFloat("JetPFCor_QGLikelihood")[outputAK5_MjjSorted.at(1).position_] ;

      NewElectronTree->vbf_maxMjj_j2_isPileUpLoose = ElectronTree->fReader->getBool("JetPFCor_isPileUpJetLoose")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_isPileUpMedium = ElectronTree->fReader->getBool("JetPFCor_isPileUpJetMedium")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_isPileUpTight = ElectronTree->fReader->getBool("JetPFCor_isPileUpJetTight")[outputAK5_MjjSorted.at(1).position_] ;


      NewElectronTree->vbf_maxMjj_j2_bDiscriminatorSSVHE = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorSSVHE")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_bDiscriminatorTCHE  = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorTCHE")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_bDiscriminatorCSV   = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorCSV")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_bDiscriminatorSSVHP = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorSSVHP")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_bDiscriminatorTCHP  = ElectronTree->fReader->getFloat("JetPFCor_bDiscriminatorTCHP")[outputAK5_MjjSorted.at(1).position_] ;

      NewElectronTree->vbf_maxMjj_j2_ChargedHadronEnergy      = ElectronTree->fReader->getFloat("JetPFCor_ChargedHadronEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_ChargedHadronEnergyFrac  = ElectronTree->fReader->getFloat("JetPFCor_ChargedHadronEnergyFrac")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_NeutralHadronEnergy      = ElectronTree->fReader->getFloat("JetPFCor_NeutralHadronEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_NeutralHadronEnergyFrac  = ElectronTree->fReader->getFloat("JetPFCor_NeutralHadronEnergyFrac")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_ChargedEmEnergy          = ElectronTree->fReader->getFloat("JetPFCor_ChargedEmEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_ChargedEmEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCor_ChargedEmEnergyFrac")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_ChargedMuEnergy          = ElectronTree->fReader->getFloat("JetPFCor_ChargedMuEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_ChargedMuEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCor_ChargedMuEnergyFrac")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_NeutralEmEnergy          = ElectronTree->fReader->getFloat("JetPFCor_NeutralEmEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_NeutralEmEnergyFrac  = ElectronTree->fReader->getFloat("JetPFCor_NeutralEmEnergyFrac")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_PhotonEnergy         = ElectronTree->fReader->getFloat("JetPFCor_PhotonEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_PhotonEnergyFraction = ElectronTree->fReader->getFloat("JetPFCor_PhotonEnergyFraction")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_ElectronEnergy       = ElectronTree->fReader->getFloat("JetPFCor_ElectronEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_ElectronEnergyFraction   = ElectronTree->fReader->getFloat("JetPFCor_ElectronEnergyFraction")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_HFHadronEnergy           = ElectronTree->fReader->getFloat("JetPFCor_HFHadronEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_HFHadronEnergyFraction   = ElectronTree->fReader->getFloat("JetPFCor_HFHadronEnergyFraction")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_HFEMEnergy               = ElectronTree->fReader->getFloat("JetPFCor_HFEMEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_HFEMEnergyFraction       = ElectronTree->fReader->getFloat("JetPFCor_HFEMEnergyFraction")[outputAK5_MjjSorted.at(1).position_] ;

      NewElectronTree->vbf_maxMjj_j2_ChargedMultiplicity        = ElectronTree->fReader->getFloat("JetPFCor_ChargedMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_NeutralMultiplicity        = ElectronTree->fReader->getFloat("JetPFCor_NeutralMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_ElectronMultiplicity           = ElectronTree->fReader->getFloat("JetPFCor_ElectronMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_ChargedHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCor_ChargedHadronMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_NeutralHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCor_NeutralHadronMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_PhotonMultiplicity         = ElectronTree->fReader->getFloat("JetPFCor_PhotonMultiplicity")[outputAK5_MjjSorted.at(1).position_] ; 
      NewElectronTree->vbf_maxMjj_j2_ElectronMultiplicity       = ElectronTree->fReader->getFloat("JetPFCor_ElectronMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_HFHadronMultiplicity       = ElectronTree->fReader->getFloat("JetPFCor_HFHadronMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
     }
     else if(outputAK5_MjjSorted.at(1).NameCollection_ == "JetPFCorVBFTag" ){

      NewElectronTree->vbf_maxMjj_j2_QGLikelihood = ElectronTree->fReader->getFloat("JetPFCorVBFTag_QGLikelihood")[outputAK5_MjjSorted.at(1).position_] ;

      NewElectronTree->vbf_maxMjj_j2_isPileUpLoose  = ElectronTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetLoose")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_isPileUpMedium = ElectronTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetMedium")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_isPileUpTight  = ElectronTree->fReader->getBool("JetPFCorVBFTag_isPileUpJetTight")[outputAK5_MjjSorted.at(1).position_] ;


      NewElectronTree->vbf_maxMjj_j2_bDiscriminatorSSVHE = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorSSVHE")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_bDiscriminatorTCHE  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorTCHE")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_bDiscriminatorCSV   = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorCSV")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_bDiscriminatorSSVHP = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorSSVHP")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_bDiscriminatorTCHP  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_bDiscriminatorTCHP")[outputAK5_MjjSorted.at(1).position_] ;

      NewElectronTree->vbf_maxMjj_j2_ChargedHadronEnergy      = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_ChargedHadronEnergyFrac  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronEnergyFrac")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_NeutralHadronEnergy      = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_NeutralHadronEnergyFrac  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronEnergyFrac")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_ChargedEmEnergy          = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedEmEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_ChargedEmEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedEmEnergyFrac")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_ChargedMuEnergy          = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedMuEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_ChargedMuEnergyFrac      = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedMuEnergyFrac")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_NeutralEmEnergy          = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralEmEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_NeutralEmEnergyFrac  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralEmEnergyFrac")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_PhotonEnergy         = ElectronTree->fReader->getFloat("JetPFCorVBFTag_PhotonEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_PhotonEnergyFraction = ElectronTree->fReader->getFloat("JetPFCorVBFTag_PhotonEnergyFraction")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_ElectronEnergy       = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ElectronEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_ElectronEnergyFraction   = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ElectronEnergyFraction")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_HFHadronEnergy           = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFHadronEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_HFHadronEnergyFraction   = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFHadronEnergyFraction")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_HFEMEnergy               = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFEMEnergy")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_HFEMEnergyFraction       = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFEMEnergyFraction")[outputAK5_MjjSorted.at(1).position_] ;

      NewElectronTree->vbf_maxMjj_j2_ChargedMultiplicity  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_NeutralMultiplicity  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_MuonMultiplicity     = ElectronTree->fReader->getFloat("JetPFCorVBFTag_MuonMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_ChargedHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_NeutralHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_PhotonMultiplicity         = ElectronTree->fReader->getFloat("JetPFCorVBFTag_PhotonMultiplicity")[outputAK5_MjjSorted.at(1).position_] ; 
      NewElectronTree->vbf_maxMjj_j2_ElectronMultiplicity       = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ElectronMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_HFHadronMultiplicity       = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFHadronMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;  
    
     }
     else { std::cerr<<" problem with High Deta Jet Name Collection "<<std::endl; continue ; }
          
    }
    NewElectronTree->fTree->Fill(); // Fill the events 
    
   } // End of Loop on the event

   // Save Results in the output
   NewElectronTree->fTree->Write(TreeName.c_str());

   int ibin = 0 ;
   nStep = 0 ;
   for( ; nstepEvents[ibin]!=0 ; ibin ++){ nStep ++ ; SelectionEvents->SetBinContent(ibin+1,nstepEvents[ibin]); }

   for(ibin =1 ; ibin<SelectionEvents->GetNbinsX() ; ibin ++) {
    if(SelectionEvents->GetBinContent(ibin+1)!=0) SelectionEfficiency->SetBinContent(ibin+1,SelectionEvents->GetBinContent(ibin+1)/SelectionEvents->GetBinContent(ibin)); 
   }

   SelectionEvents->GetXaxis()->SetRangeUser(0,nStep-1);
   SelectionEvents->Write();

   SelectionEfficiency->GetXaxis()->SetRangeUser(1,nStep-1);
   SelectionEfficiency->Write();

   std::cout << " Finish :: " << outputFile->GetName() << "    "<<  ElectronTree->fTree->GetEntries ()  << std::endl;
   outputFile->Close();
  } // End of Electron Analysis
  
 return 0 ;

}
  
