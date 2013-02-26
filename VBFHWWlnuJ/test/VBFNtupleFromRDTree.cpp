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

  double JetPtWboostedMin        = gConfigParser -> readDoubleOption("Input::JetPtWboostedMin");
  double JetPtCutMin             = gConfigParser -> readDoubleOption("Input::JetPtCutMin");
  double JetEtaCutMax            = gConfigParser -> readDoubleOption("Input::JetEtaCutMax");
  double CleaningTreshold        = gConfigParser -> readDoubleOption("Input::CleaningTreshold");

  std::cout<<"                   "<<std::endl;
  std::cout<<" Input Directory   "<<InputDirectory<<std::endl;
  std::cout<<" Input Root File   "<<InputRootFile<<std::endl;
  std::cout<<" Input TreeName    "<<TreeName<<std::endl;
  std::cout<<" Input Lepton Type "<<LeptonType<<std::endl;
  std::cout<<"                   "<<std::endl;

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

   TFile *outputFile = new TFile((OutputRootDirectory+"/Mu"+OutputRootFile).c_str(),"RECREATE");
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

    // Control the Jet Multiplicity for VBF regime
    if( (MuonTree->fReader->getInt("numPFCorJets")[0] +MuonTree->fReader->getInt("numPFCorVBFTagJets")[0]) < NumJetMin ) continue ;

    if(std::string(SelectionEvents->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEvents->GetXaxis()->SetBinLabel(nStep,"Jet Number"); 
    if(std::string(SelectionEfficiency->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEfficiency->GetXaxis()->SetBinLabel(nStep,"Jet Number"); 
    nstepEvents[nStep-1]++;
    nStep = 4;

    // Join forward and central PFJetCor Collection    
    std::vector<TLorentzVector> GroomedJet_CA8_Collection ;

    std::vector<JetAK5> JetPFCor_AK5_Collection ;
    std::vector<JetAK5> CleanedJetPFCor_AK5_Collection ;
    std::vector<JetAK5> HadronicW_AK5_Collection ;
    
    for(int iJet = 0 ; iJet < JetCollectionDimension ; iJet++) {
     TLorentzVector JetTemp ; std::string nameCollection ; 
     JetTemp.SetPtEtaPhiE(MuonTree->fReader->getFloat("GroomedJet_CA8_pt")[iJet],MuonTree->fReader->getFloat("GroomedJet_CA8_eta")[iJet], 
                          MuonTree->fReader->getFloat("GroomedJet_CA8_phi")[iJet],MuonTree->fReader->getFloat("GroomedJet_CA8_e")[iJet]);

     // Selection on CA8 Jets
     if(fabs(JetTemp.Eta())<JetEtaCutMax && JetTemp.Pt()>JetPtCutMin)
       GroomedJet_CA8_Collection.push_back(JetTemp);

     JetTemp.SetPtEtaPhiE(MuonTree->fReader->getFloat("JetPFCor_Pt")[iJet],MuonTree->fReader->getFloat("JetPFCor_Eta")[iJet],
 			  MuonTree->fReader->getFloat("JetPFCor_Phi")[iJet],MuonTree->fReader->getFloat("JetPFCor_E")[iJet]);
     // Selection on PF Cor Central Jets --> AK5
     if(fabs(JetTemp.Eta())<JetEtaCutMax && JetTemp.Pt()>JetPtCutMin){
       JetAK5 tempJetAK5 (iJet,"JetPFCor",JetTemp);
       JetPFCor_AK5_Collection.push_back(tempJetAK5);
     }
     
    }
   
    for(int iJet = 0 ; iJet < JetCollectionDimension ; iJet++) {
     TLorentzVector JetTemp ; 
     JetTemp.SetPtEtaPhiE(MuonTree->fReader->getFloat("JetPFCorVBFTag_Pt")[iJet],MuonTree->fReader->getFloat("JetPFCorVBFTag_Eta")[iJet], 
                          MuonTree->fReader->getFloat("JetPFCorVBFTag_Phi")[iJet],MuonTree->fReader->getFloat("JetPFCorVBFTag_E")[iJet]);

     // Selection on PF Cor Forward Jets --> AK5
     if(fabs(JetTemp.Eta())<JetEtaCutMax && JetTemp.Pt()>JetPtCutMin){
       JetAK5 tempJetAK5 (iJet,"JetPFCorVBFTag",JetTemp);
       JetPFCor_AK5_Collection.push_back(tempJetAK5);
     }
    }

    // Another Skim on Jet Counting
    if(GroomedJet_CA8_Collection.size() == 0 || JetPFCor_AK5_Collection.size() < 3) continue ;
   
    if(std::string(SelectionEvents->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEvents->GetXaxis()->SetBinLabel(nStep,"Jet Number CA8 - AK5");
    if(std::string(SelectionEfficiency->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEfficiency->GetXaxis()->SetBinLabel(nStep,"Jet Number CA8 - AK5"); 
    nstepEvents[nStep-1]++;
    nStep = 5;

    // Clean AK5 Jet Collection choosing the hard CA8 as W Hadronic

    for(size_t iJet = 0; iJet < JetPFCor_AK5_Collection.size() ; iJet ++){

      if(deltaR(JetPFCor_AK5_Collection.at(iJet).Momentum_.Phi(),GroomedJet_CA8_Collection.at(0).Phi(),
                JetPFCor_AK5_Collection.at(iJet).Momentum_.Eta(),GroomedJet_CA8_Collection.at(0).Eta()) < CleaningTreshold ){
	HadronicW_AK5_Collection.push_back(JetPFCor_AK5_Collection.at(iJet)); continue ;}

      CleanedJetPFCor_AK5_Collection.push_back(JetPFCor_AK5_Collection.at(iJet));

    }

    if(CleanedJetPFCor_AK5_Collection.size() < 2) continue ; 

    if(std::string(SelectionEvents->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEvents->GetXaxis()->SetBinLabel(nStep,"Matching CA8 - AK5");
    if(std::string(SelectionEfficiency->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEfficiency->GetXaxis()->SetBinLabel(nStep,"Matching CA8 - AK5");
    nstepEvents[nStep-1]++;
    nStep = 6 ;
    // vbf Tag Jet Selection
    
    std::vector<JetAK5> outputAK5_PtSorted;
    std::vector<JetAK5> outputAK5_DEtaSorted;
    std::vector<JetAK5> outputAK5_MjjSorted;

    // Sorting of AK5 Cleaned Collection in Pt

    std::sort(CleanedJetPFCor_AK5_Collection.begin(),CleanedJetPFCor_AK5_Collection.end(),TLVP_PtSort());
    outputAK5_PtSorted.push_back(CleanedJetPFCor_AK5_Collection.at(0));
    outputAK5_PtSorted.push_back(CleanedJetPFCor_AK5_Collection.at(1));
    if(outputAK5_PtSorted.size() < 2) continue ;
    
    // Sorting of AK5 Cleaned Collection in DeltaEta

    std::sort(CleanedJetPFCor_AK5_Collection.begin(),CleanedJetPFCor_AK5_Collection.end(),TLVP_EtaSort());
    outputAK5_DEtaSorted.push_back(CleanedJetPFCor_AK5_Collection.front());
    outputAK5_DEtaSorted.push_back(CleanedJetPFCor_AK5_Collection.back());
    if(outputAK5_DEtaSorted.size() < 2) continue ;

    // Sorting of AK5 Cleaned Collection in Mjj
    float maxMjj = 0. ;
    int iJ1 = 0 ;
    int iJ2 = 0 ;

    for (size_t iJet = 0 ; iJet < CleanedJetPFCor_AK5_Collection.size()-1 ; ++iJet){
      for (size_t jJet = iJet + 1 ; jJet < CleanedJetPFCor_AK5_Collection.size() ; ++jJet){

	TLorentzVector SumMomentum = CleanedJetPFCor_AK5_Collection.at(iJet).Momentum_ + CleanedJetPFCor_AK5_Collection.at(jJet).Momentum_ ;
        float Mjj = SumMomentum.M();
        if(Mjj > maxMjj){
         
          iJ1 = iJet ;
          iJ2 = jJet ;
          maxMjj = Mjj ;
        }
      }
    }

    outputAK5_MjjSorted.push_back (CleanedJetPFCor_AK5_Collection.at (iJ1)) ;
    outputAK5_MjjSorted.push_back (CleanedJetPFCor_AK5_Collection.at (iJ2)) ;
    if(outputAK5_MjjSorted.size() < 2) continue ;

    // Calculate Neutrino Pz
    TLorentzVector W_mu, W_Met, W_neutrino; 
   
    W_mu.SetPxPyPzE(MuonTree->fReader->getFloat("W_muon_px")[0],MuonTree->fReader->getFloat("W_muon_py")[0],
                    MuonTree->fReader->getFloat("W_muon_pz")[0],MuonTree->fReader->getFloat("W_muon_e")[0]);
    W_Met.SetPxPyPzE(MuonTree->fReader->getFloat("event_met_pfmet")[0] * TMath::Cos(MuonTree->fReader->getFloat("event_met_pfmetPhi")[0]),
                     MuonTree->fReader->getFloat("event_met_pfmet")[0] * TMath::Sin(MuonTree->fReader->getFloat("event_met_pfmetPhi")[0]),0.,
                     fabs(MuonTree->fReader->getFloat("event_met_pfmet")[0]));

    if(W_mu.Pt()<=0 || W_Met.Pt() <= 0 ){ std::cerr<<" Negative Lepton - Neutrino Pt "<<std::endl; continue ; }

    if(std::string(SelectionEvents->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEvents->GetXaxis()->SetBinLabel(nStep,"Leptonic 4V");
    if(std::string(SelectionEfficiency->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEfficiency->GetXaxis()->SetBinLabel(nStep,"Leptonic 4V");
    nstepEvents[nStep-1]++;
    nStep = 7;

    METzCalculator<TLorentzVector> NeutrinoPz;
    NeutrinoPz.SetMET(W_Met);
    NeutrinoPz.SetLepton(W_mu);
    NeutrinoPz.SetLeptonType("muon");
    double pz = NeutrinoPz.Calculate(); // Default one
    W_neutrino.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz,sqrt(W_Met.Pt()*W_Met.Pt()+pz*pz));
    if (NeutrinoPz.IsComplex()) {// if this is a complix, change MET                                                                                                                 
           
     double nu_pt1 = NeutrinoPz.getPtneutrino(1);
     double nu_pt2 = NeutrinoPz.getPtneutrino(2);
   
     TLorentzVector W_neutrino_1;
     W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(MuonTree->fReader->getFloat("event_met_pfmetPhi")[0]), 
                             nu_pt1 * TMath::Sin(MuonTree->fReader->getFloat("event_met_pfmetPhi")[0]), pz, sqrt(nu_pt1*nu_pt1 + pz*pz) );
     TLorentzVector W_neutrino_2;
     W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(MuonTree->fReader->getFloat("event_met_pfmetPhi")[0]),
                             nu_pt2 * TMath::Sin(MuonTree->fReader->getFloat("event_met_pfmetPhi")[0]), pz, sqrt(nu_pt2*nu_pt2 + pz*pz) );

     if ( fabs((W_mu+W_neutrino_1).M()-80.4) < fabs((W_mu+W_neutrino_2).M()-80.4) )  W_neutrino = W_neutrino_1;
     else W_neutrino = W_neutrino_2;

    }
    
    TLorentzVector W_subjet1, W_subjet2 ;  // take the two subjet of the hardest CA8 after pruning
   
    W_subjet1.SetPxPyPzE(MuonTree->fReader->getFloat("GroomedJet_CA8_prsubjet1_px")[0],MuonTree->fReader->getFloat("GroomedJet_CA8_prsubjet1_py")[0],
         		 MuonTree->fReader->getFloat("GroomedJet_CA8_prsubjet1_pz")[0],MuonTree->fReader->getFloat("GroomedJet_CA8_prsubjet1_e")[0] );
    W_subjet2.SetPxPyPzE(MuonTree->fReader->getFloat("GroomedJet_CA8_prsubjet2_px")[0],MuonTree->fReader->getFloat("GroomedJet_CA8_prsubjet2_py")[0],
  		  	 MuonTree->fReader->getFloat("GroomedJet_CA8_prsubjet2_pz")[0],MuonTree->fReader->getFloat("GroomedJet_CA8_prsubjet2_e")[0] );

    if(W_subjet1.Pt() <= 0 || W_subjet2.Pt() <= 0){ std::cerr<<" Problem with subjets "<<std::endl; continue ; }

    if(std::string(SelectionEvents->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEvents->GetXaxis()->SetBinLabel(nStep,"Subjets 4V");
    if(std::string(SelectionEfficiency->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEfficiency->GetXaxis()->SetBinLabel(nStep,"Subjets 4V");
    nstepEvents[nStep-1]++;
    nStep = 8;


    TLorentzVector W_GroomedJet_CA8_pr; 
    W_GroomedJet_CA8_pr.SetPtEtaPhiE(MuonTree->fReader->getFloat("GroomedJet_CA8_pt_pr")[0], MuonTree->fReader->getFloat("GroomedJet_CA8_eta_pr")[0],
                                     MuonTree->fReader->getFloat("GroomedJet_CA8_phi_pr")[0], MuonTree->fReader->getFloat("GroomedJet_CA8_e_pr")[0]);

    if(W_GroomedJet_CA8_pr.Pt() <=0){ std::cerr<<" Problem with pruned CA8 "<<std::endl; continue ;}

    if(std::string(SelectionEvents->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEvents->GetXaxis()->SetBinLabel(nStep,"Pruned CA8 4V");
    if(std::string(SelectionEfficiency->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEfficiency->GetXaxis()->SetBinLabel(nStep,"Pruned CA8 4V");
    nstepEvents[nStep-1]++;
    nStep = 9;

   
    // Kinematic Fit                                                                                                                                                                 
    TLorentzVector fit_muon(0,0,0,0), fit_neutrino(0,0,0,0), fit_W_subjet1(0,0,0,0), fit_W_subjet2(0,0,0,0) ;

    doKinematicFit(1, W_mu, W_neutrino, W_subjet1, W_subjet2,  fit_muon, fit_neutrino, fit_W_subjet1, fit_W_subjet2, NewMuonTree->fit_chi2, 
                   NewMuonTree->fit_NDF, NewMuonTree->fit_status, LeptonType);
    
    if(fit_muon.Pt() >0 && fit_neutrino.Pt()>0 && fit_W_subjet1.Pt()>0 && fit_W_subjet2.Pt()>0){

    NewMuonTree->fit_mu_px = fit_muon.Px();
    NewMuonTree->fit_mu_py = fit_muon.Py(); 
    NewMuonTree->fit_mu_pz = fit_muon.Pz(); 
    NewMuonTree->fit_mu_e  = fit_muon.E();
    NewMuonTree->fit_nv_px = fit_neutrino.Px(); 
    NewMuonTree->fit_nv_py = fit_neutrino.Py(); 
    NewMuonTree->fit_nv_pz = fit_neutrino.Pz(); 
    NewMuonTree->fit_nv_e  = fit_neutrino.E();

    NewMuonTree->fit_subjet1_px = fit_W_subjet1.Px();  NewMuonTree->fit_subjet2_px = fit_W_subjet2.Px();
    NewMuonTree->fit_subjet1_py = fit_W_subjet1.Py();  NewMuonTree->fit_subjet2_py = fit_W_subjet2.Py();  
    NewMuonTree->fit_subjet1_pz = fit_W_subjet1.Pz();  NewMuonTree->fit_subjet2_pz = fit_W_subjet2.Pz();
    NewMuonTree->fit_subjet1_e  = fit_W_subjet1.E();   NewMuonTree->fit_subjet2_e  = fit_W_subjet2.E();
    NewMuonTree->fit_subjet1_m  = fit_W_subjet1.M();   NewMuonTree->fit_subjet2_m  = fit_W_subjet2.M();

    NewMuonTree->fit_lvj_m   = (W_mu+W_neutrino+fit_W_subjet1+fit_W_subjet2).M();
    NewMuonTree->fit_lv_m    = (W_mu+W_neutrino).M();
    NewMuonTree->fit_j_m     = (fit_W_subjet1+fit_W_subjet2).M();
    NewMuonTree->fit_lvj_pt  = (fit_muon+fit_neutrino+fit_W_subjet1+fit_W_subjet2).M();
    NewMuonTree->fit_lvj_eta = (fit_muon+fit_neutrino+fit_W_subjet1+fit_W_subjet2).Eta();
    NewMuonTree->fit_lvj_phi = (fit_muon+fit_neutrino+fit_W_subjet1+fit_W_subjet2).Phi();
    NewMuonTree->fit_lvj_e   = (fit_muon+fit_neutrino+fit_W_subjet1+fit_W_subjet2).E();

    NewMuonTree->boosted_lvj_m   = (W_mu+W_neutrino+W_subjet1+W_subjet2).M();
    NewMuonTree->boosted_lv_m    = (W_mu+W_neutrino).M();
    NewMuonTree->boosted_j_m     = (W_subjet1+W_subjet2).M();
    NewMuonTree->boosted_lvj_pt  = (W_mu+W_neutrino+W_subjet1+W_subjet2).M ();
    NewMuonTree->boosted_lvj_eta = (W_mu+W_neutrino+W_subjet1+W_subjet2).Eta();
    NewMuonTree->boosted_lvj_phi = (W_mu+W_neutrino+W_subjet1+W_subjet2).Phi();
    NewMuonTree->boosted_lvj_e   = (W_mu+W_neutrino+W_subjet1+W_subjet2).E();

    NewMuonTree->boostedW_lvj_m   = (W_mu+W_neutrino+W_GroomedJet_CA8_pr).M();
    NewMuonTree->boostedW_lv_m    = (W_mu+W_neutrino).M();
    NewMuonTree->boostedW_j_m     = (W_GroomedJet_CA8_pr).M();
    NewMuonTree->boostedW_lvj_pt  = (W_mu+W_neutrino+W_GroomedJet_CA8_pr).M ();
    NewMuonTree->boostedW_lvj_eta = (W_mu+W_neutrino+W_GroomedJet_CA8_pr).Eta();
    NewMuonTree->boostedW_lvj_phi = (W_mu+W_neutrino+W_GroomedJet_CA8_pr).Phi();
    NewMuonTree->boostedW_lvj_e   = (W_mu+W_neutrino+W_GroomedJet_CA8_pr).E();
    // Angles for the central Higgs Kinematics
    }

    double costheta1, costheta2, phi, costhetastar, phistar1, phistar2;

    //Use the Subjet in the Boosted W Analyisis                                                                                                                                             
    if (MuonTree->fReader->getFloat("W_muon_charge")[0] < 0) calculateAngles(W_mu, W_neutrino,W_subjet1,W_subjet2,costheta1,costheta2,phi,costhetastar,phistar1,phistar2);
    else calculateAngles(W_neutrino, W_mu, W_subjet1, W_subjet2, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);
   
    NewMuonTree->boosted_wjj_ang_ha   = costheta1;
    NewMuonTree->boosted_wjj_ang_hb   = fabs(costheta2); 
    NewMuonTree->boosted_wjj_ang_hs   = costhetastar;
    NewMuonTree->boosted_wjj_ang_phi  = phi;
    NewMuonTree->boosted_wjj_ang_phia = phistar1;																		    NewMuonTree->boosted_wjj_ang_phib = phistar2;
   
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
    NewMuonTree->vbf_maxpt_jj_deta = fabs(outputAK5_PtSorted.at(0).Momentum_.Phi() - outputAK5_PtSorted.at(1).Momentum_.Phi()) ;

    if(outputAK5_PtSorted.at(0).NameCollection_ == "JetPFCor" && outputAK5_PtSorted.at(1).NameCollection_ == "JetPFCor"){
      NewMuonTree->vbf_maxpt_jj_type = 1 ;

      int nexcj = 0 , nexfj = 0; 
      std::vector<JetAK5>::const_iterator itVec = outputAK5_PtSorted.begin();
      for( ; itVec != outputAK5_PtSorted.end() ; itVec++){

	if(itVec->NameCollection_ == "JetPFCor") nexcj ++ ;
   	if(itVec->NameCollection_ == "JetPFCorVBFTag") nexfj ++ ;
   
      }

      NewMuonTree->vbf_maxpt_n_excj = nexcj ;
      NewMuonTree->vbf_maxpt_n_exfj = nexfj ;
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
    nStep = 10;

    if( NewMuonTree->vbf_maxpt_jj_type < 0 || NewMuonTree->vbf_maxpt_n_excj < 0 || NewMuonTree->vbf_maxpt_n_exfj < 0 ) continue ;
    
    if(outputAK5_PtSorted.at(0).NameCollection_ == "JetPFCor") {
    
      NewMuonTree->vbf_maxpt_j1_QGLikelihood = MuonTree->fReader->getFloat("JetPFCor_QGLikelihood")[outputAK5_PtSorted.at(0).position_] ;

      NewMuonTree->vbf_maxpt_j1_isPileUpLoose  = MuonTree->fReader->getFloat("JetPFCor_isPileUpJetLoose")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_isPileUpMedium = MuonTree->fReader->getFloat("JetPFCor_isPileUpJetMedium")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_isPileUpTight  = MuonTree->fReader->getFloat("JetPFCor_isPileUpJetTight")[outputAK5_PtSorted.at(0).position_] ;


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

      NewMuonTree->vbf_maxpt_j1_isPileUpLoose  = MuonTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetLoose")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_isPileUpMedium = MuonTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetMedium")[outputAK5_PtSorted.at(0).position_] ;
      NewMuonTree->vbf_maxpt_j1_isPileUpTight  = MuonTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetTight")[outputAK5_PtSorted.at(0).position_] ;


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

      NewMuonTree->vbf_maxpt_j2_isPileUpLoose  = MuonTree->fReader->getFloat("JetPFCor_isPileUpJetLoose")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_isPileUpMedium = MuonTree->fReader->getFloat("JetPFCor_isPileUpJetMedium")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_isPileUpTight  = MuonTree->fReader->getFloat("JetPFCor_isPileUpJetTight")[outputAK5_PtSorted.at(1).position_] ;


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

      NewMuonTree->vbf_maxpt_j2_isPileUpLoose  = MuonTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetLoose")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_isPileUpMedium = MuonTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetMedium")[outputAK5_PtSorted.at(1).position_] ;
      NewMuonTree->vbf_maxpt_j2_isPileUpTight  = MuonTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetTight")[outputAK5_PtSorted.at(1).position_] ;


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
    NewMuonTree->vbf_maxDeta_jj_deta = fabs(outputAK5_DEtaSorted.at(0).Momentum_.Phi() - outputAK5_DEtaSorted.at(1).Momentum_.Phi()) ;

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
    nStep = 11;

    if(outputAK5_DEtaSorted.at(0).NameCollection_ == "JetPFCor") {
    
      NewMuonTree->vbf_maxDeta_j1_QGLikelihood = MuonTree->fReader->getFloat("JetPFCor_QGLikelihood")[outputAK5_DEtaSorted.at(0).position_] ;

      NewMuonTree->vbf_maxDeta_j1_isPileUpLoose  = MuonTree->fReader->getFloat("JetPFCor_isPileUpJetLoose")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_isPileUpMedium = MuonTree->fReader->getFloat("JetPFCor_isPileUpJetMedium")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_isPileUpTight  = MuonTree->fReader->getFloat("JetPFCor_isPileUpJetTight")[outputAK5_DEtaSorted.at(0).position_] ;


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

      NewMuonTree->vbf_maxDeta_j1_isPileUpLoose  = MuonTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetLoose")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_isPileUpMedium = MuonTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetMedium")[outputAK5_DEtaSorted.at(0).position_] ;
      NewMuonTree->vbf_maxDeta_j1_isPileUpTight  = MuonTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetTight")[outputAK5_DEtaSorted.at(0).position_] ;


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

      NewMuonTree->vbf_maxDeta_j2_isPileUpLoose  = MuonTree->fReader->getFloat("JetPFCor_isPileUpJetLoose")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_isPileUpMedium = MuonTree->fReader->getFloat("JetPFCor_isPileUpJetMedium")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_isPileUpTight  = MuonTree->fReader->getFloat("JetPFCor_isPileUpJetTight")[outputAK5_DEtaSorted.at(1).position_] ;


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

      NewMuonTree->vbf_maxDeta_j2_isPileUpLoose = MuonTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetLoose")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_isPileUpMedium = MuonTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetMedium")[outputAK5_DEtaSorted.at(1).position_] ;
      NewMuonTree->vbf_maxDeta_j2_isPileUpTight = MuonTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetTight")[outputAK5_DEtaSorted.at(1).position_] ;


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
    NewMuonTree->vbf_maxMjj_jj_deta = fabs(outputAK5_MjjSorted.at(0).Momentum_.Phi() - outputAK5_MjjSorted.at(1).Momentum_.Phi()) ;

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

    if(outputAK5_MjjSorted.at(0).NameCollection_ == "JetPFCor") {
    
      NewMuonTree->vbf_maxMjj_j1_QGLikelihood = MuonTree->fReader->getFloat("JetPFCor_QGLikelihood")[outputAK5_MjjSorted.at(0).position_] ;

      NewMuonTree->vbf_maxMjj_j1_isPileUpLoose  = MuonTree->fReader->getFloat("JetPFCor_isPileUpJetLoose")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_isPileUpMedium = MuonTree->fReader->getFloat("JetPFCor_isPileUpJetMedium")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_isPileUpTight  = MuonTree->fReader->getFloat("JetPFCor_isPileUpJetTight")[outputAK5_MjjSorted.at(0).position_] ;


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

      NewMuonTree->vbf_maxMjj_j1_isPileUpLoose  = MuonTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetLoose")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_isPileUpMedium = MuonTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetMedium")[outputAK5_MjjSorted.at(0).position_] ;
      NewMuonTree->vbf_maxMjj_j1_isPileUpTight  = MuonTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetTight")[outputAK5_MjjSorted.at(0).position_] ;


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

      NewMuonTree->vbf_maxMjj_j2_isPileUpLoose = MuonTree->fReader->getFloat("JetPFCor_isPileUpJetLoose")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_isPileUpMedium = MuonTree->fReader->getFloat("JetPFCor_isPileUpJetMedium")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_isPileUpTight = MuonTree->fReader->getFloat("JetPFCor_isPileUpJetTight")[outputAK5_MjjSorted.at(1).position_] ;


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

      NewMuonTree->vbf_maxMjj_j2_isPileUpLoose  = MuonTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetLoose")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_isPileUpMedium = MuonTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetMedium")[outputAK5_MjjSorted.at(1).position_] ;
      NewMuonTree->vbf_maxMjj_j2_isPileUpTight  = MuonTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetTight")[outputAK5_MjjSorted.at(1).position_] ;


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
    
    NewMuonTree->fTree->Fill();
    
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
   std::cout<<"                                "<<std::endl;

   // create and open the output file 

   TFile *outputFile = new TFile((OutputRootDirectory+"/El"+OutputRootFile).c_str(),"RECREATE");
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

    // Control the Jet Multiplicity for VBF regime
    if( (ElectronTree->fReader->getInt("numPFCorJets")[0] +ElectronTree->fReader->getInt("numPFCorVBFTagJets")[0]) < NumJetMin ) continue ;

    if(std::string(SelectionEvents->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEvents->GetXaxis()->SetBinLabel(nStep,"Jet Number"); 
    if(std::string(SelectionEfficiency->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEfficiency->GetXaxis()->SetBinLabel(nStep,"Jet Number"); 
    nstepEvents[nStep-1]++;
    nStep = 4;

    // Join forward and central PFJetCor Collection    
    std::vector<TLorentzVector> GroomedJet_CA8_Collection ;

    std::vector<JetAK5> JetPFCor_AK5_Collection ;
    std::vector<JetAK5> CleanedJetPFCor_AK5_Collection ;
    std::vector<JetAK5> HadronicW_AK5_Collection ;
    
    for(int iJet = 0 ; iJet < JetCollectionDimension ; iJet++) {
     TLorentzVector JetTemp ; std::string nameCollection ; 
     JetTemp.SetPtEtaPhiE(ElectronTree->fReader->getFloat("GroomedJet_CA8_pt")[iJet],ElectronTree->fReader->getFloat("GroomedJet_CA8_eta")[iJet], 
                          ElectronTree->fReader->getFloat("GroomedJet_CA8_phi")[iJet],ElectronTree->fReader->getFloat("GroomedJet_CA8_e")[iJet]);

     // Selection on CA8 Jets
     if(fabs(JetTemp.Eta())<JetEtaCutMax && JetTemp.Pt()>JetPtCutMin)
       GroomedJet_CA8_Collection.push_back(JetTemp);

     JetTemp.SetPtEtaPhiE(ElectronTree->fReader->getFloat("JetPFCor_Pt")[iJet],ElectronTree->fReader->getFloat("JetPFCor_Eta")[iJet],
 			  ElectronTree->fReader->getFloat("JetPFCor_Phi")[iJet],ElectronTree->fReader->getFloat("JetPFCor_E")[iJet]);
     // Selection on PF Cor Central Jets --> AK5
     if(fabs(JetTemp.Eta())<JetEtaCutMax && JetTemp.Pt()>JetPtCutMin){
       JetAK5 tempJetAK5 (iJet,"JetPFCor",JetTemp);
       JetPFCor_AK5_Collection.push_back(tempJetAK5);
     }
     
    }
   
    for(int iJet = 0 ; iJet < JetCollectionDimension ; iJet++) {
     TLorentzVector JetTemp ; 
     JetTemp.SetPtEtaPhiE(ElectronTree->fReader->getFloat("JetPFCorVBFTag_Pt")[iJet],ElectronTree->fReader->getFloat("JetPFCorVBFTag_Eta")[iJet], 
                          ElectronTree->fReader->getFloat("JetPFCorVBFTag_Phi")[iJet],ElectronTree->fReader->getFloat("JetPFCorVBFTag_E")[iJet]);

     // Selection on PF Cor Forward Jets --> AK5
     if(fabs(JetTemp.Eta())<JetEtaCutMax && JetTemp.Pt()>JetPtCutMin){
       JetAK5 tempJetAK5 (iJet,"JetPFCorVBFTag",JetTemp);
       JetPFCor_AK5_Collection.push_back(tempJetAK5);
     }
    }

    // Another Skim on Jet Counting
    if(GroomedJet_CA8_Collection.size() == 0 || JetPFCor_AK5_Collection.size() < 3) continue ;
   
    if(std::string(SelectionEvents->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEvents->GetXaxis()->SetBinLabel(nStep,"Jet Number CA8 - AK5");
    if(std::string(SelectionEfficiency->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEfficiency->GetXaxis()->SetBinLabel(nStep,"Jet Number CA8 - AK5"); 
    nstepEvents[nStep-1]++;
    nStep = 5;

    // Clean AK5 Jet Collection choosing the hard CA8 as W Hadronic

    for(size_t iJet = 0; iJet < JetPFCor_AK5_Collection.size() ; iJet ++){

      if(deltaR(JetPFCor_AK5_Collection.at(iJet).Momentum_.Phi(),GroomedJet_CA8_Collection.at(0).Phi(),
                JetPFCor_AK5_Collection.at(iJet).Momentum_.Eta(),GroomedJet_CA8_Collection.at(0).Eta()) < CleaningTreshold ){
	HadronicW_AK5_Collection.push_back(JetPFCor_AK5_Collection.at(iJet)); continue ;}

      CleanedJetPFCor_AK5_Collection.push_back(JetPFCor_AK5_Collection.at(iJet));

    }

    if(CleanedJetPFCor_AK5_Collection.size() < 2) continue ; 

    if(std::string(SelectionEvents->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEvents->GetXaxis()->SetBinLabel(nStep,"Matching CA8 - AK5");
    if(std::string(SelectionEfficiency->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEfficiency->GetXaxis()->SetBinLabel(nStep,"Matching CA8 - AK5");
    nstepEvents[nStep-1]++;
    nStep = 6 ;
    // vbf Tag Jet Selection
    
    std::vector<JetAK5> outputAK5_PtSorted;
    std::vector<JetAK5> outputAK5_DEtaSorted;
    std::vector<JetAK5> outputAK5_MjjSorted;

    // Sorting of AK5 Cleaned Collection in Pt

    std::sort(CleanedJetPFCor_AK5_Collection.begin(),CleanedJetPFCor_AK5_Collection.end(),TLVP_PtSort());
    outputAK5_PtSorted.push_back(CleanedJetPFCor_AK5_Collection.at(0));
    outputAK5_PtSorted.push_back(CleanedJetPFCor_AK5_Collection.at(1));
    if(outputAK5_PtSorted.size() < 2) continue ;
    
    // Sorting of AK5 Cleaned Collection in DeltaEta

    std::sort(CleanedJetPFCor_AK5_Collection.begin(),CleanedJetPFCor_AK5_Collection.end(),TLVP_EtaSort());
    outputAK5_DEtaSorted.push_back(CleanedJetPFCor_AK5_Collection.front());
    outputAK5_DEtaSorted.push_back(CleanedJetPFCor_AK5_Collection.back());
    if(outputAK5_DEtaSorted.size() < 2) continue ;

    // Sorting of AK5 Cleaned Collection in Mjj
    float maxMjj = 0. ;
    int iJ1 = 0 ;
    int iJ2 = 0 ;

    for (size_t iJet = 0 ; iJet < CleanedJetPFCor_AK5_Collection.size()-1 ; ++iJet){
      for (size_t jJet = iJet + 1 ; jJet < CleanedJetPFCor_AK5_Collection.size() ; ++jJet){

	TLorentzVector SumMomentum = CleanedJetPFCor_AK5_Collection.at(iJet).Momentum_ + CleanedJetPFCor_AK5_Collection.at(jJet).Momentum_ ;
        float Mjj = SumMomentum.M();
        if(Mjj > maxMjj){
         
          iJ1 = iJet ;
          iJ2 = jJet ;
          maxMjj = Mjj ;
        }
      }
    }

    outputAK5_MjjSorted.push_back (CleanedJetPFCor_AK5_Collection.at (iJ1)) ;
    outputAK5_MjjSorted.push_back (CleanedJetPFCor_AK5_Collection.at (iJ2)) ;
    if(outputAK5_MjjSorted.size() < 2) continue ;

    // Calculate Neutrino Pz
    TLorentzVector W_electron, W_Met, W_neutrino; 
   
    W_electron.SetPxPyPzE(ElectronTree->fReader->getFloat("W_electron_px")[0],ElectronTree->fReader->getFloat("W_electron_py")[0],
                          ElectronTree->fReader->getFloat("W_electron_pz")[0],ElectronTree->fReader->getFloat("W_electron_e")[0]);
    W_Met.SetPxPyPzE(ElectronTree->fReader->getFloat("event_met_pfmet")[0] * TMath::Cos(ElectronTree->fReader->getFloat("event_met_pfmetPhi")[0]),
                     ElectronTree->fReader->getFloat("event_met_pfmet")[0] * TMath::Sin(ElectronTree->fReader->getFloat("event_met_pfmetPhi")[0]),0.,
                     fabs(ElectronTree->fReader->getFloat("event_met_pfmet")[0]));

    if(W_electron.Pt()<=0 || W_Met.Pt() <= 0 ){ std::cerr<<" Negative Lepton - Neutrino Pt "<<std::endl; continue ; }

    if(std::string(SelectionEvents->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEvents->GetXaxis()->SetBinLabel(nStep,"Leptonic 4V");
    if(std::string(SelectionEfficiency->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEfficiency->GetXaxis()->SetBinLabel(nStep,"Leptonic 4V");
    nstepEvents[nStep-1]++;
    nStep = 7;

    METzCalculator<TLorentzVector> NeutrinoPz;
    NeutrinoPz.SetMET(W_Met);
    NeutrinoPz.SetLepton(W_electron);
    NeutrinoPz.SetLeptonType("electron");
    double pz = NeutrinoPz.Calculate(); // Default one
    W_neutrino.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz,sqrt(W_Met.Pt()*W_Met.Pt()+pz*pz));
    if (NeutrinoPz.IsComplex()) {// if this is a complix, change MET                                                                                                                 
           
     double nu_pt1 = NeutrinoPz.getPtneutrino(1);
     double nu_pt2 = NeutrinoPz.getPtneutrino(2);
   
     TLorentzVector W_neutrino_1;
     W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(ElectronTree->fReader->getFloat("event_met_pfmetPhi")[0]), 
                             nu_pt1 * TMath::Sin(ElectronTree->fReader->getFloat("event_met_pfmetPhi")[0]), pz, sqrt(nu_pt1*nu_pt1 + pz*pz) );
     TLorentzVector W_neutrino_2;
     W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(ElectronTree->fReader->getFloat("event_met_pfmetPhi")[0]),
                             nu_pt2 * TMath::Sin(ElectronTree->fReader->getFloat("event_met_pfmetPhi")[0]), pz, sqrt(nu_pt2*nu_pt2 + pz*pz) );

     if ( fabs((W_electron+W_neutrino_1).M()-80.4) < fabs((W_electron+W_neutrino_2).M()-80.4) )  W_neutrino = W_neutrino_1;
     else W_neutrino = W_neutrino_2;

    }
    
    TLorentzVector W_subjet1, W_subjet2 ;  // take the two subjet of the hardest CA8 after pruning
   
    W_subjet1.SetPxPyPzE(ElectronTree->fReader->getFloat("GroomedJet_CA8_prsubjet1_px")[0],ElectronTree->fReader->getFloat("GroomedJet_CA8_prsubjet1_py")[0],
         		 ElectronTree->fReader->getFloat("GroomedJet_CA8_prsubjet1_pz")[0],ElectronTree->fReader->getFloat("GroomedJet_CA8_prsubjet1_e")[0] );
    W_subjet2.SetPxPyPzE(ElectronTree->fReader->getFloat("GroomedJet_CA8_prsubjet2_px")[0],ElectronTree->fReader->getFloat("GroomedJet_CA8_prsubjet2_py")[0],
  		  	 ElectronTree->fReader->getFloat("GroomedJet_CA8_prsubjet2_pz")[0],ElectronTree->fReader->getFloat("GroomedJet_CA8_prsubjet2_e")[0] );

    if(W_subjet1.Pt() <= 0 || W_subjet2.Pt() <= 0){ std::cerr<<" Problem with subjets "<<std::endl; continue ; }

    if(std::string(SelectionEvents->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEvents->GetXaxis()->SetBinLabel(nStep,"Subjets 4V");
    if(std::string(SelectionEfficiency->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEfficiency->GetXaxis()->SetBinLabel(nStep,"Subjets 4V");
    nstepEvents[nStep-1]++;
    nStep = 8;


    TLorentzVector W_GroomedJet_CA8_pr; 
    W_GroomedJet_CA8_pr.SetPtEtaPhiE(ElectronTree->fReader->getFloat("GroomedJet_CA8_pt_pr")[0], ElectronTree->fReader->getFloat("GroomedJet_CA8_eta_pr")[0],
                                     ElectronTree->fReader->getFloat("GroomedJet_CA8_phi_pr")[0], ElectronTree->fReader->getFloat("GroomedJet_CA8_e_pr")[0]);

    if(W_GroomedJet_CA8_pr.Pt() <=0){ std::cerr<<" Problem with pruned CA8 "<<std::endl; continue ;}

    if(std::string(SelectionEvents->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEvents->GetXaxis()->SetBinLabel(nStep,"Pruned CA8 4V");
    if(std::string(SelectionEfficiency->GetXaxis()->GetBinLabel(nStep)) =="") SelectionEfficiency->GetXaxis()->SetBinLabel(nStep,"Pruned CA8 4V");
    nstepEvents[nStep-1]++;
    nStep = 9;

   
    // Kinematic Fit                                                                                                                                                                 
    TLorentzVector fit_electron(0,0,0,0), fit_neutrino(0,0,0,0), fit_W_subjet1(0,0,0,0), fit_W_subjet2(0,0,0,0) ;

    doKinematicFit(1, W_electron, W_neutrino, W_subjet1, W_subjet2,  fit_electron, fit_neutrino, fit_W_subjet1, fit_W_subjet2, NewElectronTree->fit_chi2, 
                   NewElectronTree->fit_NDF, NewElectronTree->fit_status, LeptonType);
    
    if(fit_electron.Pt() >0 && fit_neutrino.Pt()>0 && fit_W_subjet1.Pt()>0 && fit_W_subjet2.Pt()>0){

    NewElectronTree->fit_el_px = fit_electron.Px();
    NewElectronTree->fit_el_py = fit_electron.Py(); 
    NewElectronTree->fit_el_pz = fit_electron.Pz(); 
    NewElectronTree->fit_el_e  = fit_electron.E();
    NewElectronTree->fit_nv_px = fit_neutrino.Px(); 
    NewElectronTree->fit_nv_py = fit_neutrino.Py(); 
    NewElectronTree->fit_nv_pz = fit_neutrino.Pz(); 
    NewElectronTree->fit_nv_e  = fit_neutrino.E();

    NewElectronTree->fit_subjet1_px = fit_W_subjet1.Px();  NewElectronTree->fit_subjet2_px = fit_W_subjet2.Px();
    NewElectronTree->fit_subjet1_py = fit_W_subjet1.Py();  NewElectronTree->fit_subjet2_py = fit_W_subjet2.Py();  
    NewElectronTree->fit_subjet1_pz = fit_W_subjet1.Pz();  NewElectronTree->fit_subjet2_pz = fit_W_subjet2.Pz();
    NewElectronTree->fit_subjet1_e  = fit_W_subjet1.E();   NewElectronTree->fit_subjet2_e  = fit_W_subjet2.E();
    NewElectronTree->fit_subjet1_m  = fit_W_subjet1.M();   NewElectronTree->fit_subjet2_m  = fit_W_subjet2.M();

    NewElectronTree->fit_lvj_m   = (W_electron+W_neutrino+fit_W_subjet1+fit_W_subjet2).M();
    NewElectronTree->fit_lv_m    = (W_electron+W_neutrino).M();
    NewElectronTree->fit_j_m     = (fit_W_subjet1+fit_W_subjet2).M();
    NewElectronTree->fit_lvj_pt  = (fit_electron+fit_neutrino+fit_W_subjet1+fit_W_subjet2).M();
    NewElectronTree->fit_lvj_eta = (fit_electron+fit_neutrino+fit_W_subjet1+fit_W_subjet2).Eta();
    NewElectronTree->fit_lvj_phi = (fit_electron+fit_neutrino+fit_W_subjet1+fit_W_subjet2).Phi();
    NewElectronTree->fit_lvj_e   = (fit_electron+fit_neutrino+fit_W_subjet1+fit_W_subjet2).E();

    NewElectronTree->boosted_lvj_m   = (W_electron+W_neutrino+W_subjet1+W_subjet2).M();
    NewElectronTree->boosted_lv_m    = (W_electron+W_neutrino).M();
    NewElectronTree->boosted_j_m     = (W_subjet1+W_subjet2).M();
    NewElectronTree->boosted_lvj_pt  = (W_electron+W_neutrino+W_subjet1+W_subjet2).M ();
    NewElectronTree->boosted_lvj_eta = (W_electron+W_neutrino+W_subjet1+W_subjet2).Eta();
    NewElectronTree->boosted_lvj_phi = (W_electron+W_neutrino+W_subjet1+W_subjet2).Phi();
    NewElectronTree->boosted_lvj_e   = (W_electron+W_neutrino+W_subjet1+W_subjet2).E();

    NewElectronTree->boostedW_lvj_m   = (W_electron+W_neutrino+W_GroomedJet_CA8_pr).M();
    NewElectronTree->boostedW_lv_m    = (W_electron+W_neutrino).M();
    NewElectronTree->boostedW_j_m     = (W_GroomedJet_CA8_pr).M();
    NewElectronTree->boostedW_lvj_pt  = (W_electron+W_neutrino+W_GroomedJet_CA8_pr).M ();
    NewElectronTree->boostedW_lvj_eta = (W_electron+W_neutrino+W_GroomedJet_CA8_pr).Eta();
    NewElectronTree->boostedW_lvj_phi = (W_electron+W_neutrino+W_GroomedJet_CA8_pr).Phi();
    NewElectronTree->boostedW_lvj_e   = (W_electron+W_neutrino+W_GroomedJet_CA8_pr).E();
    // Angles for the central Higgs Kinematics
    }

    double costheta1, costheta2, phi, costhetastar, phistar1, phistar2;

    //Use the Subjet in the Boosted W Analyisis                                                                                                                                             
    if (ElectronTree->fReader->getFloat("W_electron_charge")[0] < 0) calculateAngles(W_electron, W_neutrino,W_subjet1,W_subjet2,costheta1,costheta2,phi,costhetastar,phistar1,phistar2);
    else calculateAngles(W_neutrino, W_electron, W_subjet1, W_subjet2, costheta1, costheta2, phi, costhetastar, phistar1, phistar2);
   
    NewElectronTree->boosted_wjj_ang_ha   = costheta1;
    NewElectronTree->boosted_wjj_ang_hb   = fabs(costheta2); 
    NewElectronTree->boosted_wjj_ang_hs   = costhetastar;
    NewElectronTree->boosted_wjj_ang_phi  = phi;
    NewElectronTree->boosted_wjj_ang_phia = phistar1;																		    NewElectronTree->boosted_wjj_ang_phib = phistar2;
   
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
    NewElectronTree->vbf_maxpt_jj_deta = fabs(outputAK5_PtSorted.at(0).Momentum_.Phi() - outputAK5_PtSorted.at(1).Momentum_.Phi()) ;

    if(outputAK5_PtSorted.at(0).NameCollection_ == "JetPFCor" && outputAK5_PtSorted.at(1).NameCollection_ == "JetPFCor"){
      NewElectronTree->vbf_maxpt_jj_type = 1 ;

      int nexcj = 0 , nexfj = 0; 
      std::vector<JetAK5>::const_iterator itVec = outputAK5_PtSorted.begin();
      for( ; itVec != outputAK5_PtSorted.end() ; itVec++){

	if(itVec->NameCollection_ == "JetPFCor") nexcj ++ ;
   	if(itVec->NameCollection_ == "JetPFCorVBFTag") nexfj ++ ;
   
      }

      NewElectronTree->vbf_maxpt_n_excj = nexcj ;
      NewElectronTree->vbf_maxpt_n_exfj = nexfj ;
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
    nStep = 10;

    if( NewElectronTree->vbf_maxpt_jj_type < 0 || NewElectronTree->vbf_maxpt_n_excj < 0 || NewElectronTree->vbf_maxpt_n_exfj < 0 ) continue ;
    
    if(outputAK5_PtSorted.at(0).NameCollection_ == "JetPFCor") {
    
      NewElectronTree->vbf_maxpt_j1_QGLikelihood = ElectronTree->fReader->getFloat("JetPFCor_QGLikelihood")[outputAK5_PtSorted.at(0).position_] ;

      NewElectronTree->vbf_maxpt_j1_isPileUpLoose  = ElectronTree->fReader->getFloat("JetPFCor_isPileUpJetLoose")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_isPileUpMedium = ElectronTree->fReader->getFloat("JetPFCor_isPileUpJetMedium")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_isPileUpTight  = ElectronTree->fReader->getFloat("JetPFCor_isPileUpJetTight")[outputAK5_PtSorted.at(0).position_] ;


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

      NewElectronTree->vbf_maxpt_j1_isPileUpLoose  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetLoose")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_isPileUpMedium = ElectronTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetMedium")[outputAK5_PtSorted.at(0).position_] ;
      NewElectronTree->vbf_maxpt_j1_isPileUpTight  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetTight")[outputAK5_PtSorted.at(0).position_] ;


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

      NewElectronTree->vbf_maxpt_j2_isPileUpLoose  = ElectronTree->fReader->getFloat("JetPFCor_isPileUpJetLoose")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_isPileUpMedium = ElectronTree->fReader->getFloat("JetPFCor_isPileUpJetMedium")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_isPileUpTight  = ElectronTree->fReader->getFloat("JetPFCor_isPileUpJetTight")[outputAK5_PtSorted.at(1).position_] ;


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

      NewElectronTree->vbf_maxpt_j2_isPileUpLoose  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetLoose")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_isPileUpMedium = ElectronTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetMedium")[outputAK5_PtSorted.at(1).position_] ;
      NewElectronTree->vbf_maxpt_j2_isPileUpTight  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetTight")[outputAK5_PtSorted.at(1).position_] ;


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
    NewElectronTree->vbf_maxDeta_jj_deta = fabs(outputAK5_DEtaSorted.at(0).Momentum_.Phi() - outputAK5_DEtaSorted.at(1).Momentum_.Phi()) ;

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
    nStep = 11;

    if(outputAK5_DEtaSorted.at(0).NameCollection_ == "JetPFCor") {
    
      NewElectronTree->vbf_maxDeta_j1_QGLikelihood = ElectronTree->fReader->getFloat("JetPFCor_QGLikelihood")[outputAK5_DEtaSorted.at(0).position_] ;

      NewElectronTree->vbf_maxDeta_j1_isPileUpLoose  = ElectronTree->fReader->getFloat("JetPFCor_isPileUpJetLoose")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_isPileUpMedium = ElectronTree->fReader->getFloat("JetPFCor_isPileUpJetMedium")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_isPileUpTight  = ElectronTree->fReader->getFloat("JetPFCor_isPileUpJetTight")[outputAK5_DEtaSorted.at(0).position_] ;


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

      NewElectronTree->vbf_maxDeta_j1_isPileUpLoose  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetLoose")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_isPileUpMedium = ElectronTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetMedium")[outputAK5_DEtaSorted.at(0).position_] ;
      NewElectronTree->vbf_maxDeta_j1_isPileUpTight  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetTight")[outputAK5_DEtaSorted.at(0).position_] ;


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

      NewElectronTree->vbf_maxDeta_j2_isPileUpLoose  = ElectronTree->fReader->getFloat("JetPFCor_isPileUpJetLoose")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_isPileUpMedium = ElectronTree->fReader->getFloat("JetPFCor_isPileUpJetMedium")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_isPileUpTight  = ElectronTree->fReader->getFloat("JetPFCor_isPileUpJetTight")[outputAK5_DEtaSorted.at(1).position_] ;


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

      NewElectronTree->vbf_maxDeta_j2_isPileUpLoose = ElectronTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetLoose")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_isPileUpMedium = ElectronTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetMedium")[outputAK5_DEtaSorted.at(1).position_] ;
      NewElectronTree->vbf_maxDeta_j2_isPileUpTight = ElectronTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetTight")[outputAK5_DEtaSorted.at(1).position_] ;


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
    NewElectronTree->vbf_maxMjj_jj_deta = fabs(outputAK5_MjjSorted.at(0).Momentum_.Phi() - outputAK5_MjjSorted.at(1).Momentum_.Phi()) ;

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

    if(outputAK5_MjjSorted.at(0).NameCollection_ == "JetPFCor") {
    
      NewElectronTree->vbf_maxMjj_j1_QGLikelihood = ElectronTree->fReader->getFloat("JetPFCor_QGLikelihood")[outputAK5_MjjSorted.at(0).position_] ;

      NewElectronTree->vbf_maxMjj_j1_isPileUpLoose  = ElectronTree->fReader->getFloat("JetPFCor_isPileUpJetLoose")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_isPileUpMedium = ElectronTree->fReader->getFloat("JetPFCor_isPileUpJetMedium")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_isPileUpTight  = ElectronTree->fReader->getFloat("JetPFCor_isPileUpJetTight")[outputAK5_MjjSorted.at(0).position_] ;


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

      NewElectronTree->vbf_maxMjj_j1_isPileUpLoose  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetLoose")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_isPileUpMedium = ElectronTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetMedium")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_isPileUpTight  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetTight")[outputAK5_MjjSorted.at(0).position_] ;


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
      NewElectronTree->vbf_maxMjj_j1_MuonMultiplicity           = ElectronTree->fReader->getFloat("JetPFCorVBFTag_MuonMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_ChargedHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ChargedHadronMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_NeutralHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_NeutralHadronMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_PhotonMultiplicity         = ElectronTree->fReader->getFloat("JetPFCorVBFTag_PhotonMultiplicity")[outputAK5_MjjSorted.at(0).position_] ; 
      NewElectronTree->vbf_maxMjj_j1_ElectronMultiplicity       = ElectronTree->fReader->getFloat("JetPFCorVBFTag_ElectronMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;
      NewElectronTree->vbf_maxMjj_j1_HFHadronMultiplicity       = ElectronTree->fReader->getFloat("JetPFCorVBFTag_HFHadronMultiplicity")[outputAK5_MjjSorted.at(0).position_] ;  

    }
    else { std::cerr<<" problem with High Deta Jet Name Collection "<<std::endl; continue ; }

    if(outputAK5_MjjSorted.at(1).NameCollection_ == "JetPFCor") {
    
      NewElectronTree->vbf_maxMjj_j2_QGLikelihood = ElectronTree->fReader->getFloat("JetPFCor_QGLikelihood")[outputAK5_MjjSorted.at(1).position_] ;

      NewElectronTree->vbf_maxMjj_j2_isPileUpLoose = ElectronTree->fReader->getFloat("JetPFCor_isPileUpJetLoose")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_isPileUpMedium = ElectronTree->fReader->getFloat("JetPFCor_isPileUpJetMedium")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_isPileUpTight = ElectronTree->fReader->getFloat("JetPFCor_isPileUpJetTight")[outputAK5_MjjSorted.at(1).position_] ;


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
      NewElectronTree->vbf_maxMjj_j2_MuonMultiplicity           = ElectronTree->fReader->getFloat("JetPFCor_MuonMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_ChargedHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCor_ChargedHadronMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_NeutralHadronMultiplicity  = ElectronTree->fReader->getFloat("JetPFCor_NeutralHadronMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_PhotonMultiplicity         = ElectronTree->fReader->getFloat("JetPFCor_PhotonMultiplicity")[outputAK5_MjjSorted.at(1).position_] ; 
      NewElectronTree->vbf_maxMjj_j2_ElectronMultiplicity       = ElectronTree->fReader->getFloat("JetPFCor_ElectronMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_HFHadronMultiplicity       = ElectronTree->fReader->getFloat("JetPFCor_HFHadronMultiplicity")[outputAK5_MjjSorted.at(1).position_] ;
    }
    else if(outputAK5_MjjSorted.at(1).NameCollection_ == "JetPFCorVBFTag" ){

      NewElectronTree->vbf_maxMjj_j2_QGLikelihood = ElectronTree->fReader->getFloat("JetPFCorVBFTag_QGLikelihood")[outputAK5_MjjSorted.at(1).position_] ;

      NewElectronTree->vbf_maxMjj_j2_isPileUpLoose  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetLoose")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_isPileUpMedium = ElectronTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetMedium")[outputAK5_MjjSorted.at(1).position_] ;
      NewElectronTree->vbf_maxMjj_j2_isPileUpTight  = ElectronTree->fReader->getFloat("JetPFCorVBFTag_isPileUpJetTight")[outputAK5_MjjSorted.at(1).position_] ;


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
    
    NewElectronTree->fTree->Fill();
    
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
  
