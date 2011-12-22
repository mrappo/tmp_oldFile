///============== Ntuple Producer at Gen particle level 
///======== lunch with 
///== ls test/Latinos/dir_cfg_skimmed_MC_GenLevel/ | grep NtupleProducerNT | awk '{ print" ./bin/NtupleProducerFromNT_Gen.exe test/Latinos/dir_cfg_skimmed_MC_GenLevel/"$1" \n "}' | /bin/sh

#include "treeReader.h"
#include "hFactory.h"
#include "hFunctions.h"
#include "stdHisto.h"
#include "ConfigParser.h"
#include "ntpleUtils.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMinuit.h"
#include "Math/GenVector/VectorUtil.h"
#include "TRandom3.h"
#include <time.h>
#include <sstream>
#include "qqHWWlnulnuUtils.h"
#include "Variables_Gen.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#endif

int GetNumList(std::vector<int> &list){
 int result = 0;
 for (int it = 0; it<list.size(); it++) result += list.at(it);
 return result;
}


int main(int argc, char** argv)
{ 
 //Check if all nedeed arguments to parse are there                                                                                                                               
 if(argc != 2)
 {
  std::cerr << ">>>>> analysis.cpp::usage: " << argv[0] << " configFileName" << std::endl ;
  return 1;
 }
 
 /// -------- Parse the config file -----------                                                                                                                                                          
 parseConfigFile (argv[1]) ;
 
 /// --------- tree name, input file acquisition ------ 
 
 std::string treeName  = gConfigParser -> readStringOption("Input::treeName");
 std::string inputFile = gConfigParser -> readStringOption("Input::inputFile"); 

 ///----- Cross section + efficeincy variables
 
 double XSection = gConfigParser -> readDoubleOption("Options::XSection");
 double XSectionErrorUp = gConfigParser -> readDoubleOption("Options::XSectionErrorUp");
 double XSectionErrorDown = gConfigParser -> readDoubleOption("Options::XSectionErrorDown");
 
 std::string histoNameEvents      = gConfigParser -> readStringOption("Input::histoNameEvents"); 
 
 std::cout << ">>>>> Input::inputFile                 " << inputFile  << std::endl;  
 std::cout << ">>>>> Input::histoNameEvents      " << histoNameEvents  << std::endl;  
 std::cout<<">>>>> Options::XSection  " << XSection << std::endl;
 std::cout<<">>>>> Options::XSectionErrorUp  " << XSectionErrorUp << std::endl;
 std::cout<<">>>>> Options::XSectionErrorDown  " << XSectionErrorDown << std::endl;
 
 int entryMAX = gConfigParser -> readIntOption("Input::entryMAX");
 int entryMIN = gConfigParser -> readIntOption("Input::entryMIN");
 int entryMOD = gConfigParser -> readIntOption("Input::entryMOD");
 
 
 std::cout << ">>>>> input::entryMIN  " << entryMIN  << std::endl;  
 std::cout << ">>>>> input::entryMAX  " << entryMAX  << std::endl;  
 std::cout << ">>>>> input::entryMOD  " << entryMOD  << std::endl;  
 
 /// -------- Check number of step to do  ---------
 int nStepToDo = 1000;
 try {
  nStepToDo = gConfigParser -> readIntOption("Input::nStepToDo");
 }
 catch (char const* exceptionString){
  std::cerr << " exception = " << exceptionString << std::endl;
  nStepToDo = 1000;
 }
 std::cout << ">>>>> input::nStepToDo  " << nStepToDo  << std::endl;  
 
 
 /// ------- Open ntple and save it in a TreeReader object, one ntuple for time -----
 
 TChain* chain = new TChain(treeName.c_str());
 chain->Add(inputFile.c_str());
 treeReader reader((TTree*)(chain));
 
 // Open ntples
 TFile File(inputFile.c_str()) ; 
 TH1F* histoEvents = (TH1F*) File.Get(TString(histoNameEvents.c_str()));

 
 /// ---- Acquisition of debug variable ---
 bool  debug = false; 
 try {
  debug = gConfigParser -> readBoolOption("Input::debug");
 }
 catch (char const* exceptionString){
  std::cerr << " exception = " << exceptionString << std::endl;
 }
 std::cout << ">>>>> input::debug  " << debug  << std::endl;  
 
 ///----  Aquisition output file name -------
 std::string outFileName    = gConfigParser -> readStringOption("Output::outFileName");
 std::cout << ">>>>> Output::outFileName  " << outFileName  << std::endl;  
 
  
 /// ----- Definition of gen Variable container 

 int numEntriesBefore;
 double preselection_efficiency = 1.;

 
 Variables_Gen vars;
 InitializeTree(vars, outFileName);
 
 if (entryMAX == -1) entryMAX = reader.GetEntries();
 else if (reader.GetEntries() < entryMAX) entryMAX = reader.GetEntries();
 numEntriesBefore = entryMAX - entryMIN;
 
 if (histoEvents) preselection_efficiency = numEntriesBefore / (1. * histoEvents->GetBinContent(1));
 else preselection_efficiency = 1;

  
 vars.XSection = XSection;
 vars.numEntriesBefore = numEntriesBefore;
 vars.preselection_efficiency = preselection_efficiency;
 vars.XSectionErrorUp = XSectionErrorUp;
 vars.XSectionErrorDown = XSectionErrorDown;
 
 
 FillEfficiencyTree(vars);
 
 double start, end;
 std::cout << ">>>>> analysis::entryMIN " << entryMIN << " ==> entryMAX " << entryMAX << ":" << reader.GetEntries() << std::endl;   
  
 int step = 0;
 start = clock();
 for(int iEvent = entryMIN ; iEvent < entryMAX ; ++iEvent) {
 
  reader.GetEntry(iEvent);
  if((iEvent%entryMOD) == 0) std::cout << ">>>>> analysis::GetEntry " << iEvent << " : " << entryMAX - entryMIN << std::endl;   
   
  ///==== define variables ==== 
  std::vector<ROOT::Math::XYZTVector>* jets = reader.Get4V("mcJet");
  std::vector<ROOT::Math::XYZTVector>* muons = reader.Get4V("mcMu");
  std::vector<ROOT::Math::XYZTVector>* electrons = reader.Get4V("mcEle");
    
  ///*********************************************************************************************
  ///*********************************************************************************************
  
  ///=============================
  ///==== fill MC information for Higgs, Tag quarks, weak bosons and leptons ====
  if (debug) std::cout << " SetMCVariables " << std::endl;
  SetMCVariables(vars, reader);
  
  ///================================
  ///==== fill Event information ====
  if (debug) std::cout << " SetEventVariables " << std::endl;
  SetEventVariables(vars, reader);
  
  
  ///****************************
  ///**** STEP 0 - Ntuplizer ****
  ///************* no additional selections applied 

  step = 0;  
  if (step > nStepToDo) {
   FillTree(vars);
   continue;
  }
  if (debug) std::cout << ">>> STEP 0 <<<" << std::endl;
   
   ///*********************************************
  ///**** STEP 1 - Jet cleaning + min pT ****
  ///************* it's performed another time here to make sure that the cleaning worked well
  ///************* Jet - electrons (pT > 5)
  ///************* Jet - muons     (pT > 5)
  ///************ No selections are applied here
   
  step = 1;
  if (step > nStepToDo) {
   continue;
  }  
  if (debug) std::cout << ">>> STEP 1 <<<" << std::endl;
  
  std::vector<ROOT::Math::XYZTVector> leptons_jetCleaning;   
  
  ///==== CLEANING WITH ELECTRONS ====
  for(unsigned int iEle = 0; iEle < (reader.Get4V("mcEle")->size()); ++iEle)
  {
   if( reader.Get4V("mcEle")->at(iEle).pt() < 5.0 ) continue;
   if( fabs(reader.Get4V("mcEle")->at(iEle).eta()) > 2.5 ) continue;
   
   leptons_jetCleaning.push_back( reader.Get4V("mcEle")->at(iEle) );
  }
  
  ///==== CLEANING WITH MUONS ====
  for (int iMu = 0; iMu < reader.Get4V("mcMu")->size(); iMu++){    
   if (reader.Get4V("mcMu")->at(iMu).pt() < 5.0) continue;
   if (fabs(reader.Get4V("mcMu")->at(iMu).Eta()) > 2.5) continue;
  
   leptons_jetCleaning.push_back( reader.Get4V("mcMu")->at(iMu) );
  }
  
  
  ///==== now clean jet collection ====
  
  int nJets = reader.Get4V("mcJet")->size();
  std::vector<int> whitelistJet; ///~~~~ all jets, 0 if rejected, 1 if accepted
  std::vector<int> blacklistJet; ///~~~~ list of numbers of jets that are "rejected"
  std::vector<int> blacklistJet_forCJV;
  
  
  for (int iJet = 0; iJet < nJets; iJet++){
   bool skipJet = false;   
   if (reader.Get4V("mcJet")->at(iJet).Et() < 15.0) skipJet = true; 
  
   for(unsigned int iLep = 0; iLep < leptons_jetCleaning.size(); ++iLep) {
    ROOT::Math::XYZTVector lep = leptons_jetCleaning.at(iLep);
    if (ROOT::Math::VectorUtil::DeltaR(reader.Get4V("mcJet")->at(iJet),lep) < 0.5 ) skipJet = true;
   }
   if (skipJet) {
    whitelistJet.push_back(0); ///---- reject
    blacklistJet.push_back(iJet); ///---- reject ///== black list is in a different format
    blacklistJet_forCJV.push_back(iJet); ///---- reject ///== black list is in a different format
   }
   else {
         whitelistJet.push_back(1); ///---- select
   }
  }
  
     
     
   ///**************************************
   ///**** STEP 2 - Super-Preselections ****
   ///************* tighter preselections to start the analysis from the same point
   ///==== construct considered objets
   ///    Objects considered and selections
   
   ///   Muon
   ///   Pt>10GeV, eta<2.5
   ///   MuonId & Iso
   ///
   ///   Electron
   ///   Pt>10GeV & |eta|<2.5
   ///   eleId & Iso
   ///
   ///    At least two leptons 
   
   ///   Jet
   ///   Antikt5, L2L3 correction jets
   ///   At least two calo jets or two pf jets with pt>15 GeV
   
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  step = 2;
  if (step > nStepToDo) {
//    FillTree(vars);
   continue;
  }  
  if (debug) std::cout << ">>> STEP 2 <<<" << std::endl;
  
   ///   Electron
   ///   Pt>10GeV & |eta|<2.5
   ///   IsoTr / pTele <0.5
  ///   eleId & Iso
  std::vector<int> whitelistEle;
  std::vector<int> blacklistEle;
  int nEles = reader.Get4V("mcEle")->size();
  for (int iEle = 0; iEle < nEles; iEle++){    
   bool skipEle = false;
   if (reader.Get4V("mcEle")->at(iEle).pt() < 10.0) skipEle = true;
   if (fabs(reader.Get4V("mcEle")->at(iEle).Eta()) > 2.5) skipEle = true;
   
   if (skipEle) {
    whitelistEle.push_back(0); ///---- reject
    blacklistEle.push_back(iEle); ///---- reject ///== black list is in a different format
   }
   else {
    whitelistEle.push_back(1); ///---- select
   }
  }
       
   ///   Muon
   ///   Pt>10GeV, eta<2.5
   ///   MuonId & Iso
   std::vector<int> whitelistMu;
   std::vector<int> blacklistMu;
   int nMus = reader.Get4V("mcMu")->size();
   for (int iMu = 0; iMu < nMus; iMu++){    
    bool skipMu = false;
    if (reader.Get4V("mcMu")->at(iMu).pt() < 10.0) skipMu = true;
    if (fabs(reader.Get4V("mcMu")->at(iMu).Eta()) > 2.5) skipMu = true;    
   
    if (skipMu) {
     whitelistMu.push_back(0); ///---- reject
     blacklistMu.push_back(iMu); ///---- reject ///== black list is in a different format
    }
    else {
     whitelistMu.push_back(1); ///---- select
    }
   }
      
   ///   At least 2 leptons
   
   int numMus_Accepted = GetNumList(whitelistMu);
   int numEles_Accepted = GetNumList(whitelistEle);
   
   int numLeptons_Accepted = numMus_Accepted + numEles_Accepted;
   if (numLeptons_Accepted < 2) continue;
   
   ///   Jet
   
   int numJets_Accepted = GetNumList(whitelistJet);
   if (numJets_Accepted < 2) continue; ///==== at least 2 jets "isolated"
   

  ///*************************
  ///**** STEP 3 - Jet ID ****
  ///************* Identification of two tag jets

  step = 3;
  if (step > nStepToDo) {
//    FillTree(vars);
   continue;
  }  
  if (debug) std::cout << ">>> STEP 3 <<<" << std::endl;
  
  
  std::vector<int> itSelJet;
  double maxPt_jets_selected = SelectJets(itSelJet,*jets,"maxSumPt",-1.,&blacklistJet);
  
  int q1 = itSelJet.at(0);
  int q2 = itSelJet.at(1);
  ///---- check Pt order ----
  if (jets->at(q1).Pt() < jets->at(q2).Pt()) {
   int tempq = q1;
   q1 = q2;
   q2 = tempq;
  }

if (debug) std::cerr << " q1 = " << q1 << " : q2 = " << q2 << std::endl;


  ///---- update white/black list jets ----
  for (int iJet = 0; iJet < nJets; iJet++){
   if (q1 == iJet || q2 == iJet) {
    whitelistJet.at(iJet) = 1;
    blacklistJet.push_back(iJet); ///===>  blacklistJet used for CJV => no 2 tag jets to be considered!
    blacklistJet_forCJV.push_back(iJet); ///===>  blacklistJet used for CJV => no 2 tag jets to be considered!
   }
   else {
    whitelistJet.at(iJet) = 0;
   }
  }
 
 /// Set information on leading jets
 SetQJetVariables(vars, reader, q1, q2, blacklistJet_forCJV);
 
 
 ///=== Third jet search and information setting
 
 std::vector<ROOT::Math::XYZTVector> Jet_Candidate;
 FindAddJet (reader,q1,q2,&blacklistJet_forCJV,Jet_Candidate,0);
 
 SetThirdJetVariables(vars,Jet_Candidate);


  ///********************************
  ///**** STEP 4 - Lepton ID ****
  ///************* Identification of the two leptons
  step = 4;
  if (step > nStepToDo) {
   FillTree(vars);
   continue;
  }  
  if (debug) std::cout << ">>> STEP 4 <<<" << std::endl;
  

  std::vector<ROOT::Math::XYZTVector> leptons;
  std::vector<int> leptonFlavours;    
  std::vector<int> leptonILep;    
  
  for(unsigned int iEle = 0; iEle < reader.Get4V("mcEle")->size(); iEle++){
   if (whitelistEle.at(iEle) == 1){
    leptons.push_back( reader.Get4V("mcEle")->at(iEle) );  
    leptonFlavours.push_back(11);
    leptonILep.push_back(iEle);
   }
  }
  
  for(unsigned int iMu = 0; iMu < reader.Get4V("mcMu")->size(); iMu++){
   if (whitelistMu.at(iMu) == 1){
    leptons.push_back( reader.Get4V("mcMu")->at(iMu) );      
    leptonFlavours.push_back(13);
    leptonILep.push_back(iMu);
   }
  }
  
  std::vector<int> itSelLep;
  double maxPt_lept_selected = SelectJets(itSelLep,leptons,"maxSumPt",-1.,0);
  
  int l1 = itSelLep.at(0);
  int l2 = itSelLep.at(1);
  ///---- check Pt order ----
  if (leptons.at(l1).Pt() < leptons.at(l2).Pt()) {
   int templ = l1;
   l1 = l2;
   l2 = templ;
  }
  
  if (debug) std::cerr << " l1 = " << l1 << " : l2 = " << l2 << " leptonILep.at(l1) " <<leptonILep.at(l1)<< " leptonILep.at(l2) " << leptonILep.at(l2)<<
   " leptonFlavours.at(l1) "  <<leptonFlavours.at(l1)<< " leptonFlavours.at(l2) "<<leptonFlavours.at(l2)<<std::endl;
  
   SetLeptonsVariables(vars, reader, leptonILep.at(l1), leptonILep.at(l2),leptonFlavours.at(l1), leptonFlavours.at(l2));
  
  if (debug) std::cerr << ">> Lepton variables set" << std::endl;
  
  
  SetMetVariables(vars, reader, leptonILep.at(l1), leptonILep.at(l2),leptonFlavours.at(l1), leptonFlavours.at(l2));
  if (debug) std::cerr << ">> MET variables set" << std::endl;
   
   
 
  ///--- MT Higgs ----
  SetMTVariable(vars, reader, leptonILep.at(l1), leptonILep.at(l2),leptonFlavours.at(l1), leptonFlavours.at(l2));
  if (debug) std::cerr << ">> MT Higgs variables set" << std::endl;
   
  ///--- Dphi leptons and jet(s) ----
  SetDPhiJetll(vars,reader, leptonILep.at(l1), leptonILep.at(l2), leptonFlavours.at(l1),leptonFlavours.at(l2),q1,q2);
  if (debug) std::cerr << ">> DPhi Jet-leptons variables set" << std::endl;
   
      
   
  //---- lepton veto
  std::vector<int> blacklistLepton;
  blacklistLepton.push_back(l1);
  blacklistLepton.push_back(l2);
  
  vars.Nleptons_pT5  = getNumberPTThreshold(leptons,  5, &blacklistLepton);
  vars.Nleptons_pT10 = getNumberPTThreshold(leptons, 10, &blacklistLepton);
  vars.Nleptons_pT15 = getNumberPTThreshold(leptons, 15, &blacklistLepton);
  vars.Nleptons_pT20 = getNumberPTThreshold(leptons, 20, &blacklistLepton);
  vars.Nleptons_pT25 = getNumberPTThreshold(leptons, 25, &blacklistLepton);
  vars.Nleptons_pT30 = getNumberPTThreshold(leptons, 30, &blacklistLepton);

  if (debug) std::cerr << ">> Lepton multiplicity set" << std::endl;
  
  
  ///*********************************
  ///**** STEP 5 - Jet Selections ****
  ///************* Loose selections of tag jets
  step = 5;
  if (step > nStepToDo) {
   FillTree(vars);
   continue;
  }  
  if (debug) std::cout << ">>> STEP 5 <<<" << std::endl;
  
  ///==== save Jetdeau variables ====
//   SetJDVariables(vars, reader, leptonILep.at(l1),leptonILep.at(l2),leptonFlavours.at(l1), leptonFlavours.at(l2), q1, q2);
  if (debug) std::cerr << ">> SetJDVariables variables set" << std::endl;

  ///************************************
  ///**** STEP 6 - Final Production *****
  ///************************************
  ///**** No more selections applied ****

  step = 6;
  if (step > nStepToDo) {
   FillTree(vars);
   continue;
  }  
  if (debug) std::cout << ">>> STEP 6 <<<" << std::endl;

  ///==== if not already filled ... ====
  FillTree(vars);  
  ///=================================================
  
 }
 end = clock();
 std::cout <<"Time = " <<  ((double) (end - start)) << " (a.u.)" << std::endl;  
 
 SaveTree(vars);
 
 std::cerr << " === end === " << std::endl;
}



