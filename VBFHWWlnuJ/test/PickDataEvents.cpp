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

#include "treeReader.h"

class EventParameters{

 public:

  EventParameters(const int & eventId, const int & runId, const int & lumiId, const float & l_pt, const float & v_pt, const float & ungroomed_jet_pt, const float & met,
                  const float & mJ, const float & jet_tau2tau1, const float jet_eta, const float & mWW, const float & jet_pt_1, const float & jet_pt_2, const float & deta_jj, const float & mjj){

    eventId_ = eventId;
    runId_   = runId;
    lumiId_  = lumiId; 
    l_pt_    = l_pt;
    v_pt_    = v_pt;
    ungroomed_jet_pt_ = ungroomed_jet_pt; 
    met_     = met;
    mJ_      = mJ;
    tau2tau1_= jet_tau2tau1;
    jet_eta_ = jet_eta ;
    mWW_     = mWW;
    jet_pt_1_ = jet_pt_1;
    jet_pt_2_ = jet_pt_2;
    deta_jj_  = deta_jj;
    mjj_      = mjj;
  }

  bool operator < (const EventParameters & eventB) const {

    return (*this).mWW_ < eventB.mWW_ ;

  }
 
  int eventId_, runId_, lumiId_ ; 
  float l_pt_ , v_pt_, ungroomed_jet_pt_, met_, mJ_, mWW_, jet_pt_1_, jet_pt_2_, deta_jj_, mjj_, tau2tau1_, jet_eta_ ;
};

int main (int arcg, char** argv){

  std::string nameFileDataMuon     = "/afs/cern.ch/user/r/rgerosa/work/ElectroWeakPlusJetFramework/SQWatFramework/CMSSW_6_1_1/src/boostedWWAnalysis/trainingtrees_mu/ofile_data.root" ;
  std::string nameFileDataElectron = "/afs/cern.ch/user/r/rgerosa/work/ElectroWeakPlusJetFramework/SQWatFramework/CMSSW_6_1_1/src/boostedWWAnalysis/trainingtrees_el/ofile_data.root" ;


  TFile* inputFileMuon     = new TFile(nameFileDataMuon.c_str(),"READ");
  TFile* inputFileElectron = new TFile(nameFileDataElectron.c_str(),"READ");

  TTree* inputTreeMuon     = (TTree*) inputFileMuon->Get("otree");
  TTree* inputTreeElectron = (TTree*) inputFileElectron->Get("otree");

  treeReader* fReaderMuon     = new treeReader((TTree*)(inputTreeMuon), false);
  treeReader* fReaderElectron = new treeReader((TTree*)(inputTreeElectron), false);

  std::vector<EventParameters> SelectedEventList;
  
  std::ofstream outFile("output/PickEventList.txt",std::ios::out);

  outFile<<" Selection:"<<std::endl;
  outFile<<" standard lvJ selectio Wmass cut in the signal region [65,105] + VBF selections applied in the analysis"<<std::endl;
  outFile<<" Pz neutrino type-0 "<<std::endl;
  outFile<<" mWW calculated using MET "<<std::endl;
  outFile<<std::endl;
  outFile<<std::endl;

  TH1F* data_mWW    = new TH1F("data_mWW","",19,550,1500);
  TH1F* data_jet_pt = new TH1F("data_jet_pt","",100,200,800);
  TH1F* data_met    = new TH1F("data_met","",100,50,500);
  TH1F* data_mJJ    = new TH1F("data_mJJ","",100,250,800);
  TH1F* data_detaJJ = new TH1F("data_detaJJ","",100,3,9);
  TH1F* data_hadronicTop = new TH1F("data_hadronicTop","",100,200,500);
  TH1F* data_leptonicTop = new TH1F("data_leptonicTop","",100,200,500);


  outFile<<"##############"<<std::endl;
  outFile<<" MUON CHANNEL "<<std::endl;
  outFile<<"##############"<<std::endl;

  outFile<<std::endl;

  std::cout<<" muon data -> entries "<<inputTreeMuon->GetEntries()<<std::endl;
  for( int iEntry = 0 ; iEntry < inputTreeMuon->GetEntries() ; iEntry++){

    inputTreeMuon->GetEntry(iEntry);
    if(iEntry%10000 == 0) std::cout<<" reading muon data: "<<iEntry<<std::endl; 
    if(fReaderMuon->getInt("issignal")[0] == 1 && fReaderMuon->getFloat("v_pt")[0] > 200 && fReaderMuon->getFloat("pfMET")[0] > 50 && fReaderMuon->getFloat("l_pt")[0]>30 && 
       fReaderMuon->getFloat("ungroomed_jet_pt")[0]>200 && (fReaderMuon->getFloat("jet_mass_pr")[0] > 65 && fReaderMuon->getFloat("jet_mass_pr")[0] < 105) && 
       fReaderMuon->getFloat("jet_tau2tau1")[0] < 0.5 && fReaderMuon->getFloat("vbf_maxpt_j1_bDiscriminatorCSV")[0] < 0.679 && 
       fReaderMuon->getFloat("vbf_maxpt_j2_bDiscriminatorCSV")[0] < 0.679 && fReaderMuon->getFloat("mass_ungroomedjet_closerjet")[0] > 200 && 
       fReaderMuon->getFloat("mass_leptonic_closerjet")[0] > 200 && fReaderMuon->getFloat("mass_lvj_type0_met")[0] > 650 && fReaderMuon->getFloat("mass_lvj_type0_met")[0] < 950 &&
       fReaderMuon->getInt("numberJetBin")[0] >= 2 && fReaderMuon->getFloat("vbf_maxpt_jj_m")[0] > 250 && fabs(fReaderMuon->getFloat("vbf_maxpt_j1_eta")[0]-fReaderMuon->getFloat("vbf_maxpt_j2_eta")[0]) > 3.0){
        EventParameters tempEvent(fReaderMuon->getInt("event")[0],fReaderMuon->getInt("event_runNo")[0],fReaderMuon->getInt("event_lumi")[0],fReaderMuon->getFloat("l_pt")[0],fReaderMuon->getFloat("v_pt")[0], fReaderMuon->getFloat("ungroomed_jet_pt")[0], fReaderMuon->getFloat("pfMET")[0],fReaderMuon->getFloat("jet_mass_pr")[0],fReaderMuon->getFloat("jet_tau2tau1")[0],fReaderMuon->getFloat("ungroomed_jet_eta")[0],fReaderMuon->getFloat("mass_lvj_type0_met")[0],fReaderMuon->getFloat("vbf_maxpt_j1_pt")[0],fReaderMuon->getFloat("vbf_maxpt_j2_pt")[0],fabs(fReaderMuon->getFloat("vbf_maxpt_j1_eta")[0]-fReaderMuon->getFloat("vbf_maxpt_j2_eta")[0]),fReaderMuon->getFloat("vbf_maxpt_jj_m")[0]);
        SelectedEventList.push_back(tempEvent);

        data_mWW->Fill(tempEvent.mWW_);
        data_jet_pt->Fill(tempEvent.ungroomed_jet_pt_);
        data_met->Fill(tempEvent.met_);
        data_mJJ->Fill(tempEvent.mjj_);
        data_detaJJ->Fill(tempEvent.deta_jj_);
        data_hadronicTop->Fill(fReaderMuon->getFloat("mass_ungroomedjet_closerjet")[0]);
        data_leptonicTop->Fill(fReaderMuon->getFloat("mass_leptonic_closerjet")[0]);
       
    }
  }

  std::sort(SelectedEventList.begin(),SelectedEventList.end());
  for( unsigned int iEvent = 0 ; iEvent < SelectedEventList.size() ; iEvent++){
       
    outFile<<SelectedEventList.at(iEvent).runId_<<":"<<SelectedEventList.at(iEvent).lumiId_<<":"<<SelectedEventList.at(iEvent).eventId_<<std::endl;
    outFile<<" l_pt = "<<SelectedEventList.at(iEvent).l_pt_<<" * v_pt = "<<SelectedEventList.at(iEvent).v_pt_<<" *  W pt  = "<<SelectedEventList.at(iEvent).ungroomed_jet_pt_<<
             " * met =  "<< SelectedEventList.at(iEvent).met_<<" * W mass = "<<SelectedEventList.at(iEvent).mJ_<<" * tau2tau1 = "<<SelectedEventList.at(iEvent).tau2tau1_<<
             " * W eta = "<<SelectedEventList.at(iEvent).jet_eta_<<" * mWW = "<<SelectedEventList.at(iEvent).mWW_<<" * jet pt 1 = "<<SelectedEventList.at(iEvent).jet_pt_1_<<
             " * jet pt 2 = "<<SelectedEventList.at(iEvent).jet_pt_2_<<" * dEta jj = "<<SelectedEventList.at(iEvent).deta_jj_<<" * mJJ = "<<SelectedEventList.at(iEvent).mjj_<<std::endl;
    outFile<<std::endl;

  }

  outFile<<std::endl;
  outFile<<std::endl;
  outFile<<std::endl;
  outFile<<std::endl;

  outFile<<"##################"<<std::endl;
  outFile<<" ELECTRON CHANNEL "<<std::endl;
  outFile<<"##################"<<std::endl;

  outFile<<std::endl;
  SelectedEventList.clear();

  std::cout<<" electron data -> entries "<<inputTreeElectron->GetEntries()<<std::endl;
  for( int iEntry = 0 ; iEntry < inputTreeElectron->GetEntries() ; iEntry++){

    inputTreeElectron->GetEntry(iEntry);
    if(iEntry%10000 == 0) std::cout<<" reading electron data: "<<iEntry<<std::endl; 

    if(fReaderElectron->getInt("issignal")[0] == 1 && fReaderElectron->getFloat("v_pt")[0] > 200 && fReaderElectron->getFloat("pfMET")[0] > 50 && 
       fReaderElectron->getFloat("l_pt")[0]>30 && 
       fReaderElectron->getFloat("ungroomed_jet_pt")[0]>200 && (fReaderElectron->getFloat("jet_mass_pr")[0] > 65 && fReaderElectron->getFloat("jet_mass_pr")[0] < 105) && 
       fReaderElectron->getFloat("jet_tau2tau1")[0] < 0.5 && fReaderElectron->getFloat("vbf_maxpt_j1_bDiscriminatorCSV")[0] < 0.679 && 
       fReaderElectron->getFloat("vbf_maxpt_j2_bDiscriminatorCSV")[0] < 0.679 && fReaderElectron->getFloat("mass_ungroomedjet_closerjet")[0] > 200 && 
       fReaderElectron->getFloat("mass_leptonic_closerjet")[0] > 200 && fReaderElectron->getFloat("mass_lvj_type0_met")[0] > 650 && fReaderElectron->getFloat("mass_lvj_type0_met")[0] < 950 &&
       fReaderElectron->getInt("numberJetBin")[0] >= 2 && fReaderElectron->getFloat("vbf_maxpt_jj_m")[0] > 250 && fabs(fReaderElectron->getFloat("vbf_maxpt_j1_eta")[0]-fReaderElectron->getFloat("vbf_maxpt_j2_eta")[0]) > 3.0){
         
        EventParameters tempEvent(fReaderElectron->getInt("event")[0],fReaderElectron->getInt("event_runNo")[0],fReaderElectron->getInt("event_lumi")[0],fReaderElectron->getFloat("l_pt")[0],fReaderElectron->getFloat("v_pt")[0], fReaderElectron->getFloat("ungroomed_jet_pt")[0], fReaderElectron->getFloat("pfMET")[0],fReaderElectron->getFloat("jet_mass_pr")[0],fReaderElectron->getFloat("jet_tau2tau1")[0],fReaderElectron->getFloat("ungroomed_jet_eta")[0],fReaderElectron->getFloat("mass_lvj_type0_met")[0],fReaderElectron->getFloat("vbf_maxpt_j1_pt")[0],fReaderElectron->getFloat("vbf_maxpt_j2_pt")[0],fabs(fReaderElectron->getFloat("vbf_maxpt_j1_eta")[0]-fReaderElectron->getFloat("vbf_maxpt_j2_eta")[0]),fReaderElectron->getFloat("vbf_maxpt_jj_m")[0]);
        SelectedEventList.push_back(tempEvent);

        data_mWW->Fill(tempEvent.mWW_);
        data_jet_pt->Fill(tempEvent.ungroomed_jet_pt_);
        data_met->Fill(tempEvent.met_);
        data_mJJ->Fill(tempEvent.mjj_);
        data_detaJJ->Fill(tempEvent.deta_jj_);
        data_hadronicTop->Fill(fReaderElectron->getFloat("mass_ungroomedjet_closerjet")[0] );
        data_leptonicTop->Fill(fReaderElectron->getFloat("mass_leptonic_closerjet")[0]);

    }
  }

  std::sort(SelectedEventList.begin(),SelectedEventList.end());
  for( unsigned int iEvent = 0 ; iEvent < SelectedEventList.size() ; iEvent++){
       
    outFile<<SelectedEventList.at(iEvent).runId_<<":"<<SelectedEventList.at(iEvent).lumiId_<<":"<<SelectedEventList.at(iEvent).eventId_<<std::endl;
    outFile<<" l_pt = "<<SelectedEventList.at(iEvent).l_pt_<<" * v_pt = "<<SelectedEventList.at(iEvent).v_pt_<<" * W pt = "<<SelectedEventList.at(iEvent).ungroomed_jet_pt_<<
             " * met = "<<SelectedEventList.at(iEvent).met_<<" * W mass = "<<SelectedEventList.at(iEvent).mJ_<<" * tau2tau1 = "<<SelectedEventList.at(iEvent).tau2tau1_<<
             " * W eta "<<SelectedEventList.at(iEvent).jet_eta_<<" * mWW = "<<SelectedEventList.at(iEvent).mWW_<<" * jet pt 1 = "<<SelectedEventList.at(iEvent).jet_pt_1_<<
             " * jet pt 2 = "<<SelectedEventList.at(iEvent).jet_pt_2_<<" * deta JJ = "<<SelectedEventList.at(iEvent).deta_jj_<<" * mJJ "<<SelectedEventList.at(iEvent).mjj_<<std::endl;
    outFile<<std::endl;

  }

  TFile* output = new TFile("output/outputDataPlots.root","RECREATE");
  output->cd();

  data_mWW->Write();
  data_jet_pt->Write();
  data_met->Write();   
  data_mJJ->Write();   
  data_detaJJ->Write();
  data_hadronicTop->Write();
  data_leptonicTop->Write();

  output->Close();
   
  return 0 ;

}
