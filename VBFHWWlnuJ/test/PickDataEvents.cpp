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

int main (int arcg, char** argv){

  std::string nameFileDataMuon     = "/data2/rgerosa/OTREES/otrees_ExoLeptonID_Analysis/trainingtrees_mu/ofile_data.root" ;
  std::string nameFileDataElectron = "/data2/rgerosa/OTREES/otrees_ExoLeptonID_Analysis/trainingtrees_el/ofile_data.root" ;


  TFile* inputFileMuon     = new TFile(nameFileDataMuon.c_str(),"READ");
  TFile* inputFileElectron = new TFile(nameFileDataElectron.c_str(),"READ");

  TTree* inputTreeMuon     = (TTree*) inputFileMuon->Get("otree");
  TTree* inputTreeElectron = (TTree*) inputFileElectron->Get("otree");

  treeReader* fReaderMuon     = new treeReader((TTree*)(inputTreeMuon), false);
  treeReader* fReaderElectron = new treeReader((TTree*)(inputTreeElectron), false);
  
  std::ofstream outFile("output/PickEventList.txt",std::ios::out);

  outFile<<" Selection:"<<std::endl;
  outFile<<" standard lvJ selectio + mWW>1700 GeV, no MT cut, Wmass cut in the signal region [65,105]"<<std::endl;
  outFile<<" Pz neutrino type-2 (smallest of the two root)"<<std::endl;
  outFile<<" mWW calculated using MET "<<std::endl;
  outFile<<std::endl;
  outFile<<std::endl;


  outFile<<"##############"<<std::endl;
  outFile<<" MUON CHANNEL "<<std::endl;
  outFile<<"##############"<<std::endl;

  outFile<<std::endl;


  for( int iEntry = 0 ; iEntry < inputTreeMuon->GetEntries() ; iEntry++){

    inputTreeMuon->GetEntry(iEntry);

    if(fReaderMuon->getInt("issignal")[0] && fReaderMuon->getFloat("v_pt")[0] > 200 && fReaderMuon->getFloat("pfMET")[0] > 40 &&
       fReaderMuon->getFloat("l_pt")[0]>50 && fReaderMuon->getFloat("ungroomed_jet_pt")[0]>200 && abs(fReaderMuon->getFloat("l_eta")[0]) < 2.1 && 
       // ((fReaderMuon->getFloat("jet_mass_pr")[0] > 40  && fReaderMuon->getFloat("jet_mass_pr")[0] < 65 ) || 
	//        (fReaderMuon->getFloat("jet_mass_pr")[0] > 105  && fReaderMuon->getFloat("jet_mass_pr")[0] < 130 )) && 
       (fReaderMuon->getFloat("jet_mass_pr")[0] > 65 && fReaderMuon->getFloat("jet_mass_pr")[0] < 105) &&
       fReaderMuon->getFloat("nbjets_csvm_veto")[0] == 0 && fReaderMuon->getFloat("mass_lvj")[0] > 1700 && fabs(fReaderMuon->getFloat("ungroomed_jet_eta")[0])<2.4){

       outFile<<fReaderMuon->getInt("event_runNo")[0]<<":"<<fReaderMuon->getInt("event_lumi")[0]<<":"<<fReaderMuon->getInt("event")[0]<<std::endl;

       outFile<<fReaderMuon->getInt("event")[0]<<"*"<<fReaderMuon->getInt("event_runNo")[0]<<"*"<<fReaderMuon->getFloat("mass_lvj")[0]<<"*"<<
	        fReaderMuon->getFloat("l_pt")[0]<<"*"<<fReaderMuon->getFloat("l_eta")[0]<<"*"<<fReaderMuon->getFloat("l_phi")[0]<<"*"<<fReaderMuon->getFloat("pfMET")[0]<<"*"<<
                fReaderMuon->getFloat("pfMET_Phi")[0]<<"*"<<fReaderMuon->getFloat("v_mt")[0]<<"*"<<fReaderMuon->getFloat("v_pt")[0]<<"*"<<fReaderMuon->getFloat("ungroomed_jet_pt")[0]<<"*"<<
       	        fReaderMuon->getFloat("ungroomed_jet_eta")[0]<<"*"<<fReaderMuon->getFloat("ungroomed_jet_phi")[0]<<"*"<<fReaderMuon->getFloat("jet_mass_pr")[0]<<"*"<<
	 fReaderMuon->getFloat("jet_tau2tau1")[0]<<std::endl;

       outFile<<std::endl;
    
    }

  }

  outFile<<std::endl;
  outFile<<std::endl;
  outFile<<std::endl;
  outFile<<std::endl;

  outFile<<"##############"<<std::endl;
  outFile<<" ELECTRON CHANNEL "<<std::endl;
  outFile<<"##############"<<std::endl;

  outFile<<std::endl;

  for( int iEntry = 0 ; iEntry < inputTreeElectron->GetEntries() ; iEntry++){

    inputTreeElectron->GetEntry(iEntry);

    if(fReaderElectron->getInt("issignal")[0] && fReaderElectron->getFloat("v_pt")[0] > 200 && fReaderElectron->getFloat("pfMET")[0] > 80 &&
       fReaderElectron->getFloat("l_pt")[0]>90 && fReaderElectron->getFloat("ungroomed_jet_pt")[0]>200 && abs(fReaderElectron->getFloat("l_eta")[0]) < 2.4 && 
       //       ((fReaderElectron->getFloat("jet_mass_pr")[0] > 40 && fReaderElectron->getFloat("jet_mass_pr")[0] < 65) || 
       //        (fReaderElectron->getFloat("jet_mass_pr")[0] > 105 && fReaderElectron->getFloat("jet_mass_pr")[0] < 130)) && 
       (fReaderElectron->getFloat("jet_mass_pr")[0] > 65 && fReaderElectron->getFloat("jet_mass_pr")[0] < 105) &&
       fReaderElectron->getFloat("nbjets_csvm_veto")[0] == 0 && fReaderElectron->getFloat("mass_lvj")[0] > 1700 && fabs(fReaderElectron->getFloat("ungroomed_jet_eta")[0])<2.4){
    
       outFile<<fReaderElectron->getInt("event_runNo")[0]<<":"<<fReaderElectron->getInt("event_lumi")[0]<<":"<<fReaderElectron->getInt("event")[0]<<std::endl;

       outFile<<fReaderElectron->getInt("event")[0]<<"*"<<fReaderElectron->getInt("event_runNo")[0]<<"*"<<fReaderElectron->getFloat("mass_lvj")[0]<<"*"<<
	        fReaderElectron->getFloat("l_pt")[0]<<"*"<<fReaderElectron->getFloat("l_eta")[0]<<"*"<<fReaderElectron->getFloat("l_phi")[0]<<"*"<<fReaderElectron->getFloat("pfMET")[0]<<"*"<<
                fReaderElectron->getFloat("pfMET_Phi")[0]<<"*"<<fReaderElectron->getFloat("v_mt")[0]<<"*"<<fReaderElectron->getFloat("v_pt")[0]<<"*"<<
                fReaderElectron->getFloat("ungroomed_jet_pt")[0]<<"*"<<
       	        fReaderElectron->getFloat("ungroomed_jet_eta")[0]<<"*"<<fReaderElectron->getFloat("ungroomed_jet_phi")[0]<<"*"<<fReaderElectron->getFloat("jet_mass_pr")[0]<<"*"<<
                fReaderElectron->getFloat("jet_tau2tau1")[0]<<std::endl;

       outFile<<std::endl;
    
    }

  }


  return 0 ;

}
