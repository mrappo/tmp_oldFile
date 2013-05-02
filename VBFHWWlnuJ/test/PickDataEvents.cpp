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

  std::string nameFileData = "/data2/rgerosa/otrees_ExoLeptonID_v2/trainingtrees_mu/ofile_data.root" ;

  TFile* inputFile = new TFile(nameFileData.c_str(),"READ");

  TTree* inputTree = (TTree*) inputFile->Get("otree");

  treeReader* fReader = new treeReader((TTree*)(inputTree), false);
  
  std::ofstream outFile("EventList.txt",std::ios::out);

  outFile<<" Selection:"<<std::endl;
  outFile<<" standard lvJ selectio + mWW>1800 GeV, no MT cut, no Wmass cut"<<std::endl;
  outFile<<" Pz neutrino type-2 (smallest of the two root)"<<std::endl;
  outFile<<" mWW calculated using corrected neutrino pT "<<std::endl;
  outFile<<std::endl;
  outFile<<std::endl;


  outFile<<"##############"<<std::endl;
  outFile<<" MUON CHANNEL "<<std::endl;
  outFile<<"##############"<<std::endl;

  outFile<<std::endl;


  for( int iEntry = 0 ; iEntry < inputTree->GetEntries() ; iEntry++){

    inputTree->GetEntry(iEntry);

    if(fReader->getInt("issignal")[0] && fReader->getFloat("v_pt")[0] > 200 && fReader->getFloat("pfMET")[0] > 40 &&
       fReader->getFloat("l_pt")[0]>50 && fReader->getFloat("ungroomed_jet_pt")[0]>200 && abs(fReader->getFloat("l_eta")[0]) < 2.1 && 
       ((fReader->getFloat("jet_mass_pr")[0] > 40 && fReader->getFloat("jet_mass_pr")[0] < 65 ) || (fReader->getFloat("jet_mass_pr")[0] > 105 && fReader->getFloat("jet_mass_pr")[0] < 130 ))
       && fReader->getInt("nbjets_csvm_veto")[0] == 0 && fReader->getFloat("mass_lvj_type2")[0] > 1800 ){

       outFile<<fReader->getInt("event_runNo")[0]<<":"<<fReader->getInt("event_lumi")[0]<<":"<<fReader->getInt("event")[0]<<std::endl;

       outFile<<fReader->getInt("event")[0]<<"*"<<fReader->getInt("event_runNo")[0]<<"*"<<fReader->getFloat("mass_lvj_type2")[0]<<"*"<<
	        fReader->getFloat("l_pt")[0]<<"*"<<fReader->getFloat("l_eta")[0]<<"*"<<fReader->getFloat("l_phi")[0]<<"*"<<fReader->getFloat("pfMET")[0]<<"*"<<
                fReader->getFloat("pfMET_Phi")[0]<<"*"<<fReader->getFloat("v_mt")[0]<<"*"<<fReader->getFloat("v_pt")[0]<<"*"<<fReader->getFloat("ungroomed_jet_pt")[0]<<"*"<<
       	        fReader->getFloat("ungroomed_jet_eta")[0]<<"*"<<fReader->getFloat("ungroomed_jet_phi")[0]<<"*"<<fReader->getFloat("jet_mass_pr")[0]<<"*"<<
                fReader->getFloat("jet_tau2tau1")[0]<<std::endl;

       outFile<<std::endl;
    
    }

  }

  return 0 ;

}
