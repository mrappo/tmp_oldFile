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

  std::string nameWJetsPT100  = "/data2/rgerosa/otrees_ExoLeptonID_v2/trainingtrees_el/ofile_WJets_Pythia.root";
  std::string nameWJetsPT180  = "/data2/rgerosa/otrees_ExoLeptonID_v2/trainingtrees_el/ofile_WJets_Pythia180.root";
  std::string nameWJetsMerge1 = "/data2/rgerosa/otrees_ExoLeptonID_v2/trainingtrees_el/ofile_WJets_PythiaMerged1.root";
  std::string nameWJetsMerge2 = "/data2/rgerosa/otrees_ExoLeptonID_v2/trainingtrees_el/ofile_WJets_PythiaMerged2.root";
  std::string nameWJetsMerge  = "/data2/rgerosa/otrees_ExoLeptonID_v2/trainingtrees_el/ofile_WJets_PythiaMerged.root";

  double WJets100WeightMuon = 34.29*(8955318/289.2)*(1/(9739464+(34.29/289.2)*8955318)) ;
  double WJets180WeightMuon = 9739464/(9739464+(34.29/289.2)*8955318) ;


  TFile* inputFileWJetsPT100 = new TFile(nameWJetsPT100.c_str(),"READ");

  TTree* inputTreeWJetsPT100 = (TTree*) inputFileWJetsPT100->Get("otree");

  treeReader* fReaderWJetsPT100 = new treeReader((TTree*)(inputTreeWJetsPT100), false);

  TFile* ouputFile1 = new TFile(nameWJetsMerge1.c_str(),"RECREATE");
  ouputFile1->cd();

  inputTreeWJetsPT100->SetBranchStatus("wSampleWeight",0);

  TTree* outputTree1 = inputTreeWJetsPT100->CloneTree(0);

  float wSampleWeight = 0.;  
  TBranch* branch_wSampleWeight = outputTree1->Branch("wSampleWeight",&wSampleWeight,"wSampleWeight/F");

  inputTreeWJetsPT100->SetBranchStatus("wSampleWeight",1);
 
	
  for( int iEntry = 0 ; iEntry < inputTreeWJetsPT100->GetEntries() ; iEntry++){

    if(iEntry%100000 == 0) std::cout<<" Wjets Pt 100 -> entry = "<<iEntry<<std::endl;

    inputTreeWJetsPT100->GetEntry(iEntry);
    if(fReaderWJetsPT100->getFloat("W_pt_gen")[0] <= 180.){ outputTree1->Fill();
 	                                                    wSampleWeight=fReaderWJetsPT100->getFloat("wSampleWeight")[0];
                                                            branch_wSampleWeight->Fill();
    }
    else{
         
	 wSampleWeight=WJets100WeightMuon*fReaderWJetsPT100->getFloat("wSampleWeight")[0];
         outputTree1->Fill();
         branch_wSampleWeight->Fill();
	 //  std::cout<<" reader branch new "<<wSampleWeight<<" old "<<fReaderWJetsPT100->getFloat("wSampleWeight")[0]<<std::endl;
       }
    
  }

  outputTree1->Write("otree");
  ouputFile1->Close();

  
  TFile* inputFileWJetsPT180 = new TFile(nameWJetsPT180.c_str(),"READ");

  TTree* inputTreeWJetsPT180 = (TTree*) inputFileWJetsPT180->Get("otree");

  treeReader* fReaderWJetsPT180 = new treeReader((TTree*)(inputTreeWJetsPT180), false);

  TFile* ouputFile2 = new TFile(nameWJetsMerge2.c_str(),"RECREATE");
  ouputFile2->cd();

  inputTreeWJetsPT180->SetBranchStatus("wSampleWeight",0);

  TTree* outputTree2 = inputTreeWJetsPT180->CloneTree(0);

  TBranch* branch_wSampleWeight2 = outputTree2->Branch("wSampleWeight",&wSampleWeight,"wSampleWeight/F");

  inputTreeWJetsPT180->SetBranchStatus("wSampleWeight",1);

  for( int iEntry = 0 ; iEntry < inputTreeWJetsPT180->GetEntries() ; iEntry++){

    if(iEntry%100000 == 0) std::cout<<" Wjets Pt 180 -> entry = "<<iEntry<<std::endl;

    inputTreeWJetsPT180->GetEntry(iEntry);

    wSampleWeight=WJets180WeightMuon*fReaderWJetsPT180->getFloat("wSampleWeight")[0];
    outputTree2->Fill();
    branch_wSampleWeight2->Fill();
    //   std::cout<<" reader branch new "<<wSampleWeight2<<" old "<<fReaderWJetsPT180->getFloat("wSampleWeight")[0]<<std::endl;
  }
  
  outputTree2->Write("otree");
  ouputFile2->Close();
    
  system((std::string("hadd -f out.root "+nameWJetsMerge1+" "+nameWJetsMerge2)).c_str());
  system((std::string("mv out.root "+nameWJetsMerge)).c_str());
  system(std::string("rm "+nameWJetsMerge1).c_str());
  system(std::string("rm "+nameWJetsMerge2).c_str());

  return 0;

}

