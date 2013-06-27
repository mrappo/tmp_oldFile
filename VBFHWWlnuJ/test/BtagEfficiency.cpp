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

#include "METzCalculator.h"
#include "treeReader.h"




/// Main Programme                                                                                                                                                                               
int main (int argc, char** argv){

  if (argc != 2){

    std::cerr << ">>> Usage:   " << argv[1] << "   treeFile.root" << std::endl;
    return -1;
  }


  // parse config file parameter                                                                                                                                                                
  parseConfigFile(argv[1]);

  std::string InputFileName  = gConfigParser -> readStringOption("Input::InputFileName");
  std::string TreeName       = gConfigParser -> readStringOption("Input::TreeName");
  std::string LeptonType     = gConfigParser -> readStringOption("Input::LeptonType");

  std::cout<<"                     "<<std::endl;
  std::cout<<" Input : FileList    "<<InputFileName<<std::endl;
  std::cout<<" Input : TreeName    "<<TreeName<<std::endl;
  std::cout<<" Input : Lepton Type "<<LeptonType<<std::endl;
  std::cout<<"                     "<<std::endl;

  std::string bTagAlgorithm = gConfigParser -> readStringOption("Option::bTagAlgorithm");

  int    JetCollectionSize       = gConfigParser -> readIntOption("Option::JetCollectionSize");

  double JetPtCutMin             = gConfigParser -> readDoubleOption("Option::JetPtCutMin");
  double JetEtaCutMax            = gConfigParser -> readDoubleOption("Option::JetEtaCutMax");
  double bTagWorkingPoint        = gConfigParser -> readDoubleOption("Option::bTagWorkingPoint");

  std::cout<<"                               "<<std::endl;
  std::cout<<" Option : bTag  Algorithm      "<<bTagAlgorithm<<std::endl;
  std::cout<<" Option : Jet Collection Size  "<<JetCollectionSize<<std::endl;
  std::cout<<" Option : Jet Pt Cut           "<<JetPtCutMin<<std::endl;
  std::cout<<" Option : Jet Eta Cut          "<<JetEtaCutMax<<std::endl;
  std::cout<<" Option : B tag Cut            "<<bTagWorkingPoint<<std::endl;
  std::cout<<"                               "<<std::endl;

  int NbinsX        = gConfigParser -> readDoubleOption("Option::NbinsX");
  int NbinsY        = gConfigParser -> readDoubleOption("Option::NbinsY");

  std::vector<double> BinXEdges       = gConfigParser -> readDoubleListOption("Option::BinXEdges");
  std::vector<double> BinYEdges       = gConfigParser -> readDoubleListOption("Option::BinYEdges");


  std::cout<<"                      "<<std::endl;
  std::cout<<" Option : NbinsX      "<<NbinsX<<std::endl;
  std::cout<<" Option : NbinsY      "<<NbinsY<<std::endl;
  std::cout<<"                      "<<std::endl;
  

  std::string OuputFilePath = gConfigParser -> readStringOption("Output::OuputFilePath");
  std::string OuputFileName = gConfigParser -> readStringOption("Output::OuputFileName");

  std::cout<<"                         "<<std::endl;
  std::cout<<" Output : OuputFilePath  "<<OuputFilePath<<std::endl;
  std::cout<<" Outpur : OuputFileName  "<<OuputFileName<<std::endl;
  std::cout<<"                         "<<std::endl;


  // TChain To be created taking all the input file of  InputFileList

  TFile* inputFile = new TFile (InputFileName.c_str(),"READ");

  TTree *inputTreeList = (TTree*) inputFile->Get(TreeName.c_str());
  
  std::cout << " Input Entries  : " << inputTreeList->GetEntries() << " entries in  MC  sample" << std::endl;

  // Read the InpuTree --> set all the branches
  treeReader* fReaderTree     = new treeReader((TTree*)(inputTreeList), false);

  // Create output root file

  TFile* outputFileRoot = new TFile ((OuputFilePath+"/"+OuputFileName+".root").c_str(),"RECREATE");

  double *EdgesX = &BinXEdges.at(0) ;
  double *EdgesY = &BinYEdges.at(0) ;

  TH2F* Numerator_b   = new TH2F("Numerator_b","Numerator_b",NbinsX,EdgesX,NbinsY,EdgesY);
  Numerator_b->Sumw2();
  TH2F* Denominator_b = new TH2F("Denominator_b","Denominator_b",NbinsX,EdgesX,NbinsY,EdgesY);
  Denominator_b->Sumw2();
  TH2F* Efficiency_b  = new TH2F("Efficiency_b","Efficiency_b",NbinsX,EdgesX,NbinsY,EdgesY);
  Efficiency_b->Sumw2();

  TH2F* Numerator_c   = new TH2F("Numerator_c","Numerator_c",NbinsX,EdgesX,NbinsY,EdgesY);
  Numerator_c->Sumw2();
  TH2F* Denominator_c = new TH2F("Denominator_c","Denominator_c",NbinsX,EdgesX,NbinsY,EdgesY);
  Denominator_c->Sumw2();
  TH2F* Efficiency_c  = new TH2F("Efficiency_c","Efficiency_c",NbinsX,EdgesX,NbinsY,EdgesY);
  Efficiency_c->Sumw2();

  TH2F* Numerator_udsg   = new TH2F("Numerator_udsg","Numerator_udsg",NbinsX,EdgesX,NbinsY,EdgesY);
  Numerator_udsg->Sumw2();
  TH2F* Denominator_udsg = new TH2F("Denominator_udsg","Denominator_udsg",NbinsX,EdgesX,NbinsY,EdgesY);
  Denominator_udsg->Sumw2();
  TH2F* Efficiency_udsg  = new TH2F("Efficiency_udsg","Efficiency_udsg",NbinsX,EdgesX,NbinsY,EdgesY);
  Efficiency_udsg->Sumw2();


  std::cout<<"                      "<<std::endl;
  std::cout<<" Start Loop on the Event "<<std::endl;
  std::cout<<"                      "<<std::endl;
  
  for( int iEvent = 0; iEvent < inputTreeList->GetEntries()/10000 ; iEvent ++){

    if(iEvent%10000 ==0) std::cout<<" Event "<<iEvent<<std::endl;
  
    inputTreeList->GetEntry(iEvent);
    
    for( int iJet = 0; iJet < JetCollectionSize ; iJet++){


      if(fReaderTree->getFloat("JetPFCor_Pt")[iJet] < JetPtCutMin ) continue ;
      if(fabs(fReaderTree->getFloat("JetPFCor_Eta")[iJet]) > JetEtaCutMax ) continue ;
      //      std::cout<<" Pt "<<fReaderTree->getFloat("JetPFCor_Pt")[iJet]<<" eta "<<fabs(fReaderTree->getFloat("JetPFCor_Eta")[iJet])<<std::endl;

      // Fill Numerator and denominator for udsg
      //std::cout<<" Flavor "<<fabs(fReaderTree->getInt("JetPFCor_partonFlavour")[iJet])<<std::endl;
          
      if(fabs(fReaderTree->getInt("JetPFCor_partonFlavour")[iJet]) == 1 || fabs(fReaderTree->getInt("JetPFCor_partonFlavour")[iJet]) == 2 ||
         fabs(fReaderTree->getInt("JetPFCor_partonFlavour")[iJet]) == 3 || fabs(fReaderTree->getInt("JetPFCor_partonFlavour")[iJet]) == 21){
         
	if(bTagAlgorithm == "CSV") {
          Denominator_udsg->Fill(fReaderTree->getInt("JetPFCor_Pt")[iJet],fabs(fReaderTree->getInt("JetPFCor_Eta")[iJet]));
	  if( fReaderTree->getFloat("JetPFCor_bDiscriminatorCSV")[iJet] >  bTagWorkingPoint ) 
	    Numerator_udsg -> Fill(fReaderTree->getInt("JetPFCor_Pt")[iJet],fabs(fReaderTree->getInt("JetPFCor_Eta")[iJet]));

	}

         
	if(bTagAlgorithm == "SSV") {
          Denominator_udsg->Fill(fReaderTree->getInt("JetPFCor_Pt")[iJet],fabs(fReaderTree->getInt("JetPFCor_Eta")[iJet]));
	  if( fReaderTree->getFloat("JetPFCor_bDiscriminator")[iJet] >  bTagWorkingPoint ) 
	    Numerator_udsg -> Fill(fReaderTree->getInt("JetPFCor_Pt")[iJet],fabs(fReaderTree->getInt("JetPFCor_Eta")[iJet]));

	}
	    
    }

      // Fill Numerator and denominator for c
      
    else if(fabs(fReaderTree->getInt("JetPFCor_partonFlavour")[iJet]) == 4){

         
	if(bTagAlgorithm == "CSV") {
          Denominator_c->Fill(fReaderTree->getInt("JetPFCor_Pt")[iJet],fabs(fReaderTree->getInt("JetPFCor_Eta")[iJet]));
	  if( fReaderTree->getFloat("JetPFCor_bDiscriminatorCSV")[iJet] >  bTagWorkingPoint ) 
	    Numerator_c -> Fill(fReaderTree->getInt("JetPFCor_Pt")[iJet],fabs(fReaderTree->getInt("JetPFCor_Eta")[iJet]));

	}

         
	if(bTagAlgorithm == "SSV") {
          Denominator_c->Fill(fReaderTree->getInt("JetPFCor_Pt")[iJet],fabs(fReaderTree->getInt("JetPFCor_Eta")[iJet]));
	  if( fReaderTree->getFloat("JetPFCor_bDiscriminator")[iJet] >  bTagWorkingPoint ) 
	    Numerator_c -> Fill(fReaderTree->getInt("JetPFCor_Pt")[iJet],fabs(fReaderTree->getInt("JetPFCor_Eta")[iJet]));

	}

    }

   else if(fabs(fReaderTree->getInt("JetPFCor_partonFlavour")[iJet]) == 5){
         
	if(bTagAlgorithm == "CSV") {
          Denominator_b->Fill(fReaderTree->getInt("JetPFCor_Pt")[iJet],fabs(fReaderTree->getInt("JetPFCor_Eta")[iJet]));
	  if( fReaderTree->getFloat("JetPFCor_bDiscriminatorCSV")[iJet] >  bTagWorkingPoint ) 
	    Numerator_b -> Fill(fReaderTree->getInt("JetPFCor_Pt")[iJet],fabs(fReaderTree->getInt("JetPFCor_Eta")[iJet]));

	}

         
	if(bTagAlgorithm == "SSV") {
          Denominator_b->Fill(fReaderTree->getInt("JetPFCor_Pt")[iJet],fabs(fReaderTree->getInt("JetPFCor_Eta")[iJet]));
	  if( fReaderTree->getFloat("JetPFCor_bDiscriminator")[iJet] >  bTagWorkingPoint ) 
	    Numerator_b -> Fill(fReaderTree->getInt("JetPFCor_Pt")[iJet],fabs(fReaderTree->getInt("JetPFCor_Eta")[iJet]));

	}

    }
      
   }
    
  }
  
  // calculate the efficiency with binomial errors
  Efficiency_b->Divide(Numerator_b,Denominator_b,1.,1.,"B");
  Efficiency_c->Divide(Numerator_c,Denominator_c,1.,1.,"B");
  Efficiency_udsg->Divide(Numerator_udsg,Denominator_udsg,1.,1.,"B");

  outputFileRoot->cd();


  Numerator_b->Write("Numerator_b");
  Denominator_b->Write("Denomiantor_b");
  Efficiency_b->Write("Efficiency_b");

  Numerator_c->Write("Numerator_c");
  Denominator_c->Write("Denominator_c");
  Efficiency_c->Write("Efficiency_c");

  Numerator_udsg->Write("Numerator_udsg");
  Denominator_udsg->Write("Denominator_udsg");
  Efficiency_udsg->Write("Efficiency_udsg");

  outputFileRoot->Close();

  
  return 0 ;

}

