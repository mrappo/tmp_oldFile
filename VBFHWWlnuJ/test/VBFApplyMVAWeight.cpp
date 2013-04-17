#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <functional>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TMath.h"
#include "TCanvas.h"

#include "TSystem.h"
#include "TROOT.h"


#include "ntpleUtils.h"
#include "ConfigParser.h"

#include "ReadInputFile.h"

#include "ApplyMVAWeightClass.h"

/// Main Programme                                                                                                                                                                              
                                                                                                                                                                                                
int main (int argc, char** argv){

  if (argc != 2){

    std::cerr << ">>> Usage:   " << argv[1] << "   cfg file" << std::endl;
    return -1;
  }

  // Load TTree Lybrary

  gSystem->Load("libTree.so");

  TMVA::Tools::Instance();

  // parse config file                                                                                                                                                                                              
  parseConfigFile(argv[1]);

  std::string InputDirectory        = gConfigParser -> readStringOption("Input::InputDirectory");
  std::string InputSampleName       = gConfigParser -> readStringOption("Input::InputSampleName");
  std::string InputVariableList     = gConfigParser -> readStringOption("Input::InputVariableList");
  std::string InputSpectatorList    = gConfigParser -> readStringOption("Input::InputSpectatorList");

  std::string LeptonType    = gConfigParser -> readStringOption("Input::LeptonType");

  std::string InputWeightFilePath  = gConfigParser -> readStringOption("Input::InputWeightFilePath");

  std::string TreeName ;
  try{ TreeName = gConfigParser -> readStringOption("Input::TreeName"); }
  catch(char const* exceptionString){
    TreeName = "WJet";
    std::cerr<<" Tree Name set to --> WJet default "<<std::endl;
  }

  
  std::cout<<std::endl;
  std::cout<<" Input Directory      = "<<InputDirectory<<std::endl;
  std::cout<<" Input Sample Name    = "<<InputSampleName<<std::endl;
  std::cout<<" Input Variable  List = "<<InputVariableList<<std::endl;
  std::cout<<" Input Spectator List = "<<InputSpectatorList<<std::endl;
  std::cout<<" Input TreeName       = "<<TreeName<<std::endl;
  std::cout<<" Input LeptonType     = "<<LeptonType<<std::endl;
  std::cout<<std::endl;

  std::string PreselectionCutType ;
  try{ PreselectionCutType  = gConfigParser -> readStringOption("Option::PreselectionCutType"); }
  catch(char const* exceptionString){
    PreselectionCutType = "none";
    std::cerr<<" No preselection cut is set --> all the events are used "<<std::endl;
  }


  std::cout<<" Option Preselection Cut = "<<PreselectionCutType<<std::endl;
  std::cout<<std::endl;

  std::vector<std::string> MethodLabelName;
  try{ MethodLabelName = gConfigParser -> readStringListOption("Option::MethodLabelName"); }
  catch(char const* exceptionString){ MethodLabelName.push_back("MLP");
    MethodLabelName.push_back("BDT");
    std::cerr<<" Default Method are MLP and BDT --> Set By Default  "<<std::endl;
  }
 
  
  std::vector<std::string> UseMethodName;
  try{ UseMethodName = gConfigParser -> readStringListOption("Option::UseMethodName"); }
  catch(char const* exceptionString){ UseMethodName.push_back("Cuts");
    UseMethodName.push_back("BDT");
    std::cerr<<" Default Method are MLP and BDT --> Set By Default  "<<std::endl;
  }

  std::cout << std::endl;
  std::cout << " >>>>> Option::UseMethodName size = " << UseMethodName.size() << std::endl;

  std::vector<double> JetPtBinOfTraining;
  try{ JetPtBinOfTraining = gConfigParser -> readDoubleListOption("Option::JetPtBinOfTraining"); }
  catch(char const* exceptionString){ JetPtBinOfTraining.push_back(0);
    JetPtBinOfTraining.push_back(2000);
    std::cerr<<" Default Method just 1 bin  "<<std::endl;
  }

  std::cout << std::endl;
  std::cout << " >>>>> Option::JetPtBinOfTraining size = " << JetPtBinOfTraining.size()/2 +1 << std::endl;

  for (unsigned int iCat = 0; iCat+1 < JetPtBinOfTraining.size(); iCat++){
    std::cout << " bin min =  " << JetPtBinOfTraining.at(iCat) << " ;  bin max =  "<<JetPtBinOfTraining.at(iCat+1) <<std::endl;
  }
  std::cout << std::endl;
   
  // Read Input File Variables for Training                                                                                                                                                      
  std::vector<std::string> mapTrainingVariables;

  std::cout<<" Read Training Variables File List "<<std::endl;
  std::cout<<std::endl;

  if(ReadInputVariableFile(InputVariableList,mapTrainingVariables) <= 0){
   std::cerr<<" Empty Input Variable List File or not Exisisting --> Exit "<<std::endl; return -1;}

  // Read Spectator Variables for Training                                                                                                                                                    
  std::vector<std::string> mapSpectatorVariables;

  std::cout<<" Read Spectator File List "<<std::endl;
  std::cout<<std::endl;

  if(ReadInputVariableFile(InputSpectatorList,mapSpectatorVariables) <= 0){
    std::cerr<<" Empty Spectator Variable List File or not Exisisting --> Exit "<<std::endl; return -1;}

  std::vector <TFile*> SampleFileList;
  std::vector <TTree*> SampleTreeList;

  std::cout<<std::endl;

  TString NameFile = Form("%s/%s.root",InputDirectory.c_str(),InputSampleName.c_str());
  std::cout<<" Input File : "<< NameFile.Data()<<std::endl;

  SampleFileList.push_back ( new TFile (NameFile.Data(),"UPDATE") );
  if(SampleFileList.back()!=0) SampleTreeList.push_back( (TTree*) SampleFileList.back()->Get(TreeName.c_str()));

  std::cout<<std::endl;
 
  // Book MVA Reader  Object --> one for each pT bin                                                                                                                                            
  std::vector<ApplyMVAWeightClass*> WWReaderVector ;

  std::string weightFile ;
  
  for(size_t pTBin = 0; pTBin+1 < JetPtBinOfTraining.size() ; pTBin++){

    std::cout<<" pT bin of Training: Min = "<<JetPtBinOfTraining.at(pTBin)<<" Max = "<<JetPtBinOfTraining.at(pTBin+1)<<std::endl;
    std::cout<<std::endl;
    for( size_t iMethod = 0 ; iMethod < UseMethodName.size() && iMethod < MethodLabelName.size() ; iMethod ++){

      TString Label_ = Form("otree_%s_PTBin_%d_%d_%s",MethodLabelName.at(iMethod).c_str(),int(JetPtBinOfTraining.at(pTBin)),
                                                      int(JetPtBinOfTraining.at(pTBin+1)),UseMethodName.at(iMethod).c_str()) ;
      std::string buffer ;
      if(UseMethodName.at(iMethod)!="MLP"){

	system((" ls "+InputWeightFilePath+"/ | grep "+std::string(Label_)+" | grep weights | grep xml > ./tmp.txt ").c_str());
	//	std::cout<<" ls "+InputWeightFilePath+"/ | grep "+std::string(Label_)+" | grep weights | grep xml > ./tmp.txt "<<std::endl;
	std::ifstream inputFile ("./tmp.txt");
	while(!inputFile.eof()){ getline(inputFile,buffer); 
	                         if(buffer.empty() || !buffer.find("#") || buffer==" " ) continue;
                                 weightFile = buffer ;
        }
	
        system("rm ./tmp.txt");

      }
 
      else{
            system((" ls "+InputWeightFilePath+"/ | grep "+std::string(Label_)+" | grep weights | grep xml | grep BFGS > ./tmp.txt").c_str());
	    //      	    std::cout<<" ls "+InputWeightFilePath+"/ | grep "+std::string(Label_)+" | grep weights | grep xml | grep BFGS > ./tmp.txt"<<std::endl;
	    std::ifstream inputFile ("tmp.txt");
   	    while(!inputFile.eof()){ getline(inputFile,buffer); 
	                             if(buffer.empty() || !buffer.find("#") || buffer==" " ) continue;
                                     weightFile = buffer ;
	    }
	    system("rm ./tmp.txt");
      }

      std::cout<<" Weight File Name " <<weightFile<<std::endl;
      std::cout<<std::endl;

      WWReaderVector.push_back(new ApplyMVAWeightClass(SampleTreeList, TreeName,InputWeightFilePath, std::string(Label_)));
           
      WWReaderVector.back()->AddTrainingVariables(mapTrainingVariables, mapSpectatorVariables);
      
      WWReaderVector.back()->AddPrepareReader(LeptonType,PreselectionCutType,&JetPtBinOfTraining,pTBin);

      TString NameBranch = Form("%s_PTBin_%d_%d",UseMethodName.at(iMethod).c_str(),int(JetPtBinOfTraining.at(pTBin)),int(JetPtBinOfTraining.at(pTBin+1))) ;

      WWReaderVector.back()->BookMVAWeight(UseMethodName.at(iMethod),weightFile, std::string(NameBranch));

      WWReaderVector.back()->FillMVAWeight(LeptonType,PreselectionCutType);
      
   }

 }

  //Print Output Plots                                                                                                                                                                          
  std::cout<<std::endl;
  std::cout<<" Save Output Root File  ..  "<<std::endl;
  std::cout<<std::endl;

  return 0 ;


}
