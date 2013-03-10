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

#include "ntpleUtils.h"
#include "ConfigParser.h"

#include "ReadInputFile.h"

#include "TrainingMVAClass.h"

/// Main Programme                                                                                                                                                                                
int main (int argc, char** argv){

  if (argc != 2){

    std::cerr << ">>> Usage:   " << argv[1] << "   cfg file" << std::endl;
    return -1;
  }

  // Load TTree Lybrary                                                                                                                                                                           
  gSystem->Load("libTree.so");

  // parse config file parameter                                                                                                                                                                  
  parseConfigFile(argv[1]);

  std::string InputDirectory        = gConfigParser -> readStringOption("Input::InputDirectory");
  std::string InputSampleList       = gConfigParser -> readStringOption("Input::InputSampleList");
  std::string InputVariableList     = gConfigParser -> readStringOption("Input::InputVaribleList");
  std::string InputSpectatorList    = gConfigParser -> readStringOption("Input::InputSpectatorList");

  std::string TreeName ; 
  try{ TreeName = gConfigParser -> readStringOption("Input::TreeName"); }
  catch(char const* exceptionString){ 
    TreeName = "WJet";
    std::cerr<<" Tree Name set to --> WJet default "<<std::endl;
  }


  std::string Label         = gConfigParser -> readStringOption("Input::Label");

  std::cout<<"                      "<<std::endl;
  std::cout<<" Input Directory      "<<InputDirectory<<std::endl;
  std::cout<<" Input Sample List    "<<InputSampleList<<std::endl;
  std::cout<<" Input Variable  List "<<InputVariableList<<std::endl;
  std::cout<<" Input Spectator List "<<InputSpectatorList<<std::endl;
  std::cout<<" Input TreeName       "<<TreeName<<std::endl;
  std::cout<<" Input Label          "<<Label<<std::endl;
  std::cout<<"                      "<<std::endl;


  std::string SignalName ;
  try{ SignalName = gConfigParser -> readStringOption("Option::SignalName"); }
  catch(char const* exceptionString){
    SignalName = "qqH";
    std::cerr<<" Signal Name set to --> qqH default "<<std::endl;
  }

  std::cout<<" Option Signal Name "<<SignalName<<std::endl;
  std::cout<<std::endl;

  std::string EventWeight ;

  try{ EventWeight = gConfigParser -> readStringOption("Option::EventWeight"); }
  catch(char const* exceptionString){
    EventWeight = "puwt*effwt";
    std::cerr<<" Event Weight set to --> puwt*effwt default "<<std::endl;
  }

  std::cout<<" Option Event Weight "<<EventWeight<<std::endl;
  std::cout<<std::endl;
 
  std::string PreselectionCut ;
  try{ PreselectionCut  = gConfigParser -> readStringOption("Option::PreselectionCut"); }
  catch(char const* exceptionString){
    PreselectionCut = "";
    std::cerr<<" No preselection cut is set --> all the events are used "<<std::endl;
  }


  std::cout<<" Option Preselection Cut "<<PreselectionCut<<std::endl;
  std::cout<<"                         "<<std::endl;


  std::vector<std::string> UseMethodName;
  try{ UseMethodName = gConfigParser -> readStringListOption("Option::UseMethodName"); }
  catch(char const* exceptionString){ UseMethodName.push_back("Cuts");
                                      UseMethodName.push_back("BDT");
				      std::cerr<<" Default Method are Cuts and BDT --> Set By Default  "<<std::endl;
  }

  std::cout << std::endl;
  std::cout << " >>>>> Option::UseMethodName size = " << UseMethodName.size() << std::endl;
  std::cout << " >>>>> >>>>>  ";
  for (unsigned int iCat = 0; iCat < UseMethodName.size(); iCat++){
    std::cout << " " << UseMethodName.at(iCat) << ", ";
  }
  std::cout << std::endl;


  std::string outputFileDirectory = gConfigParser -> readStringOption("Output::outputFileDirectory"); 
  std::string outputFileName      = gConfigParser -> readStringOption("Output::outputFileName"); 

  std::cout<<"                      "<<std::endl;
  std::cout<<" Outout Directory     "<<outputFileDirectory<<std::endl;
  std::cout<<" Input Sample List    "<<outputFileName<<std::endl;
  std::cout<<"                      "<<std::endl;

  // read sample input file list to Plot                                                                                                                                                          

  std::vector <std::string> NameSample;
  std::vector <std::string> NameReducedSample;
  std::vector <int> ColorSample;
  std::vector <double> SampleCrossSection;
  std::vector <int> NumEntriesBefore;


  if(ReadInputSampleFile(InputSampleList,NameSample,NameReducedSample,ColorSample,SampleCrossSection,NumEntriesBefore) <= 0){
    std::cerr<<" Empty Input Sample File or not Exisisting --> Exit "<<std::endl; return -1;}


  // Read Input File Variables for Training 
  
  std::vector<std::string> mapTrainingVariables;
  
  if(ReadInputVariableFile(InputVariableList,mapTrainingVariables) <= 0){
    std::cerr<<" Empty Input Variable List File or not Exisisting --> Exit "<<std::endl; return -1;}

  // Read Spectator Variables for Training 
 
  std::vector<std::string> mapSpectatorVariables;
  
  if(ReadInputVariableFile(InputSpectatorList,mapSpectatorVariables) <= 0){
    std::cerr<<" Empty Spectator Variable List File or not Exisisting --> Exit "<<std::endl; return -1;}

  // Import Sample and signal - background collections

  std::vector <TFile*> signalFileList;
  std::vector <TFile*> backgroundFileList;
  std::vector <TTree*> signalTreeList;
  std::vector <TTree*> backgroundTreeList;


  for (size_t iSample=0; iSample<NameSample.size(); iSample++){

	TString NameFile = Form("%s/%s.root",InputDirectory.c_str(),NameSample.at(iSample).c_str());
	std::cout<<" Input File : "<< NameFile.Data()<<std::endl;

	if(NameSample.at(iSample) == SignalName ){
          signalFileList.push_back ( new TFile (NameFile.Data(),"READ") );
    	  if(signalFileList.back()!=0) signalTreeList.push_back( (TTree*) signalFileList.at(iSample)->Get(TreeName.c_str()));
        }
        else {
               backgroundFileList.push_back ( new TFile (NameFile.Data(),"READ") );
               if(backgroundFileList.back()!=0) backgroundTreeList.push_back( (TTree*) backgroundFileList.at(iSample)->Get(TreeName.c_str()));
	}
 
  }

  // Book MVA Training Object

  TrainingMVAClass* WWTraining = new TrainingMVAClass(signalTreeList, backgroundTreeList, TreeName, outputFileDirectory, outputFileName, Label);

  // Set Global Weight and signal + background Tree for MVA Training

  std::vector<double> signalGlobalWeight;
  std::vector<double> backgroundGlobalWeight;

  int isSignal = 0;
  int isBackground = 0;
 
  for(size_t iSample =0; iSample<NameSample.size() ; iSample++){

    if( NameSample.at(iSample) == SignalName ) {
       
      if(signalTreeList.size() ==1 ) { signalGlobalWeight.at(isSignal) = 1.0 ; isSignal ++ ; }
      else{
           signalGlobalWeight.at(isSignal) = SampleCrossSection.at(iSample)/NumEntriesBefore.at(iSample);
           isSignal ++;
      }

    }
    else{
         if(backgroundTreeList.size() ==1 ) { backgroundGlobalWeight.at(isBackground) = 1.0 ; isBackground ++ ; }
         else{
               backgroundGlobalWeight.at(isBackground) = SampleCrossSection.at(iSample)/NumEntriesBefore.at(iSample);
               isBackground ++;    
	 }
    }
  }


  WWTraining->BookMVATrees(signalGlobalWeight, backgroundGlobalWeight);

  // Set Input and Spectator Variables

  WWTraining->AddTrainingVariables(mapTrainingVariables, mapSpectatorVariables, EventWeight);

  // Prepare and Set the MVA Factory

  WWTraining->AddPrepareTraining ( PreselectionCut ) ;

  // Book and Run TMVA Training and testing for the selected methods

  for(size_t iMethod =0; iMethod<UseMethodName.size(); iMethod++){

    // Rectangular Cuts
    if(UseMethodName.at(iMethod) == "CutsMC" ) WWTraining->BookandTrainRectangularCuts("MC");
    if(UseMethodName.at(iMethod) == "CutsGA" ) WWTraining->BookandTrainRectangularCuts("GA");
    if(UseMethodName.at(iMethod) == "CutsSA" ) WWTraining->BookandTrainRectangularCuts("SA");
 
    // Likelihood 
    if(UseMethodName.at(iMethod) == "Likelihood")    WWTraining->BookandTrainLikelihood(); 
    if(UseMethodName.at(iMethod) == "LikelihoodKDE") WWTraining->BookandTrainLikelihood("LikelihoodKDE"); 
    if(UseMethodName.at(iMethod) == "PDERS")         WWTraining->BookandTrainLikelihood("PDERS"); 

    // Fisher Discriminant
    if(UseMethodName.at(iMethod) == "Fisher") WWTraining->BookandTrainFisherDiscriminant(); 

    // Linear Discriminant
    if(UseMethodName.at(iMethod) == "LD")     WWTraining->BookandTrainLinearDiscriminant();

    // MLP
    if(UseMethodName.at(iMethod) == "MLP")     WWTraining->BookandTrainMLP();

    // BDT
    if(UseMethodName.at(iMethod) == "BDT")     WWTraining->BookandTrainBDT();

    // BDTG
    if(UseMethodName.at(iMethod) == "BDTG")     WWTraining->BookandTrainBDTG();

    // BDTF
    if(UseMethodName.at(iMethod) == "BDTF")     WWTraining->BookandTrainBDTF();

  }

  // Print Output Plots

  WWTraining->PrintTrainingResults ();
    

  return 0 ;


}
