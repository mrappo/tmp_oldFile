#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <istream>
#include <sstream>

#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "THStack.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TMath.h"
#include "TF1.h"
#include "TH2F.h"
#include "TList.h"

#include "ntpleUtils.h"
#include "ConfigParser.h"
#include "ReadInputFile.h"
#include "DataMCPlotTool.h"

#include "TMVAGlob.h"

/// Main programme 
int main (int argc, char **argv){

  if(argc<2){ std::cout<<" Not correct number of input parameter --> Need Just one cfg file exit "<<std::endl; return -1; }

  // Load TTree Lybrary                                                                                                                                                                   
  gSystem->Load("libTree.so");

  // Set Root style from global enviroment path                                                                                                                                           
  std::string ROOTStyle;
  if(getenv ("ROOTStyle")!=NULL){
    ROOTStyle = getenv ("ROOTStyle");
    gROOT->ProcessLine((".x "+ROOTStyle+"/rootLogon.C").c_str());
    gROOT->ProcessLine((".x "+ROOTStyle+"/rootPalette.C").c_str());
    gROOT->ProcessLine((".x "+ROOTStyle+"/rootColors.C").c_str());
    gROOT->ProcessLine((".x "+ROOTStyle+"/setTDRStyle.C").c_str());
  }
  

  gStyle->SetOptStat(0);
  gStyle->SetPadTopMargin(0.09);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetErrorX(0.5);
  
  parseConfigFile(argv[1]);

  std::vector<std::string> InputFileName     = gConfigParser -> readStringListOption("Input::InputFileName"); // get the input file list for the TMVA root file output after training
  std::vector<std::string> InputVariableOrMethodName     = gConfigParser -> readStringListOption("Input::InputVariableOrMethodName"); // get the input file list for the TMVA root file output after training
  
  std::cout<<"      "<<std::endl;
  for(size_t iFile = 0; iFile < InputFileName.size(); iFile++)   
    std::cout<<" InputFileName: iFile "<<iFile<<" Name : "<<InputFileName.at(iFile)<<std::endl;
  std::cout<<"      "<<std::endl;

  std::cout<<"      "<<std::endl;
  for(size_t iName = 0; iName < InputVariableOrMethodName.size(); iName++)   
    std::cout<<" InputMethodName: iName "<<iName<<" Name : "<<InputVariableOrMethodName.at(iName)<<std::endl;
  std::cout<<"      "<<std::endl;

  std::string outputPlotDirectory = gConfigParser -> readStringOption("Output::outputPlotDirectory"); 

  std::cout<<" Outout Directory     "<<outputPlotDirectory<<std::endl;
  std::cout<<"                      "<<std::endl;

  std::vector<double> jetPTBinofTraining = gConfigParser -> readDoubleListOption("Input::jetPTBinofTraining");
  for(size_t iPTbin = 0; iPTbin < jetPTBinofTraining.size(); iPTbin++)   
    std::cout<<" jetPTBinofTraining: iPTbin "<<iPTbin<<" Value : "<<jetPTBinofTraining.at(iPTbin)<<std::endl;
  std::cout<<"      "<<std::endl;

  if(jetPTBinofTraining.size()!=2){std::cerr<<" Plot one PTbin training for each time ---> exit "<<std::endl; return -1 ; }

  std::string command;
  command = "if [ ! -e "+outputPlotDirectory+" ] ; then mkdir "+outputPlotDirectory+" ; fi";
  std::cout<<" command = "<<command<<std::endl;
  std::cout<<"           "<<std::endl;

  system(command.c_str());

  command = "if [ ! -f "+outputPlotDirectory+" ] ; then rm "+outputPlotDirectory+"/* ; fi";
  std::cout<<" command = "<<command<<std::endl;
  std::cout<<"           "<<std::endl;

  system(command.c_str());

  
  // Declare the object for the manipolation of the TMVA ROOT file
  TMVAGlob* TMVATraining = new TMVAGlob();
  TMVATraining->Initialize();
  
  TMVATraining->openFileInput(InputFileName);
  std::vector<TFile*> inputFile = TMVATraining->GetInputFile();

  TMVATraining->SetMethodName(InputVariableOrMethodName);  

  // Loop on the inputFile
  for(size_t iFile = 0;  iFile < inputFile.size() ; iFile ++){ 

   TIter nextKey(inputFile.at(iFile)->GetListOfKeys()); // iterator to the list of keys in the memory map of the file  
   TKey *key = 0 ; // loop over the keys

   while ( (key = (TKey*) nextKey())) {

     TClass* classType = gROOT->GetClass(key->GetClassName()); // take the class type of each key inside the root file to check what is inside
     if (!classType->InheritsFrom("TDirectory")) continue;     // if it don't herit from TDirectory it is neglet
     TDirectory *dir = (TDirectory*)key->ReadObj(); 
     TString path(dir->GetPath());
     if (path.Contains("multicutMVA")){ TMVATraining->plotEfficiency(inputFile.at(iFile),dir,jetPTBinofTraining.at(0),jetPTBinofTraining.at(1)); // call the plot efficiency function     
                                        TMVATraining->PrintImage(dir,outputPlotDirectory);}

   }

   TMVATraining->plotEfficiency(inputFile.at(iFile),gDirectory,jetPTBinofTraining.at(0),jetPTBinofTraining.at(1)); // call the plot efficiency function 

  }

  TMVATraining->PrintImage(gDirectory,outputPlotDirectory);


  return 0 ;
}
