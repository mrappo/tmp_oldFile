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

TString GetPreselectionCut (const std::string & LeptonType,const std::string & preselectionCutType, const double & pTJetMin_, const double & pTJetMax_);

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

  std::vector<std::string> InputFileName ;
  try{ InputFileName   = gConfigParser -> readStringListOption("Input::InputFileName"); } // get the input file list for the TMVA root file output after training 
  catch(char* exceptionString){
    std::cerr<<" No input File Name list provided --> exit from the code"<<std::endl;
    return -1 ;
  }

  std::cout<<"      "<<std::endl;
  for(size_t iFile = 0; iFile < InputFileName.size(); iFile++)   
    std::cout<<" InputFileName: iFile "<<iFile<<" Name : "<<InputFileName.at(iFile)<<std::endl;


  std::vector<std::string> InputVariableOrMethodName ;
  try{ InputVariableOrMethodName = gConfigParser -> readStringListOption("Input::InputVariableOrMethodName"); } // get the input file list for the TMVA root file output after training
  catch(char* exceptionString){ std::cout<<"No variable name provided --> assume no rectuagular cut are done"<<std::endl; }

  

  std::cout<<"      "<<std::endl;
  for(size_t iName = 0; iName < InputVariableOrMethodName.size(); iName++)   
    std::cout<<" InputMethodName: iName "<<iName<<" Name : "<<InputVariableOrMethodName.at(iName)<<std::endl;

  std::string InputDirectory;
  try{ InputDirectory  = gConfigParser -> readStringOption("Input::InputDirectory");}
  catch(char* exceptionString) {
    std::cerr<<" No input Directory specified for trees --> exit from the code "<<std::endl;
    return -1 ;
  }

  std::cout<<"      "<<std::endl;
  std::cout<<" InputDirectory "<<InputDirectory<<std::endl;

  std::string InputSampleList ;
  try{ InputSampleList  = gConfigParser -> readStringOption("Input::InputSampleList"); }
  catch(char* exceptionString){
    std::cerr<<" No input SampleList specified for trees --> exit from the code "<<std::endl;
    return -1 ;
  }

  std::cout<<"      "<<std::endl;
  std::cout<<" InputSampleList "<<InputSampleList<<std::endl;

  std::string SignalName ;
  try{ SignalName = gConfigParser -> readStringOption("Option::SignalName"); }
  catch(char* exceptionString){
    std::cerr<<" No input SignalName specified for trees --> exit from the code "<<std::endl;
    return -1 ;
  }

  std::cout<<"      "<<std::endl;
  std::cout<<" SignalName "<<SignalName<<std::endl;
  
  std::string EventWeight ;
  try{  EventWeight = gConfigParser -> readStringOption("Option::EventWeight"); }
  catch(char* exceptionString){
    std::cerr<<" No input EventWeight specified for trees --> exit from the code "<<std::endl;
    return -1 ;
  }

  std::cout<<"      "<<std::endl;
  std::cout<<" EventWeight "<<EventWeight<<std::endl;

  
  std::string TreeName ;
  try{  TreeName = gConfigParser -> readStringOption("Input::TreeName"); }
  catch(char* exceptionString){
    std::cerr<<" No input TreeName specified for trees --> set to otree "<<std::endl;
    TreeName = "otree";
  }

  
  std::cout<<"      "<<std::endl;
  std::cout<<" TreeName "<<TreeName<<std::endl;

  std::vector<double> jetPTBinofTraining ;
  try{ jetPTBinofTraining = gConfigParser -> readDoubleListOption("Option::jetPTBinofTraining");}
  catch(char* exceptionString){
    std::cerr<<" jetPTBinofTraining set to 200, 2000 by default "<<std::endl;
    jetPTBinofTraining.push_back(200);
    jetPTBinofTraining.push_back(2000);
  }

  std::cout<<"      "<<std::endl;
  for(size_t iPTbin = 0; iPTbin < jetPTBinofTraining.size(); iPTbin++)   
    std::cout<<" jetPTBinofTraining: iPTbin "<<iPTbin<<" Value : "<<jetPTBinofTraining.at(iPTbin)<<std::endl;

  if(jetPTBinofTraining.size()!=2){std::cerr<<" Plot one PTbin training for each time ---> exit "<<std::endl; return -1 ; }

  std::string LeptonType ; 
  try{ LeptonType = gConfigParser -> readStringOption("Option::LeptonType"); }
  catch(char* exceptionString){
    std::cerr<<" LeptonType not specified --> set to MuonEle by default"<<std::endl;
    LeptonType = "MuonEle" ;
  }

  std::cout<<"                "<<std::endl;
  std::cout<<" LeptonType     "<<LeptonType<<std::endl;

  std::string PreselectionCutType ;
  try{ PreselectionCutType  = gConfigParser -> readStringOption("Output::PreselectionCutType"); }
  catch(char* exceptionString){
    std::cerr<<" Preselection Cut type applied to the definition of the training region don't specified --> exit from the program"<<std::endl;
    return -1 ;
  }

  std::cout<<"                         "<<std::endl;
  std::cout<<" PreselectionCutType     "<<PreselectionCutType<<std::endl;

  
  std::string outputPlotDirectory ;
  try{ outputPlotDirectory = gConfigParser -> readStringOption("Output::outputPlotDirectory"); }
  catch(char* exceptionString){
    std::cerr<<" Output plot directory don't specified --> set by default to output/TMVATrainingPlots"<<std::endl;
    outputPlotDirectory = "output/TMVATrainingPlots";
  }

  std::cout<<"                      "<<std::endl;
  std::cout<<" Outout Directory     "<<outputPlotDirectory<<std::endl;
  std::cout<<"                      "<<std::endl;

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

  // Loop on the inputFile and do the ROC plot
  for(size_t iFile = 0;  iFile < inputFile.size() ; iFile ++){ 

   TIter nextKey(inputFile.at(iFile)->GetListOfKeys()); // iterator to the list of keys in the memory map of the file  
   TKey *key = 0 ; // loop over the keys

   while ( (key = (TKey*) nextKey())) {

     TClass* classType = gROOT->GetClass(key->GetClassName()); // take the class type of each key inside the root file to check what is inside
     if (!classType->InheritsFrom("TDirectory")) continue;     // if it don't herit from TDirectory it is neglet
     TDirectory *dir = (TDirectory*)key->ReadObj(); 
     TString path(dir->GetPath());
     if (path.Contains("multicutMVA")){ TMVATraining->plotEfficiency(inputFile.at(iFile),dir,jetPTBinofTraining.at(0),jetPTBinofTraining.at(1)); // call the plot efficiency function     
                                        TMVATraining->PrintImageROC(dir,outputPlotDirectory);
     }

   }

   TMVATraining->plotEfficiency(inputFile.at(iFile),gDirectory,jetPTBinofTraining.at(0),jetPTBinofTraining.at(1)); // call the plot efficiency function 
  }

  TMVATraining->PrintImageROC(gDirectory,outputPlotDirectory);

  // Read List of Input Files in order to get the number of expected signal and background events
  double numberSignalEvents = 0;
  double numberBackgroundEvents = 0;

  std::vector <std::string> NameSample;
  std::vector <std::string> NameReducedSample;
  std::vector <int> ColorSample;
  std::vector <double> SampleCrossSection;
  std::vector <int> NumEntriesBefore;

  if(ReadInputSampleFile(InputSampleList,NameSample,NameReducedSample,ColorSample,SampleCrossSection,NumEntriesBefore) <= 0){
    std::cerr<<" Empty Input Sample File or not Exisisting --> Exit "<<std::endl; return -1;}

  TString CutString = GetPreselectionCut(LeptonType,PreselectionCutType,jetPTBinofTraining.at(0),jetPTBinofTraining.at(1));

  std::vector <TTree*> TreeVect;
  std::vector <TFile*> FileVect;

  TH1F* histos[NameSample.size()];
  TString hname ;

  for (size_t iSample=0; iSample<NameSample.size(); iSample++){

    TString NameFile = Form("%s/%s.root",InputDirectory.c_str(),NameSample.at(iSample).c_str());
    std::cout<<" Input File : "<< NameFile.Data()<<std::endl;

    FileVect.push_back ( new TFile (NameFile.Data(),"READ") );
    TreeVect.push_back( (TTree*) FileVect.at(iSample)->Get(TreeName.c_str()));

    hname.Form ("%s_%d",NameSample.at(iSample).c_str(),int(iSample));
    hname.ReplaceAll("[","_");
    hname.ReplaceAll("]","_");
    hname.ReplaceAll("(","_");
    hname.ReplaceAll(")","_");

    histos[iSample] = new TH1F (hname.Data(),"",1000,0,5000);
    histos[iSample]->Sumw2();

    if( NameReducedSample.at(iSample) == "DATA" ){
      TreeVect.at(iSample)-> Draw(("l_pt >> "+std::string(hname)).c_str(), std::string(CutString).c_str() ,"goff");
      std::cout<<" Data Entries "<<histos[iSample]->GetEntries()<< " weigthed events "<<histos[iSample]->Integral(0,1000)<<std::endl;      
    }

    else if(NameReducedSample.at(iSample) == SignalName && SignalName!="NULL"){
      TreeVect.at(iSample)->Draw(("l_pt >> "+std::string(hname)).c_str(),("("+EventWeight+")*( "+std::string(CutString)+")").c_str() ,"goff");
      std::cout<<" Signal ggH Entries "<<histos[iSample]->GetEntries()<< " weigthed events "<<histos[iSample]->Integral(0,1000)<<std::endl;
      numberSignalEvents = numberSignalEvents + histos[iSample]->Integral(0,1000) ;     
    }

    else {

      TreeVect.at(iSample)->Draw(("l_pt >> "+std::string(hname)).c_str(),("("+EventWeight+") * ("+std::string(CutString)+")").c_str() ,"goff");
      std::cout<<" Bkg "<<NameSample.at(iSample)<<" Entries "<<histos[iSample]->GetEntries()<<" weighted events "<<
      histos[iSample]->Integral(0,1000)<<std::endl;
      numberBackgroundEvents = numberBackgroundEvents + histos[iSample]->Integral(0,1000) ;
    }

  }


  // Plot correlation variables
  for(size_t iFile = 0;  iFile < inputFile.size() ; iFile ++){ 

   TIter nextKey(inputFile.at(iFile)->GetListOfKeys()); // iterator to the list of keys in the memory map of the file  
   TKey *key = 0 ; // loop over the keys

   while ( (key = (TKey*) nextKey())) {

     TClass* classType = gROOT->GetClass(key->GetClassName()); // take the class type of each key inside the root file to check what is inside
     if (!classType->InheritsFrom("TDirectory")) continue;     // if it don't herit from TDirectory it is neglet
     TDirectory *dir = (TDirectory*)key->ReadObj(); 
     TString path(dir->GetPath());
     if (path.Contains("multicutMVA")){ TMVATraining->plotCorrelationMatrix(inputFile.at(iFile),iFile,outputPlotDirectory); // call the plot efficiency function                 
                                        TMVATraining->plotMVAs(inputFile.at(iFile),TMVATraining->MVAType,outputPlotDirectory);
                                        TMVATraining->plotMVAs(inputFile.at(iFile),TMVATraining->ProbaType,outputPlotDirectory);
                                        TMVATraining->plotMVAs(inputFile.at(iFile),TMVATraining->CompareType,outputPlotDirectory);

                                TMVATraining->plotSignificance(inputFile.at(iFile),TMVATraining->SoverB, numberSignalEvents, numberBackgroundEvents,true,true,outputPlotDirectory);
                                TMVATraining->plotSignificance(inputFile.at(iFile),TMVATraining->SoverSqrtB, numberSignalEvents, numberBackgroundEvents,true,true,outputPlotDirectory);
                                TMVATraining->plotSignificance(inputFile.at(iFile),TMVATraining->SoverSqrtSB, numberSignalEvents, numberBackgroundEvents,true,true,outputPlotDirectory);
                                TMVATraining->plotSignificance(inputFile.at(iFile),TMVATraining->Pvalue, numberSignalEvents, numberBackgroundEvents,true,true,outputPlotDirectory);

                                TMVATraining->plotSignificance(inputFile.at(iFile),TMVATraining->SoverB, numberSignalEvents, numberBackgroundEvents,false,false,outputPlotDirectory);
                                TMVATraining->plotSignificance(inputFile.at(iFile),TMVATraining->SoverSqrtB, numberSignalEvents, numberBackgroundEvents,false,false,outputPlotDirectory);
                                TMVATraining->plotSignificance(inputFile.at(iFile),TMVATraining->SoverSqrtSB, numberSignalEvents, numberBackgroundEvents,false,false,outputPlotDirectory);
                                TMVATraining->plotSignificance(inputFile.at(iFile),TMVATraining->Pvalue, numberSignalEvents, numberBackgroundEvents,false,false,outputPlotDirectory);
     }

   }

    TMVATraining->plotCorrelationMatrix(inputFile.at(iFile),iFile,outputPlotDirectory); // call the plot efficiency function 
    TMVATraining->plotMVAs(inputFile.at(iFile),TMVATraining->MVAType,outputPlotDirectory);
    TMVATraining->plotMVAs(inputFile.at(iFile),TMVATraining->ProbaType,outputPlotDirectory);
    TMVATraining->plotMVAs(inputFile.at(iFile),TMVATraining->CompareType,outputPlotDirectory);

    TMVATraining->plotSignificance(inputFile.at(iFile),TMVATraining->SoverB, numberSignalEvents, numberBackgroundEvents,true,true,outputPlotDirectory);
    TMVATraining->plotSignificance(inputFile.at(iFile),TMVATraining->SoverSqrtB, numberSignalEvents, numberBackgroundEvents,true,true,outputPlotDirectory);
    TMVATraining->plotSignificance(inputFile.at(iFile),TMVATraining->SoverSqrtSB, numberSignalEvents, numberBackgroundEvents,true,true,outputPlotDirectory);
    TMVATraining->plotSignificance(inputFile.at(iFile),TMVATraining->Pvalue, numberSignalEvents, numberBackgroundEvents,true,true,outputPlotDirectory);

    TMVATraining->plotSignificance(inputFile.at(iFile),TMVATraining->SoverB, numberSignalEvents, numberBackgroundEvents,false,false,outputPlotDirectory);
    TMVATraining->plotSignificance(inputFile.at(iFile),TMVATraining->SoverSqrtB, numberSignalEvents, numberBackgroundEvents,false,false,outputPlotDirectory);
    TMVATraining->plotSignificance(inputFile.at(iFile),TMVATraining->SoverSqrtSB, numberSignalEvents, numberBackgroundEvents,false,false,outputPlotDirectory);
    TMVATraining->plotSignificance(inputFile.at(iFile),TMVATraining->Pvalue, numberSignalEvents, numberBackgroundEvents,false,false,outputPlotDirectory);

   }

  return 0 ;
}


TString GetPreselectionCut (const std::string & LeptonType,const std::string & preselectionCutType, const double & pTJetMin_, const double & pTJetMax_){


  if(preselectionCutType == "basicPreselectionCutEXO" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form("issignal && v_pt > 200 && pfMET > 40 && l_pt > 50 && ungroomed_jet_pt > 200 && nbjets_csvm_veto == 0 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",
                pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicPreselectionCutEXO" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt > 200 && pfMET > 80 && l_pt > 90 && ungroomed_jet_pt > 200 && nbjets_csvm_veto == 0 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",
                pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicPreselectionCutEXO" && (LeptonType == "MuEl" || LeptonType == "muel" || LeptonType == "MuonEle" || LeptonType == "muonele") )
    return Form("issignal && v_pt > 200 && pfMET > 80 && l_pt > 90 && ungroomed_jet_pt > 200 && nbjets_csvm_veto == 0 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",
                pTJetMin_,pTJetMax_);


  else if(preselectionCutType == "basicSBPreselectionCutEXO" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form("issignal && v_pt > 200 && pfMET > 40 && l_pt > 50 && ungroomed_jet_pt > 200 && ( ( jet_mass_pr >=40 && jet_mass_pr <= 65 ) || ( jet_mass_pr >=105 && jet_mass_pr <= 130 ) ) && nbjets_csvm_veto == 0 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSBPreselectionCutEXO" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt > 200 && pfMET > 80 && l_pt > 90 && ungroomed_jet_pt > 200 && ( ( jet_mass_pr >=40 && jet_mass_pr <= 65 ) || ( jet_mass_pr >=105 && jet_mass_pr <= 130 )) && nbjets_csvm_veto == 0 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSBPreselectionCutEXO" && (LeptonType == "MuEl" || LeptonType == "muel" || LeptonType == "MuonEle" || LeptonType == "muonele") )
    return Form("issignal && v_pt > 200 && pfMET > 80 && l_pt > 90 && ungroomed_jet_pt > 200 && ( ( jet_mass_pr >=40 && jet_mass_pr <= 65 ) || ( jet_mass_pr >=105 && jet_mass_pr <= 130 )) && nbjets_csvm_veto == 0 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);


  else if(preselectionCutType == "basicSRPreselectionCutEXO" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form(" issignal && v_pt > 200 && pfMET > 40 && l_pt > 50 && ungroomed_jet_pt > 200 && ( jet_mass_pr >=65 && jet_mass_pr <= 105 )"
                " && nbjets_csvm_veto == 0 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSRPreselectionCutEXO" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt > 200 && pfMET > 80 && l_pt > 90 && ungroomed_jet_pt > 200 && ( jet_mass_pr >=65 && jet_mass_pr <= 105 ) && nbjets_csvm_veto == 0 && "
                "( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSRPreselectionCutEXO" && (LeptonType == "MuEl" || LeptonType == "muel" || LeptonType == "MuonEle" || LeptonType == "muonele") )
    return Form("issignal && v_pt > 200 && pfMET > 80 && l_pt > 90 && ungroomed_jet_pt > 200 && ( jet_mass_pr >=65 && jet_mass_pr <= 105 ) && nbjets_csvm_veto == 0 && "
                "( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);


  else if(preselectionCutType == "basicSRSBPreselectionCutEXO" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form("issignal && v_pt > 200 && pfMET > 40 && l_pt > 50 && ungroomed_jet_pt > 200 && ( jet_mass_pr >=40 && jet_mass_pr <= 130 ) && nbjets_csvm_veto == 0"
                " && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSRSBPreselectionCutEXO" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt > 200 && pfMET > 80 && l_pt > 90 && ungroomed_jet_pt > 200 && ( jet_mass_pr >=40 && jet_mass_pr <= 130 ) && nbjets_csvm_veto == 0"
                " && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSRSBPreselectionCutEXO" && (LeptonType == "MuEl" || LeptonType == "muel" || LeptonType == "MuonEle" || LeptonType == "muonele") )
    return Form("issignal && v_pt > 200 && pfMET > 80 && l_pt > 90 && ungroomed_jet_pt > 200 && ( jet_mass_pr >=40 && jet_mass_pr <= 130 ) && nbjets_csvm_veto == 0"
                " && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);


  ///////////////////// Higgs Like Selection                                                                                                                        

  if(preselectionCutType == "basicPreselectionCutHiggs" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form("issignal && v_pt>200 && pfMET>50 && l_pt>30 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && numberJetBin < 2 && nbjets_csvm_veto == 0 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicPreselectionCutHiggs" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt > 200 && pfMET > 70 && l_pt > 35 && ungroomed_jet_pt > 200 && nbjets_csvm_veto == 0 && numberJetBin < 2 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicPreselectionCutHiggs" && (LeptonType == "MuEl" || LeptonType == "muel" || LeptonType == "MuonEle" || LeptonType == "muonele") )
    return Form("issignal && v_pt > 200 && pfMET > 70 && l_pt > 35 && ungroomed_jet_pt > 200 && nbjets_csvm_veto == 0 && numberJetBin < 2 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSBPreselectionCutHiggs" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form("issignal && v_pt>200 && pfMET>50 && l_pt>30 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && nbjets_csvm_veto == 0 && ( (jet_mass_pr > 40 && jet_mass_pr < 65) || (jet_mass_pr > 105 && jet_mass_pr < 130) ) && numberJetBin < 2 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSBPreselectionCutHiggs" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt > 200 && pfMET > 70 && l_pt > 35 && ungroomed_jet_pt > 200 && ( ( jet_mass_pr >=40 && jet_mass_pr <= 60 ) || ( jet_mass_pr >=100 && jet_mass_pr <= 130 )) && nbjets_csvm_veto == 0 && numberJetBin < 2 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSBPreselectionCutHiggs" && (LeptonType == "MuEl" || LeptonType == "muel" || LeptonType == "MuonEle" || LeptonType == "muonele") )
    return Form("issignal && v_pt > 200 && pfMET > 70 && l_pt > 35 && ungroomed_jet_pt > 200 && ( ( jet_mass_pr >=40 && jet_mass_pr <= 65 ) || ( jet_mass_pr >=105 && jet_mass_pr <= 130 )) && nbjets_csvm_veto == 0 && numberJetBin < 2 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);


  else if(preselectionCutType == "basicSRPreselectionCutHiggs" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form("issignal && v_pt>200 && pfMET>50 && l_pt>30 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && nbjets_csvm_veto == 0 && (jet_mass_pr > 65 && jet_mass_pr < 105) && numberJetBin < 2 ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f)",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSRPreselectionCutHiggs" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt > 200 && pfMET > 70 && l_pt > 35 && ungroomed_jet_pt > 200 && ( jet_mass_pr >=65 && jet_mass_pr <= 105 ) && nbjets_csvm_veto == 0 && numberJetBin < 2 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSRPreselectionCutHiggs" && (LeptonType == "MuEl" || LeptonType == "muel" || LeptonType == "MuonEle" || LeptonType == "muonele") )
    return Form("issignal && v_pt > 200 && pfMET > 70 && l_pt > 35 && ungroomed_jet_pt > 200 && ( jet_mass_pr >=65 && jet_mass_pr <= 105 ) && nbjets_csvm_veto == 0 && numberJetBin < 2 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);


  else if(preselectionCutType == "basicSRSBPreselectionCutHiggs" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form("issignal && v_pt>200 && pfMET>50 && l_pt>30 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && nbjets_csvm_veto == 0 && (jet_mass_pr > 40 && jet_mass_pr < 130) && numberJetBin < 2 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSRSBPreselectionCutHiggs" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt > 200 && pfMET > 70 && l_pt > 35 && ungroomed_jet_pt > 200 && ( jet_mass_pr >=40 && jet_mass_pr <= 130 ) && nbjets_csvm_veto == 0 && numberJetBin < 2 &&( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if(preselectionCutType == "basicSRSBPreselectionCutHiggs" && (LeptonType == "MuEl" || LeptonType == "muel" || LeptonType == "MuonEle" || LeptonType == "muonele") )
    return Form("issignal && v_pt > 200 && pfMET > 70 && l_pt > 35 && ungroomed_jet_pt > 200 && ( jet_mass_pr >=40 && jet_mass_pr <= 130 ) && nbjets_csvm_veto == 0 && numberJetBin < 2 &&( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);


  /// VBF Cuts                                                                                                                                                                            

  else if( preselectionCutType == "basicVBFPreselectionCutHiggs" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form("issignal && v_pt>200 && pfMET>50 && l_pt>30 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679&& numberJetBin >= 2 && jet_tau2tau1 < 0.5 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if( preselectionCutType == "basicVBFPreselectionCutHiggs" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt>200 && pfMET>70 && l_pt>35 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679&& numberJetBin >= 2 && jet_tau2tau1 < 0.5 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if( preselectionCutType == "basicVBFPreselectionCutHiggs" && (LeptonType == "MuEl" || LeptonType == "muel" || LeptonType == "MuonEle" || LeptonType == "muonele") )
    return Form("issignal && v_pt>200 && pfMET>70 && l_pt>35 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679&& numberJetBin >= 2 && jet_tau2tau1 < 0.5 && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);


  else if( preselectionCutType == "basicVBFSBHiggs" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form("issignal && v_pt>200 && pfMET>50 && l_pt>30 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679&& numberJetBin >= 2 && jet_tau2tau1 < 0.5 && ((jet_mass_pr > 40 && jet_mass_pr <65) || (jet_mass_pr > 105 && jet_mass_pr < 130)) && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if( preselectionCutType == "basicVBFSBHiggs" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt>200 && pfMET>70 && l_pt>35 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679&& numberJetBin >= 2 && jet_tau2tau1 < 0.5 && ((jet_mass_pr > 40 && jet_mass_pr <65) || (jet_mass_pr > 105 && jet_mass_pr < 130)) && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if( preselectionCutType == "basicVBFSBHiggs" && (LeptonType == "MuEl" || LeptonType == "muel" || LeptonType == "MuonEle" || LeptonType == "muonele") )
    return Form("issignal && v_pt>200 && pfMET>70 && l_pt>35 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679&& numberJetBin >= 2 && jet_tau2tau1 < 0.5 && ((jet_mass_pr > 40 && jet_mass_pr <65) || (jet_mass_pr > 105 && jet_mass_pr < 130)) && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if( preselectionCutType == "basicVBFSRHiggs" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form("issignal && v_pt>200 && pfMET>50 && l_pt>30 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679 && numberJetBin >= 2 && jet_tau2tau1 < 0.5 && (jet_mass_pr > 65 && jet_mass_pr <105)  && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if( preselectionCutType == "basicVBFSRHiggs" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt>200 && pfMET>70 && l_pt>35 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679&& numberJetBin >= 2 && jet_tau2tau1 < 0.5 && (jet_mass_pr > 65 && jet_mass_pr <105)  && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if( preselectionCutType == "basicVBFSRHiggs" && (LeptonType == "MuEl" || LeptonType == "muel" || LeptonType == "MuonEle" || LeptonType == "muonele") )
    return Form("issignal && v_pt>200 && pfMET>70 && l_pt>35 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679&& numberJetBin >= 2 && jet_tau2tau1 < 0.5 && (jet_mass_pr > 65 && jet_mass_pr <105) && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);


  else if( preselectionCutType == "basicVBFSBSRHiggs" && (LeptonType == "Mu" || LeptonType == "mu" || LeptonType == "Muon" || LeptonType == "muon") )
    return Form("issignal && v_pt>200 && pfMET>50 && l_pt>30 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679&& numberJetBin >= 2 && jet_tau2tau1 < 0.5 && (jet_mass_pr > 40 && jet_mass_pr <130)  && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if( preselectionCutType == "basicVBFSBSRHiggs" && (LeptonType == "El" || LeptonType == "el" || LeptonType == "Electron" || LeptonType == "electron") )
    return Form("issignal && v_pt>200 && pfMET>70 && l_pt>35 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679&& numberJetBin >= 2 && jet_tau2tau1 < 0.5 && (jet_mass_pr > 40 && jet_mass_pr <130)  && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else if( preselectionCutType == "basicVBFSBSRHiggs" && (LeptonType == "MuEl" || LeptonType == "muel" || LeptonType == "MuonEle" || LeptonType == "muonele") )
    return Form("issignal && v_pt>200 && pfMET>70 && l_pt>35 && ungroomed_jet_pt>200 && abs(l_eta)<2.4 && vbf_maxpt_j1_bDiscriminatorCSV <=0.679 && vbf_maxpt_j2_bDiscriminatorCSV <=0.679&& numberJetBin >= 2 && jet_tau2tau1 < 0.5 && (jet_mass_pr > 40 && jet_mass_pr <130) && ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);

  else return Form("v_pt > 200 && pfMET > 40 && l_pt > 50 && ungroomed_jet_pt > 200 && nbjets_csvm_veto == 0 ( ungroomed_jet_pt > %f  && ungroomed_jet_pt < %f )",pTJetMin_,pTJetMax_);


}
