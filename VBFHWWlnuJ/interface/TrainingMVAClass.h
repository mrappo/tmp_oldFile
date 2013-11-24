#ifndef TrainingMVAClass_h
#define TrainingMVAClass_h

#include <iostream>
#include <cstdlib>
#include <string>
#include <map>
#include <vector>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TCut.h"

#include "ntpleUtils.h"

#include "treeReader.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/MsgLogger.h"
#include "TMVA/Config.h"
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#endif

class TrainingMVAClass {

 public :

  // default constructor
  TrainingMVAClass(){};

  // constructor from list of files
  TrainingMVAClass(const std::vector<TFile*> & signalFileList, const std::vector<TFile*> & backgroundFileList, const std::string & TreeName, 
                   const std::string & outputFilePath , const std::string & outputFileName, const std::string & Label );

  // constructor from list of trees
  TrainingMVAClass(const std::vector<TTree*> & signalTreeList, const std::vector<TTree*> & backgroundTreeList,const std::string & TreeName,
                   const std::string & outputFilePath , const std::string & outputFileName, const std::string & Label);

  // default de-constructor
  ~TrainingMVAClass();

  // Set global event weight for signal and background
  void BookMVATrees (const std::vector<double> & signalGlobalWeight, const std::vector<double> & backgroundGlobalWeight) ;

  // Set Training and Spectator Variables
  void AddTrainingVariables ( const std::vector<std::string> & mapTrainingVariables, const std::vector<std::string> & mapSpectatorVariables);	       

  // Set the pre-selection cut for the trainning dependig on: lepton tye = muon,electron or MuonEle, type of preselection, weight, pTbin and splitting 
  void AddPrepareTraining (const std::string & LeptonType, const std::string & preselectionCutType, 
                           const std::string & weightString , std::vector<double> * JetPtBinOfTraining = NULL , const int & pTBin = 0,
                           const int & nTraining = 0, const int & nTesting = 0, const std::string & splitMode = "Random", 
                           const std::string & NormMode = "NumEvents");


  // Train rectangular cut methods
  void BookandTrainRectangularCuts    ( const std::string & FitMethod );
  // Train Likelihood
  void BookandTrainLikelihood         ( const std::string & LikelihoodType = "Likelihood");
  // Train Linear Discriminant
  void BookandTrainLinearDiscriminant ();
  // Train Fisher Discriminant
  void BookandTrainFisherDiscriminant ();

  // Train MLP
  void BookandTrainMLP       ( const int & nCycles = 1000, const std::string & HiddenLayers = "N+5", const std::string & NeuronType = "sigmoid",const std::string & TrainingMethod = "BP",                               const int & TestRate = 10, const int & ConvergenceTests = 10, const std::string & EstimatorType = "tanh");

  // Train Clemont Ferrand ANN
  void BookandTrainCFMlpANN           ( const int & nCycles = 1000, const std::string & HiddenLayers = "N+5") ;

  // Train TMVA ANN
  void BookandTrainTMlpANN            ( const int & nCycles = 500, const std::string & HiddenLayers = "N+5",const std::string & TrainingMethod = "BFGS", const float & ValidationFraction=0.3);

  // Train BDT
  void BookandTrainBDT                ( const int & NTrees = 500, const bool & optimizeMethod = false, const std::string & BoostType = "AdaBoost", const float & AdaBoostBeta = 0.5, 
                                        const std::string & PruneMethod = "CostComplexity", const int & PruneStrength=5, 
                                        const int & MaxDepth = 5, const std::string & SeparationType = "GiniIndex");

  // Train Gradient BDT
  void BookandTrainBDTG               ( const int & NTrees = 2000, const bool & optimizeMethod = false, const float & GradBaggingFraction = 0.5, 
                                        const std::string & PruneMethod = "CostComplexity", const int & PruneStrength = 50,
                                        const int & MaxDepth = 5, const std::string & SeparationType = "GiniIndex");
  
  // Train Mit-Fisher BDT
  void BookandTrainBDTF               ( const int & NTrees = 500, const bool & optimizeMethod = false, const std::string & BoostType = "AdaBoost", 
                                        const float & AdaBoostBeta = 0.5, const std::string & PruneMethod = "NoPruning", const int & PruneStrength = 5,
                                        const int & MaxDepth = 5, const std::string & SeparationType = "GiniIndex");

  // Get The Preselection Cut string
  TString GetPreselectionCut (const std::string & LeptonType,const std::string & preselectionCutType = "none") ;

  // In order to do TMVA plots for training + testing results
  void PrintTrainingResults ();

  // Close the output file
  void CloseTrainingAndTesting (){ outputFile_->Close();}

  // Set Signal Tree giving file
  void SetSignalTree (const std::vector<TFile*> & signalFileList,  const std::string & TreeName = "WJet");
  // Set Signal Tree giving tree
  void SetSignalTree (const std::vector<TTree*> & signalTreeList);

  // Set Background Tree giving file
  void SetBackgroundTree (const std::vector<TFile*> & backgroundFileList, const std::string & TreeName = "WJet");
  // Set Background Tree giving tree
  void SetBackgroundTree (const std::vector<TTree*> & backgroundTreeList);

  // Set the training variables name
  void SetTrainingVariables  (const std::vector<std::string> & mapTrainingVariables);

  // Set the spectator variables name
  void SetSpectatorVariables (const std::vector<std::string> & mapSpectatorVariables);

  // Set the output file Name
  void SetOutputFile ( const std::string & outputFilePath , const std::string & outputFileName );

  // Set the name of the tree
  void SetTreeName ( const std::string & TreeName );

  // Set Label for the output file
  void SetLabel ( const std::string & Label );

  // Set Global Event re-weight due to the luminosity
  void SetGlobalSampleWeight (const std::vector<double> & signalGlobalWeight, const std::vector<double> & backgroundGlobalWeight) ;

  // Set Event re-weight : pile-Up, efficiency, cps, interference, btag .. etc
  void SetEventWeight ( const std::string & weightString ) ;

 private : 

  // pT rabge of the training
  double pTJetMin_ ;
  double pTJetMax_ ;

  // list of trees for signal and background
  std::vector<TTree*> signalTreeList_ ;
  std::vector<TTree*> backgroundTreeList_ ;
  
  // list of input and spectator variables
  std::vector<std::string> mapTrainingVariables_ ;
  std::vector<std::string> mapSpectatorVariables_ ;

  // Global re-weight for luminosity
  std::vector<double> signalGlobalWeight_ ;
  std::vector<double> backgroundGlobalWeight_ ;
   
  // TreeName
  std::string TreeName_ ;
  // Label
  std::string Label_ ;

  // outputFilePath
  std::string outputFilePath_ ;
  // output Name
  std::string outputFileName_ ;
  // output Complete Name = path + Name
  std::string outputFileNameComplete_ ;

  // Name of the final file xml with the weights
  std::map<std::string,std::string> outputFileWeightName_ ;

  // preselection cut to apply
  TCut*  preselectionCut_ ;

  // output file
  TFile* outputFile_ ;

  // factory object
  TMVA::Factory* factory_ ; 

};

#endif
