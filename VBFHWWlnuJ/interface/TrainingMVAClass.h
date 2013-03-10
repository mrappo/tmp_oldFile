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

#include "TMVA/MsgLogger.h"
#include "TMVA/Config.h"
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"

class TrainingMVAClass {

 public :

  TrainingMVAClass(){};

  TrainingMVAClass(const std::vector<TFile*> & signalFileList, const std::vector<TFile*> & backgroundFileList, const std::string & TreeName, 
                   const std::string & outputFilePath , const std::string & outputFileName, const std::string & Label );

  TrainingMVAClass(const std::vector<TTree*> & signalTreeList, const std::vector<TTree*> & backgroundTreeList,const std::string & TreeName,
                   const std::string & outputFilePath , const std::string & outputFileName, const std::string & Label);

  ~TrainingMVAClass();

  void BookMVATrees (const std::vector<double> & signalGlobalWeight, const std::vector<double> & backgroundGlobalWeight) ;

  void AddTrainingVariables ( const std::vector<std::string> & mapTrainingVariables, const std::vector<std::string> & mapSpectatorVariables, const std::string & weightString  );	       

  void AddPrepareTraining (const std::string & cutString, const int & nTraining = 0, const int & nTesting = 0, const std::string & splitMode = "Random", 
                           const std::string & NormMode = "NumEvents");


  void BookandTrainRectangularCuts    ( const std::string & FitMethod );
  void BookandTrainLikelihood         ( const std::string & LikelihoodType = "Likelihood");
  void BookandTrainLinearDiscriminant ();
  void BookandTrainFisherDiscriminant ();

  void BookandTrainMLP                ( const int & nCycles = 600, const std::string & HiddenLayers = "N+5", const std::string & NeuronType = "sigmoid",
                                        const std::string & TrainingMethod = "BP", const int & TestRate = 5, const int & ConvergenceTests = 5);
  void BookandTrainBDT                ( const int & NTrees = 500, const std::string & BoostType = "AdaBoost", const float & AdaBoostBeta = 0.5, 
                                        const std::string & PruneMethod = "CostComplexity", const int & PruneStrength=30, 
                                        const int & MaxDepth = 5, const std::string & SeparationType = "GiniIndex");
  void BookandTrainBDTG               ( const int & NTrees = 500, const float & GradBaggingFraction = 0.5, const std::string & PruneMethod = "CostComplexity", const int & PruneStrength = 30,
                                        const int & MaxDepth = 5, const std::string & SeparationType = "GiniIndex");

  void BookandTrainBDTF               ( const int & NTrees = 500, const float & GradBaggingFraction = 0.5, const std::string & PruneMethod = "CostComplexity", const int & PruneStrength = 30,
                                        const int & MaxDepth = 5, const std::string & SeparationType = "GiniIndex");

  void PrintTrainingResults ();
 

  void SetSignalTree (const std::vector<TFile*> & signalFileList,  const std::string & TreeName = "WJet");
  void SetSignalTree (const std::vector<TTree*> & signalTreeList);

  void SetBackgroundTree (const std::vector<TFile*> & backgroundFileList, const std::string & TreeName = "WJet");
  void SetBackgroundTree (const std::vector<TTree*> & backgroundTreeList);

  void SetTrainingVariables  (const std::vector<std::string> & mapTrainingVariables);

  void SetSpectatorVariables (const std::vector<std::string> & mapSpectatorVariables);

  void SetOutputFile ( const std::string & outputFilePath , const std::string & outputFileName );

  void SetTreeName ( const std::string & TreeName );

  void SetLabel ( const std::string & Label );

  void SetGlobalSampleWeight (const std::vector<double> & signalGlobalWeight, const std::vector<double> & backgroundGlobalWeight) ;

  void SetEventWeight ( const std::string & weightString ) ;


 private : 

  std::vector<TTree*> signalTreeList_ ;
  std::vector<TTree*> backgroundTreeList_ ;
  
  std::vector<std::string> mapTrainingVariables_ ;
  std::vector<std::string> mapSpectatorVariables_ ;

  std::vector<double> signalGlobalWeight_ ;
  std::vector<double> backgroundGlobalWeight_ ;
   
  std::string TreeName_ ;
  std::string Label_ ;

  std::string outputFilePath_ ;
  std::string outputFileName_ ;

  std::vector<std::string> outputFileNameComplete_;

  TCut* preselectionCut_ ;

  TFile * outputFile_ ;

  TMVA::Factory* factory_ ; 
};

#endif
