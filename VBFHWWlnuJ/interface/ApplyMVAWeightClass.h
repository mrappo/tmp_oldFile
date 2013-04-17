#ifndef ApplyMVAWeightClass_h
#define ApplyMVAWeightClass_h

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

class ApplyMVAWeightClass {

 public :

  ApplyMVAWeightClass(){};

  ApplyMVAWeightClass(const std::vector<TFile*> & SampleFileList, const std::string & TreeName,
                      const std::string & InputFilePath , const std::string & Label );

  ApplyMVAWeightClass(const std::vector<TTree*> & SampleTreeList, const std::string & TreeName,
                      const std::string & InputFilePath , const std::string & Label);

  ~ApplyMVAWeightClass();

  void AddTrainingVariables ( const std::vector<std::string> & mapTrainingVariables, const std::vector<std::string> & mapSpectatorVariables);

  void AddPrepareReader (const std::string & LeptonType, const std::string & preselectionCutType,
                         std::vector<double> * JetPtBinOfTraining = NULL , const int & pTBin = 0);

  void BookMVAWeight (const std::string & methodName, const std::string & weightFile, const std::string & nameBranch) ;

  void FillMVAWeight (const std::string & LeptonType, const std::string & preselectionCutType) ;

  void SetInputTree (const std::vector<TFile*> & SampleFileList,  const std::string & TreeName = "WJet");

  void SetInputTree (const std::vector<TTree*> & SampleTreeList);

  void SetTrainingVariables  (const std::vector<std::string> & mapTrainingVariables);

  void SetSpectatorVariables (const std::vector<std::string> & mapSpectatorVariables);

  void SetInputFilePath ( const std::string & InputFilePath);

  void SetTreeName ( const std::string & TreeName );

  void SetLabel ( const std::string & Label );

  void SetWeightFile ( const std::string & weightFile ) ;

  void SetReaderTree ();
  void SetReaderTree (const std::vector<TTree*> & SampleTreeList);


 private :

  double pTJetMin_ ;

  double pTJetMax_ ;

  std::vector<TTree*> SampleTreeList_ ;

  std::vector<std::string> mapTrainingVariables_ ;
  std::vector<std::string> mapSpectatorVariables_ ;

  std::vector<Float_t> *setTrainingVariables_ ;
  std::vector<Float_t> *setSpectatorVariables_ ;

  Float_t weight_ ;

  std::string TreeName_ ;
  std::string Label_ ;

  std::string inputFilePath_ ;

  std::string methodName_ ;
  std::string weightFile_ ;
  std::string nameBranch_ ;
 
  TBranch* newBranch_ ;

  std::string preselectionCutType_ ;
  std::string LeptonType_ ;

  TMVA::Reader*  reader_ ;

  std::vector<treeReader*> treeReader_ ;

};

#endif
