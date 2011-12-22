#ifndef SimpleNtupleGen_h
#define SimpleNtupleGen_h

// cmssw include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"                                                                                                                            
#include "CommonTools/UtilAlgos/interface/TFileService.h" 

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

// user include files
#include "PhysicsTools/NtupleUtils/interface/NtupleFactory.h"
#include "PhysicsTools/MiBiCommonPAT/interface/MCDumper.h"
#include "PhysicsTools/MiBiCommonPAT/interface/MCDumperTTBar.h"
#include "PhysicsTools/MiBiCommonPAT/interface/MCDumperHiggs.h"

// root/c++ include files
#include "TTree.h"
#include <iostream>
#include <vector>






//---------------------------
//---- class declaration ----
//---------------------------

class SimpleNtupleGen : public edm::EDAnalyzer {
 
 public:
  
  //! ctor
  explicit SimpleNtupleGen(const edm::ParameterSet&);
  
  //! dtor
  ~SimpleNtupleGen();
  
  
  
 private:

  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  
  void fillGenPtHatInfo(const edm::Event & iEvent, const edm::EventSetup & iESetup) ; 
  void fillGenEleInfo  (const edm::Event & iEvent, const edm::EventSetup & iESetup) ;
  void fillGenMuInfo   (const edm::Event & iEvent, const edm::EventSetup & iESetup) ;
  void fillGenTauInfo  (const edm::Event & iEvent, const edm::EventSetup & iESetup) ;
  void fillGenTauJInfo (const edm::Event & iEvent, const edm::EventSetup & iESetup) ;
  void fillGenMetInfo  (const edm::Event & iEvent, const edm::EventSetup & iESetup) ;
  void fillGenJetInfo  (const edm::Event & iEvent, const edm::EventSetup & iESetup) ;
  void fillGenHiggsInfo(const edm::Event & iEvent, const edm::EventSetup & iESetup) ;
  void fillGenTTBarInfo(const edm::Event & iEvent, const edm::EventSetup & iESetup) ;
  
  
  
  TTree* outTree_;
  NtupleFactory* NtupleFactory_;
  
  ///---- input tag ----
  edm::InputTag GenParticlesTag_;
  edm::InputTag GenMetTag_;
  edm::InputTag GenJetTag_;
  
  
  ///---- save MC Info ----
  bool saveGenPtHat_;
  bool saveGenEle_;
  bool saveGenMu_;
  bool saveGenTau_;
  bool saveGenTauJ_;
  bool saveGenMet_;
  bool saveGenJet_;
  bool saveGenTTBar_;
  bool saveGenHiggs_;
  bool saveGenHiggsWW_;
  
  int eventType_;  //---- 0 = signal    1 = background 
  bool verbosity_; //---- true = loquacious    false = silence  
  int eventNaiveId_;
  
  MCDumper* mcAnalysis_;
  MCDumperHiggs* mcAnalysisHiggs_;
  MCDumperTTBar* mcAnalysisTTBar_;
};

#endif
