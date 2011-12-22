#include "PhysicsTools/MiBiCommonPAT/plugins/SimpleNtupleGen.h"



///--------------
///---- ctor ----

SimpleNtupleGen::SimpleNtupleGen(const edm::ParameterSet& iConfig)
{
  //---- Out file ----
  edm::Service<TFileService> fs;
  outTree_ = fs -> make<TTree>("SimpleNtupleGen","SimpleNtupleGen"); 
  NtupleFactory_ = new NtupleFactory(outTree_);  
  
  
  //---- MC dumpers ----
  mcAnalysis_      = NULL;
  mcAnalysisHiggs_ = NULL;
  mcAnalysisTTBar_ = NULL;
  
  
  //---- Input tags ---- 
  GenParticlesTag_ = iConfig.getParameter<edm::InputTag>("GenParticlesTag");
  GenMetTag_       = iConfig.getParameter<edm::InputTag>("GenMetTag");
  GenJetTag_       = iConfig.getParameter<edm::InputTag>("GenJetTag");
  
  
  //---- flags ----
  saveGenPtHat_   = iConfig.getUntrackedParameter<bool>("saveGenPtHat", true);
  saveGenEle_     = iConfig.getUntrackedParameter<bool>("saveGenEle", true);
  saveGenMu_      = iConfig.getUntrackedParameter<bool>("saveGenMu", true);
  saveGenTau_     = iConfig.getUntrackedParameter<bool>("saveGenTau", true);
  saveGenTauJ_    = iConfig.getUntrackedParameter<bool>("saveGenTauJ", true);
  saveGenMet_     = iConfig.getUntrackedParameter<bool>("saveGenMet", true);
  saveGenJet_     = iConfig.getUntrackedParameter<bool>("saveGenJet", true);
  saveGenTTBar_   = iConfig.getUntrackedParameter<bool>("saveGenTTBar", false);
  saveGenHiggs_   = iConfig.getUntrackedParameter<bool>("saveGenHiggs", true);
  saveGenHiggsWW_ = iConfig.getUntrackedParameter<bool>("saveGenHiggsWW", true);
  
  verbosity_ = iConfig.getUntrackedParameter<bool>("verbosity", false);
  eventType_ = iConfig.getUntrackedParameter<int>("eventType", 1);
  
  
  
  //---- Add branches to ntuple ----  
  NtupleFactory_ -> AddInt("runId"); 
  NtupleFactory_ -> AddInt("lumiId"); 
  NtupleFactory_ -> AddInt("BXId"); 
  NtupleFactory_ -> AddInt("eventId"); 
  NtupleFactory_ -> AddInt("eventNaiveId"); 
  eventNaiveId_ = 0;
  
  if(saveGenPtHat_)
  {
    NtupleFactory_ -> AddFloat("mcPtHat");
  }
  
  if(saveGenEle_)
  {
    NtupleFactory_ -> Add4V("mcEle");
    NtupleFactory_ -> AddFloat("mcEle_charge");
  }
  
  if(saveGenMu_)
  {
    NtupleFactory_ -> Add4V("mcMu");
    NtupleFactory_ -> AddFloat("mcMu_charge");
  }
  
  if(saveGenTau_)
  {
    NtupleFactory_ -> Add4V("mcTau");
    NtupleFactory_ -> AddFloat("mcTau_charge");
  }
  
  if(saveGenTauJ_)
  {
    NtupleFactory_ -> Add4V("mcTauJ");
    NtupleFactory_ -> AddFloat("mcTauJ_charge");
    NtupleFactory_ -> AddInt("mcTauJ_multiplicity");
  }
  
  if(saveGenMet_)
  {
    NtupleFactory_ -> Add4V("mcMet");
  }
  
  if(saveGenJet_)
  {
    NtupleFactory_ -> Add4V("mcJet");
  }
  
  if(saveGenTTBar_)
  {
    NtupleFactory_->Add4V("mcT1");    
    NtupleFactory_->AddFloat("mcT1_charge");    
    NtupleFactory_->Add4V("mcT2");    
    NtupleFactory_->AddFloat("mcT2_charge");    
    
    NtupleFactory_->Add4V("mcB1");    
    NtupleFactory_->AddFloat("mcB1_charge");    
    NtupleFactory_->Add4V("mcB2");    
    NtupleFactory_->AddFloat("mcB2_charge");   
    
    NtupleFactory_->Add4V("mcV1");         
    NtupleFactory_->AddFloat("mcV1_charge");    
    NtupleFactory_->AddFloat("mcV1_pdgId");    
    
    NtupleFactory_->Add4V("mcV2");         
    NtupleFactory_->AddFloat("mcV2_charge");    
    NtupleFactory_->AddFloat("mcV2_pdgId");  
    
    NtupleFactory_->Add4V("mcF1_fromV1");   
    NtupleFactory_->AddFloat("mcF1_fromV1_charge");    
    NtupleFactory_->AddFloat("mcF1_fromV1_pdgId");  
    
    NtupleFactory_->Add4V("mcF2_fromV1");         
    NtupleFactory_->AddFloat("mcF2_fromV1_charge");    
    NtupleFactory_->AddFloat("mcF2_fromV1_pdgId");  
    
    NtupleFactory_->Add4V("mcF1_fromV2");         
    NtupleFactory_->AddFloat("mcF1_fromV2_charge");    
    NtupleFactory_->AddFloat("mcF1_fromV2_pdgId");  
    
    NtupleFactory_->Add4V("mcF2_fromV2");         
    NtupleFactory_->AddFloat("mcF2_fromV2_charge");    
    NtupleFactory_->AddFloat("mcF2_fromV2_pdgId");    
  }
  
  if(saveGenHiggs_)
  {
    NtupleFactory_->Add4V("mcH");    
    NtupleFactory_->AddFloat("mcH_charge");    
    
    if(saveGenHiggsWW_)
    {
      NtupleFactory_->Add4V("mcQ1_tag");    
      NtupleFactory_->AddFloat("mcQ1_tag_charge");    
      NtupleFactory_->AddFloat("mcQ1_tag_pdgId");  
      
      NtupleFactory_->Add4V("mcQ2_tag");         
      NtupleFactory_->AddFloat("mcQ2_tag_charge");    
      NtupleFactory_->AddFloat("mcQ2_tag_pdgId");  
      
      NtupleFactory_->Add4V("mcV1");         
      NtupleFactory_->AddFloat("mcV1_charge");    
      NtupleFactory_->AddFloat("mcV1_pdgId");    
      
      NtupleFactory_->Add4V("mcV2");         
      NtupleFactory_->AddFloat("mcV2_charge");    
      NtupleFactory_->AddFloat("mcV2_pdgId");  
      
      NtupleFactory_->Add4V("mcF1_fromV1");   
      NtupleFactory_->AddFloat("mcF1_fromV1_charge");    
      NtupleFactory_->AddFloat("mcF1_fromV1_pdgId");  
      
      NtupleFactory_->Add4V("mcF2_fromV1");         
      NtupleFactory_->AddFloat("mcF2_fromV1_charge");    
      NtupleFactory_->AddFloat("mcF2_fromV1_pdgId");  
      
      NtupleFactory_->Add4V("mcF1_fromV2");         
      NtupleFactory_->AddFloat("mcF1_fromV2_charge");    
      NtupleFactory_->AddFloat("mcF1_fromV2_pdgId");  
      
      NtupleFactory_->Add4V("mcF2_fromV2");         
      NtupleFactory_->AddFloat("mcF2_fromV2_charge");    
      NtupleFactory_->AddFloat("mcF2_fromV2_pdgId");  
    }
  }
}






///--------------
///---- dtor ----

SimpleNtupleGen::~SimpleNtupleGen()
{
  NtupleFactory_->WriteNtuple();
  delete NtupleFactory_;
}






///-------------
///---- Gen ----

void SimpleNtupleGen::fillGenPtHatInfo(const edm::Event & iEvent, const edm::EventSetup & iESetup) 
{
  //std::cout << "SimpleNtupleGen::fillGenPtHatInfo::begin" << std::endl; 
  
  edm::Handle< GenEventInfoProduct > GenInfoHandle;
  iEvent.getByLabel( "generator", GenInfoHandle);
  float ptHat = ( GenInfoHandle->hasBinningValues() ? (GenInfoHandle->binningValues())[0] : 0.0);
  
  NtupleFactory_ -> FillFloat("mcPtHat", ptHat);
  
  //std::cout << "SimpleNtupleGen::fillGenPtHatInfo::end" << std::endl; 
}



void SimpleNtupleGen::fillGenEleInfo(const edm::Event & iEvent, const edm::EventSetup & iESetup) 
{
  //std::cout << "SimpleNtupleGen::fillGenEleInfo::begin" << std::endl; 
  
  std::vector<const reco::Candidate*> genEle = mcAnalysis_ -> GetMcE();
  
  for(unsigned int eleIt = 0; eleIt < genEle.size(); ++eleIt)
  {
    NtupleFactory_ -> Fill4V("mcEle", (genEle.at(eleIt))->p4());
    NtupleFactory_ -> FillFloat("mcEle_charge", (genEle.at(eleIt))->charge());
  }
  
  //std::cout << "SimpleNtupleGen::fillGenEleInfo::end" << std::endl; 
}



void SimpleNtupleGen::fillGenMuInfo(const edm::Event & iEvent, const edm::EventSetup & iESetup) 
{
  //std::cout << "SimpleNtupleGen::fillGenMuInfo::begin" << std::endl; 
  
  std::vector<const reco::Candidate*> genMu = mcAnalysis_ -> GetMcMu();
  
  for(unsigned int muIt = 0; muIt < genMu.size(); ++muIt)
  {
    NtupleFactory_ -> Fill4V("mcMu", (genMu.at(muIt))->p4());
    NtupleFactory_ -> FillFloat("mcMu_charge", (genMu.at(muIt))->charge());
  }
  
  //std::cout << "SimpleNtupleGen::fillGenMuInfo::end" << std::endl; 
}



void SimpleNtupleGen::fillGenTauInfo(const edm::Event & iEvent, const edm::EventSetup & iESetup) 
{
  //std::cout << "SimpleNtupleGen::fillGenTauInfo::begin" << std::endl; 
  
  std::vector<const reco::Candidate*> genTau = mcAnalysis_ -> GetMcTau();
  
  for(unsigned int tauIt = 0; tauIt < genTau.size(); ++tauIt)
  {
    NtupleFactory_ -> Fill4V("mcTau", (genTau.at(tauIt))->p4());
    NtupleFactory_ -> FillFloat("mcTau_charge", (genTau.at(tauIt))->charge());
  }
  
  //std::cout << "SimpleNtupleGen::fillGenTauInfo::end" << std::endl; 
}



void SimpleNtupleGen::fillGenTauJInfo(const edm::Event & iEvent, const edm::EventSetup & iESetup) 
{
  //std::cout << "SimpleNtupleGen::fillGenTauJInfo::begin" << std::endl; 
  
  std::vector<std::vector<const reco::Candidate*> > genTauJ = mcAnalysis_ -> GetMcTauJ();
  
  for(unsigned int tauJIt = 0; tauJIt < genTauJ.size(); ++tauJIt)
  {
    std::vector<const reco::Candidate*> components = genTauJ.at(tauJIt);
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > p4(0.,0.,0.,0.);
    float charge = 0.;
    
    for(unsigned int compIt = 0; compIt < components.size(); ++compIt)
    {
      p4 = p4 + components.at(compIt)->p4();
      charge += components.at(compIt) -> charge();
    }
    
    NtupleFactory_ -> Fill4V("mcTauJ", p4);
    NtupleFactory_ -> FillFloat("mcTauJ_charge", charge);
    NtupleFactory_ -> FillInt("mcTauJ_multiplicity", components.size());
  }
  
  //std::cout << "SimpleNtupleGen::fillGenTauJInfo::end" << std::endl; 
}



void SimpleNtupleGen::fillGenMetInfo(const edm::Event & iEvent, const edm::EventSetup & iESetup) 
{
  //std::cout << "SimpleNtuple::fillGenMetInfo::begin" << std::endl;
  
  edm::Handle<edm::View<reco::GenMET> > genMetHandle;
  iEvent.getByLabel(GenMetTag_, genMetHandle);
  edm::View<reco::GenMET> genMet = *genMetHandle;

  NtupleFactory_->Fill4V("mcMet", genMet.at(0).p4());
  
  //std::cout << "SimpleNtuple::fillGenMetInfo::end" << std::endl;
}



void SimpleNtupleGen::fillGenJetInfo(const edm::Event & iEvent, const edm::EventSetup & iESetup) 
{
  //std::cout << "SimpleNtupleGen::fillGenJetInfo::begin" << std::endl; 
  
  edm::Handle<edm::View<reco::GenJet> > genJetHandle;
  iEvent.getByLabel(GenJetTag_, genJetHandle);
  edm::View<reco::GenJet> genJets = *genJetHandle;
  
  // loop on jets
  for(unsigned int genJetIt = 0; genJetIt < genJets.size(); ++genJetIt)
  {
    reco::GenJet genJet = genJets.at(genJetIt);
    
    NtupleFactory_ -> Fill4V("mcJet",genJet.p4());
  }
  
  //std::cout << "SimpleNtupleGen::fillGenJetInfo::end" << std::endl; 
}



void SimpleNtupleGen::fillGenHiggsInfo(const edm::Event & iEvent, const edm::EventSetup & iESetup) 
{
  //std::cout << "SimpleNtupleGen::fillGenHiggsDecayInfo" << std::endl; 
  
  bool isValid = mcAnalysisHiggs_ -> isValid();
  
  if( (eventType_ == 0) && (isValid == true) )
  {
    NtupleFactory_->Fill4V("mcH",mcAnalysisHiggs_ -> mcH()->p4());
    NtupleFactory_->FillFloat("mcH_charge",mcAnalysisHiggs_ -> mcH()->charge());
    
    if(saveGenHiggsWW_)
    {  
      NtupleFactory_->Fill4V("mcQ1_tag",mcAnalysisHiggs_ -> mcQ1_tag()->p4());
      NtupleFactory_->FillFloat("mcQ1_tag_charge",mcAnalysisHiggs_ -> mcQ1_tag()->charge());
      NtupleFactory_->FillFloat("mcQ1_tag_pdgId",mcAnalysisHiggs_ -> mcQ1_tag()->pdgId());
      
      NtupleFactory_->Fill4V("mcQ2_tag",mcAnalysisHiggs_ -> mcQ2_tag()->p4());
      NtupleFactory_->FillFloat("mcQ2_tag_charge",mcAnalysisHiggs_ -> mcQ2_tag()->charge());
      NtupleFactory_->FillFloat("mcQ2_tag_pdgId",mcAnalysisHiggs_ -> mcQ2_tag()->pdgId());
      
      NtupleFactory_->Fill4V("mcV1",mcAnalysisHiggs_ -> mcV1()->p4());
      NtupleFactory_->FillFloat("mcV1_charge",mcAnalysisHiggs_ -> mcV1()->charge());
      NtupleFactory_->FillFloat("mcV1_pdgId",mcAnalysisHiggs_ -> mcV1()->pdgId());
      
      NtupleFactory_->Fill4V("mcV2",mcAnalysisHiggs_ -> mcV2()->p4());
      NtupleFactory_->FillFloat("mcV2_charge",mcAnalysisHiggs_ -> mcV2()->charge());
      NtupleFactory_->FillFloat("mcV2_pdgId",mcAnalysisHiggs_ -> mcV2()->pdgId());
      
      NtupleFactory_->Fill4V("mcF1_fromV1",mcAnalysisHiggs_ -> mcF1_fromV1()->p4());
      NtupleFactory_->FillFloat("mcF1_fromV1_charge",mcAnalysisHiggs_ -> mcF1_fromV1()->charge());
      NtupleFactory_->FillFloat("mcF1_fromV1_pdgId",mcAnalysisHiggs_ -> mcF1_fromV1()->pdgId());
      
      NtupleFactory_->Fill4V("mcF2_fromV1",mcAnalysisHiggs_ -> mcF2_fromV1()->p4());
      NtupleFactory_->FillFloat("mcF2_fromV1_charge",mcAnalysisHiggs_ -> mcF2_fromV1()->charge());
      NtupleFactory_->FillFloat("mcF2_fromV1_pdgId",mcAnalysisHiggs_ -> mcF2_fromV1()->pdgId());
      
      NtupleFactory_->Fill4V("mcF1_fromV2",mcAnalysisHiggs_ -> mcF1_fromV2()->p4());
      NtupleFactory_->FillFloat("mcF1_fromV2_charge",mcAnalysisHiggs_ -> mcF1_fromV2()->charge());
      NtupleFactory_->FillFloat("mcF1_fromV2_pdgId",mcAnalysisHiggs_ -> mcF1_fromV2()->pdgId());
      
      NtupleFactory_->Fill4V("mcF2_fromV2",mcAnalysisHiggs_ -> mcF2_fromV2()->p4());
      NtupleFactory_->FillFloat("mcF2_fromV2_charge",mcAnalysisHiggs_ -> mcF2_fromV2()->charge());
      NtupleFactory_->FillFloat("mcF2_fromV2_pdgId",mcAnalysisHiggs_ -> mcF2_fromV2()->pdgId());
    }
  } 

}
 


void SimpleNtupleGen::fillGenTTBarInfo(const edm::Event & iEvent, const edm::EventSetup & iESetup) 
{
  //std::cout << "SimpleNtupleGen::fillGenTTBarInfo" << std::endl;
  
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel(GenParticlesTag_, genParticles);
  
  bool isValid = mcAnalysisTTBar_->isValid();
   
  if( (eventType_ == 0) && (isValid == true) )
  {
    NtupleFactory_->Fill4V("mcT1",mcAnalysisTTBar_->mcT1()->p4());
    NtupleFactory_->FillFloat("mcT1_charge",mcAnalysisTTBar_->mcT1()->charge());
    
    NtupleFactory_->Fill4V("mcT2",mcAnalysisTTBar_->mcT2()->p4());
    NtupleFactory_->FillFloat("mcT2_charge",mcAnalysisTTBar_->mcT2()->charge());
    
    NtupleFactory_->Fill4V("mcB1",mcAnalysisTTBar_->mcB1()->p4());
    NtupleFactory_->FillFloat("mcB1_charge",mcAnalysisTTBar_->mcB1()->charge());
    
    NtupleFactory_->Fill4V("mcB2",mcAnalysisTTBar_->mcB2()->p4());
    NtupleFactory_->FillFloat("mcB2_charge",mcAnalysisTTBar_->mcB2()->charge());
    
    NtupleFactory_->Fill4V("mcV1",mcAnalysisTTBar_->mcV1()->p4());
    NtupleFactory_->FillFloat("mcV1_charge",mcAnalysisTTBar_->mcV1()->charge());
    NtupleFactory_->FillFloat("mcV1_pdgId",mcAnalysisTTBar_->mcV1()->pdgId());
    
    NtupleFactory_->Fill4V("mcV2",mcAnalysisTTBar_->mcV2()->p4());
    NtupleFactory_->FillFloat("mcV2_charge",mcAnalysisTTBar_->mcV2()->charge());
    NtupleFactory_->FillFloat("mcV2_pdgId",mcAnalysisTTBar_->mcV2()->pdgId());
    
    NtupleFactory_->Fill4V("mcF1_fromV1",mcAnalysisTTBar_->mcF1_fromV1()->p4());
    NtupleFactory_->FillFloat("mcF1_fromV1_charge",mcAnalysisTTBar_->mcF1_fromV1()->charge());
    NtupleFactory_->FillFloat("mcF1_fromV1_pdgId",mcAnalysisTTBar_->mcF1_fromV1()->pdgId());
    
    NtupleFactory_->Fill4V("mcF2_fromV1",mcAnalysisTTBar_->mcF2_fromV1()->p4());
    NtupleFactory_->FillFloat("mcF2_fromV1_charge",mcAnalysisTTBar_->mcF2_fromV1()->charge());
    NtupleFactory_->FillFloat("mcF2_fromV1_pdgId",mcAnalysisTTBar_->mcF2_fromV1()->pdgId());
    
    NtupleFactory_->Fill4V("mcF1_fromV2",mcAnalysisTTBar_->mcF1_fromV2()->p4());
    NtupleFactory_->FillFloat("mcF1_fromV2_charge",mcAnalysisTTBar_->mcF1_fromV2()->charge());
    NtupleFactory_->FillFloat("mcF1_fromV2_pdgId",mcAnalysisTTBar_->mcF1_fromV2()->pdgId());
    
    NtupleFactory_->Fill4V("mcF2_fromV2",mcAnalysisTTBar_->mcF2_fromV2()->p4());
    NtupleFactory_->FillFloat("mcF2_fromV2_charge",mcAnalysisTTBar_->mcF2_fromV2()->charge());
    NtupleFactory_->FillFloat("mcF2_fromV2_pdgId",mcAnalysisTTBar_->mcF2_fromV2()->pdgId());
  }

}






///-----------------
///---- analyze ----

void SimpleNtupleGen::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  ++eventNaiveId_;
  
  NtupleFactory_->FillInt("runId", iEvent.id().run());
  NtupleFactory_->FillInt("lumiId", iEvent.luminosityBlock());
  NtupleFactory_->FillInt("BXId", iEvent.bunchCrossing());
  NtupleFactory_->FillInt("eventId", iEvent.id().event());
  NtupleFactory_->FillInt("eventNaiveId", eventNaiveId_);
  
  
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel(GenParticlesTag_, genParticles);
  
  mcAnalysis_ = new MCDumper(genParticles, verbosity_);
  
  if(saveGenTTBar_)
    mcAnalysisTTBar_ = new MCDumperTTBar(genParticles, eventType_, verbosity_);
  
  if(saveGenHiggs_)
    mcAnalysisHiggs_ = new MCDumperHiggs(genParticles, eventType_, verbosity_);
  
  
  
  ///---- fill GenPtHat ---- 
  if(saveGenPtHat_) fillGenPtHatInfo(iEvent, iSetup);
  
  ///---- fill GenParticles ----
  if(saveGenEle_)   fillGenEleInfo(iEvent, iSetup);
  if(saveGenMu_)    fillGenMuInfo(iEvent, iSetup);
  if(saveGenTau_)   fillGenTauInfo(iEvent, iSetup);
  if(saveGenTauJ_)  fillGenTauJInfo(iEvent, iSetup);
  if(saveGenJet_)   fillGenJetInfo(iEvent, iSetup);
  if(saveGenMet_)   fillGenMetInfo(iEvent, iSetup);
  
  ///---- fill Higgs/TTBar ----
  if(saveGenTTBar_) fillGenTTBarInfo(iEvent, iSetup);
  if(saveGenHiggs_) fillGenHiggsInfo(iEvent, iSetup);
  if(saveGenHiggsWW_) fillGenHiggsInfo(iEvent, iSetup);
  
  ///---- save the entry of the tree ----
  NtupleFactory_->FillNtuple();
}






///===================================
DEFINE_FWK_MODULE(SimpleNtupleGen) ;
