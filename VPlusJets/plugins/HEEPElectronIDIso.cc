

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"



class HEEPElectronIDIso : public edm::EDFilter {
 public:
  explicit HEEPElectronIDIso(const edm::ParameterSet & );
  ~HEEPElectronIDIso();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


 private:
  virtual void beginJob() ;
  virtual void endJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);

  virtual bool beginRun(edm::Run const&, edm::EventSetup const&);
  virtual bool endRun(edm::Run const&, edm::EventSetup const&);
  virtual bool beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual bool endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  edm::InputTag eleLabel_;
  int           maxNumber_;
  int           minNumber_;

};

HEEPElectronIDIso::HEEPElectronIDIso(const edm::ParameterSet &  iConfig){

  if( iConfig.existsAs<edm::InputTag>("electronCollection") )
    eleLabel_=iConfig.getParameter< edm::InputTag >("electronCollection");
  else eleLabel_= edm::InputTag("heepPatElectrons");

  if( iConfig.existsAs<int>("maxNumber") )
      maxNumber_=iConfig.getUntrackedParameter<int>("maxNumber");
  else maxNumber_= 99999 ;

  if( iConfig.existsAs<int>("minNumber") )
      minNumber_=iConfig.getUntrackedParameter<int>("minNumber");
  else minNumber_= 1 ;

}

HEEPElectronIDIso::~HEEPElectronIDIso(){}


bool HEEPElectronIDIso::filter(edm::Event & iEvent, const edm::EventSetup & iSetup){


  edm::Handle<edm::View<pat::Electron> > eleHandle;
  iEvent.getByLabel(eleLabel_,eleHandle);
  const edm::View<pat::Electron>& eles = *(eleHandle.product());

  int nTightElectron = 0;

  for(size_t eleNr = 0; eleNr != eles.size(); ++eleNr) {

   const pat::Electron& ele = eles[eleNr];
   int eleCutCode = ele.userInt("HEEPId");
   if(eleCutCode == 0) nTightElectron++ ;
  }

  if (nTightElectron >= minNumber_ && nTightElectron <= maxNumber_ ) return true;
  else return false;

}



// ------------ method called once each job just before starting event loop  ------------                                                                                                         
void HEEPElectronIDIso::beginJob(){}

// ------------ method called once each job just after ending the event loop  ------------                                                                                                        
void HEEPElectronIDIso::endJob(){}

// ------------ method called when starting to processes a run  ------------                                                                                                                     

bool HEEPElectronIDIso::beginRun(edm::Run const&, edm::EventSetup const&){ return true; }


// ------------ method called when ending the processing of a run  ------------                                                                                                                   
bool HEEPElectronIDIso::endRun(edm::Run const&, edm::EventSetup const&){ return true; }


bool HEEPElectronIDIso::beginLuminosityBlock (edm::LuminosityBlock const&, edm::EventSetup const&){ return true; }
bool HEEPElectronIDIso::endLuminosityBlock   (edm::LuminosityBlock const&, edm::EventSetup const&){ return true; }


void HEEPElectronIDIso::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation                                                                                                             
  // Please change this to state exactly what you do use, even if it is no parameters                                                                                                             
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//define this as a plug-in                                                                                                                                                     

DEFINE_FWK_MODULE(HEEPElectronIDIso);
