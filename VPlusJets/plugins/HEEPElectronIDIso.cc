#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

//trigger                                                                                                                                                                            
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
  std::string   IdLabel_;
  bool          applyLooseID_;
  bool          applyTightID_;
  int           maxNumber_;
  int           minNumber_;

};

HEEPElectronIDIso::HEEPElectronIDIso(const edm::ParameterSet &  iConfig){

  if( iConfig.existsAs<edm::InputTag>("electronCollection") )
    eleLabel_=iConfig.getParameter< edm::InputTag >("electronCollection");
  else eleLabel_= edm::InputTag("heepPatElectrons");

  if( iConfig.existsAs<std::string>("eleIdType") )
    IdLabel_=iConfig.getParameter<std::string>("eleIdType");
  else IdLabel_= std::string("TightID");

  maxNumber_=iConfig.getUntrackedParameter<int>("maxNumber", 0);
  minNumber_=iConfig.getUntrackedParameter<int>("minNumber", 999);


  applyTightID_ = false ;
  applyLooseID_ = false ;

  if( IdLabel_ == "TightID" || IdLabel_ == "Tight" || IdLabel_ == "tightID" || IdLabel_ == "tightId" || IdLabel_ == "tightid" || IdLabel_ == "tight" ) applyTightID_ = true ;
  if( IdLabel_ == "LooseID" || IdLabel_ == "Loose" || IdLabel_ == "LooseID" || IdLabel_ == "looseId" || IdLabel_ == "looseid" || IdLabel_ == "loose" ) applyLooseID_ = true ;

  

}

HEEPElectronIDIso::~HEEPElectronIDIso(){}


bool HEEPElectronIDIso::filter(edm::Event & iEvent, const edm::EventSetup & iSetup){


  edm::Handle<edm::View<pat::Electron> > eleHandle;
  iEvent.getByLabel(eleLabel_,eleHandle);
  const edm::View<pat::Electron>& eles = *(eleHandle.product());

  Int_t nTightElectron = 0;
  Int_t nLooseElectron = 0;

  for(size_t eleNr = 0; eleNr != eles.size(); ++eleNr) {

   const pat::Electron& ele = eles[eleNr];
 
   int eleCutCode = ele.userInt("HEEPId");
   Double_t et = ele.caloEnergy()*sin(ele.p4().theta());
   Double_t scEta = fabs(ele.superCluster()->eta());

   if(eleCutCode == 0 && et > 90 && fabs(ele.eta()) < 2.5 && !(scEta > 1.4442 && scEta < 1.566) && fabs(ele.phi()) < 3.2 ) nTightElectron++ ;
   if(eleCutCode == 0 && et > 20 && fabs(ele.eta()) < 2.5 && !(scEta > 1.4442 && scEta < 1.566) && fabs(ele.phi()) < 3.2 ) nLooseElectron++ ;

  }

  if(applyTightID_){
                    if(nTightElectron >= minNumber_ && nTightElectron <= maxNumber_ ) return true;
                    else return false ;
  }
  
  if(applyLooseID_){
   if(nLooseElectron >= minNumber_ && nLooseElectron <= maxNumber_ ) return true;
   else return false ;
  }

  return false;

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
