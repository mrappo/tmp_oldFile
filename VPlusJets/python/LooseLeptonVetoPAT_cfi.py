import FWCore.ParameterSet.Config as cms
from ElectroWeakAnalysis.VPlusJets.AllPassFilter_cfi import AllPassFilter

def LooseLeptonVetoPAT(process,isQCD, isHEEPID, isMuonAnalyzer, looseEleIdLabel="loose"):

 ##  Define loose electron selection for veto ######

 process.looseElectrons = cms.EDProducer("PATElectronIdSelector",
                                         src = cms.InputTag( "selectedPatElectronsPFlow" ),
                                         idLabel = cms.string(looseEleIdLabel),
                                         useMVAbasedID_   = cms.bool(True),
                                         idType  = cms.string("")  
 )

 process.looseElectrons.idType = cms.string("HiggsID")
 if isHEEPID : process.looseElectrons.idType = cms.string("HEEPID")


 ##  Define loose muon selection for veto ######

 process.looseMuons = cms.EDFilter("PATMuonRefSelector",
                                   src = cms.InputTag("selectedPatMuonsPFlow"),
                                   cut = cms.string(""),
                          )

 process.looseMuons.cut = cms.string(" pt>10 &&isPFMuon && (isGlobalMuon || isTrackerMuon) && abs(eta)<2.4"
                                     " && (pfIsolationR04().sumChargedHadronPt+max(0.,pfIsolationR04().sumNeutralHadronEt+pfIsolationR04().sumPhotonEt-0.5*pfIsolationR04().sumPUPt))/pt< 0.2")

 if isHEEPID : process.looseMuons.cut = cms.string("pt>20 && abs(eta) < 2.1 && trackIso/pt < 0.1 && abs(dB) <0.2"
                                          " && globalTrack().hitPattern().numberOfValidPixelHits>0 "
                                          " && globalTrack().hitPattern().numberOfValidMuonHits>0 "
                                          " && globalTrack().hitPattern().trackerLayersWithMeasurement>8 "
                                          " && numberOfMatchedStations>1 " )

 if isMuonAnalyzer :
      nLooseElectron = 0 ;
      nLooseMuon     = 1 ;
 else :
      nLooseElectron = 1 ;
      nLooseMuon     = 0 ;
 
 process.looseElectronFilter = cms.EDFilter("PATCandViewCountFilter",
                                       minNumber = cms.uint32(nLooseElectron),
                                       maxNumber = cms.uint32(nLooseElectron),
                                       src = cms.InputTag("looseElectrons")
                                    )

 process.looseElectronStep = AllPassFilter.clone()


 process.looseMuonFilter = cms.EDFilter("PATCandViewCountFilter",
                                        minNumber = cms.uint32(nLooseMuon),
                                        maxNumber = cms.uint32(nLooseMuon),
                                        src = cms.InputTag("looseMuons")
                                  )

 process.looseMuonStep = AllPassFilter.clone()
   
 process.VetoSequence = cms.Sequence( process.looseElectrons *
                              process.looseElectronFilter *
                              process.looseElectronStep *
                              process.looseMuons *
                              process.looseMuonFilter *
                              process.looseMuonStep
                            )
