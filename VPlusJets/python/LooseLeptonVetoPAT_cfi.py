import FWCore.ParameterSet.Config as cms
from ElectroWeakAnalysis.VPlusJets.AllPassFilter_cfi import AllPassFilter
from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *

def LooseLeptonVetoPAT(process,isQCD, isHEEPID, isMuonAnalyzer, looseEleIdLabel="loose"):

 if isHEEPID :

               process.looseMuons = cms.EDFilter("PATMuonRefSelector",
                                                  src = cms.InputTag("selectedPatMuons"),
                                                  cut = cms.string(""))

               process.looseMuons.cut = cms.string("isGlobalMuon && isTrackerMuon && pt>20 && abs(eta) < 2.4 && abs(phi)<3.2 && trackIso/pt < 0.1 && abs(dB) <0.2"
                                                   " && globalTrack().hitPattern().numberOfValidPixelHits>0 "
                                                   " && globalTrack().hitPattern().numberOfValidMuonHits>0 "
                                                   " && globalTrack().hitPattern().trackerLayersWithMeasurement>8 "
                                                   " && numberOfMatchedStations>1 " )
 
               process.looseMuonFilter = cms.EDFilter("PATCandViewCountFilter",
                                           src = cms.InputTag("looseMuons")
                                         )

               process.looseMuonStep = AllPassFilter.clone()
   
               if isMuonAnalyzer:
                    process.heepPatElectrons = cms.EDProducer("HEEPAttStatusToPAT",
                                                              eleLabel = cms.InputTag("selectedPatElectronsPFlow"),
                                                              barrelCuts = cms.PSet(heepBarrelCuts),
                                                              endcapCuts = cms.PSet(heepEndcapCuts),
                                                              applyRhoCorrToEleIsol = cms.bool(True),
                                                              eleIsolEffectiveAreas = cms.PSet (heepEffectiveAreas),
                                                              eleRhoCorrLabel = cms.InputTag("kt6PFJetsPFlow","rho"),
                                                              verticesLabel = cms.InputTag("offlinePrimaryVerticesWithBS")
                                                              )
               
               process.looseElectronFilter = cms.EDFilter("HEEPElectronIDIso",
                                                          electronCollection = cms.InputTag("heepPatElectrons"),
                                                          eleIdType = cms.string("LooseID"),
                                                          maxNumber = cms.untracked.int32(1),
                                                          minNumber = cms.untracked.int32(1)                                                          
                                                          )

               process.looseElectrons = cms.EDProducer("HEEPElectronProducer",
                                                       electronCollection = cms.InputTag("heepPatElectrons"),
                                                       eleIdType = cms.string("LooseID"),
                                                     )

               if isMuonAnalyzer :
                                  process.looseMuonFilter.minNumber = cms.uint32(1)
                                  process.looseMuonFilter.maxNumber = cms.uint32(1)
                                  process.looseElectronFilter.minNumber = cms.untracked.int32(0)
                                  process.looseElectronFilter.maxNumber = cms.untracked.int32(0)
               else :
                                  process.looseMuonFilter.minNumber = cms.uint32(0)
                                  process.looseMuonFilter.maxNumber = cms.uint32(0)
                                  process.looseElectronFilter.minNumber = cms.untracked.int32(1)
                                  process.looseElectronFilter.maxNumber = cms.untracked.int32(1)
                      
               process.looseElectronStep = AllPassFilter.clone()

               if isMuonAnalyzer:
                  process.VetoSequence = cms.Sequence(process.looseMuons*process.looseMuonFilter*process.looseMuonStep*process.heepPatElectrons*process.looseElectrons*process.looseElectronFilter*process.looseElectronStep)
               else: process.VetoSequence = cms.Sequence(process.looseElectrons*process.looseElectronFilter*process.looseElectronStep*process.looseMuons*process.looseMuonFilter*process.looseMuonStep)
                                                
 else:

         ##  Define loose electron selection for veto ######

         process.looseElectrons = cms.EDProducer("PATElectronIdSelector",
                                                  src = cms.InputTag( "selectedPatElectronsPFlow" ),
                                                  idLabel = cms.string(looseEleIdLabel),
                                                  useMVAbasedID_   = cms.bool(True),
                                                  idType  = cms.string("HiggsID")  
                                                )

         ##  Define loose muon selection for veto ######

         process.looseMuons = cms.EDFilter("PATMuonRefSelector",
                                   src = cms.InputTag("selectedPatMuonsPFlow"),
                                   cut = cms.string(""))

         process.looseMuons.cut = cms.string(" pt>10 &&isPFMuon && (isGlobalMuon || isTrackerMuon) && abs(eta)<2.4"
                                             " && (pfIsolationR04().sumChargedHadronPt+max(0.,pfIsolationR04().sumNeutralHadronEt+pfIsolationR04().sumPhotonEt-0.5*pfIsolationR04().sumPUPt))/pt< 0.2")

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
   
         process.VetoSequence = cms.Sequence(process.looseElectrons*process.looseElectronFilter*process.looseElectronStep*process.looseMuons*process.looseMuonFilter*process.looseMuonStep)
