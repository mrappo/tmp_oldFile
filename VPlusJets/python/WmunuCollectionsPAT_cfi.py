import FWCore.ParameterSet.Config as cms

from ElectroWeakAnalysis.VPlusJets.AllPassFilter_cfi import AllPassFilter
from ElectroWeakAnalysis.VPlusJets.LooseLeptonVetoPAT_cfi import *

isMuonAnalyzer = True


def WmunuCollectionsPAT(process,isQCD, isHEEPID,isTransverseMassCut):

 isolationCutString = cms.string("")

 process.tightMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("selectedPatMuonsPFlow"),
    cut = cms.string("")
 )


 if isQCD:
    isolationCutString = "(pfIsolationR04().sumChargedHadronPt+max(0.,pfIsolationR04().sumNeutralHadronEt+pfIsolationR04().sumPhotonEt-0.5*pfIsolationR04().sumPUPt))/pt> 0.12" 
 else:
     if isHEEPID : isolationCutString = "trackIso()< 0.1"
     else : isolationCutString = "(pfIsolationR04().sumChargedHadronPt+max(0.,pfIsolationR04().sumNeutralHadronEt+pfIsolationR04().sumPhotonEt-0.5*pfIsolationR04().sumPUPt))/pt< 0.12"


 if isHEEPID : process.tightMuons.cut = cms.string(" isGlobalMuon && pt > 50 && abs(dB)<0.2 && globalTrack().hitPattern().numberOfValidPixelHits() >0 "
                                          " && globalTrack().hitPattern().numberOfValidMuonHits() >0 && globalTrack().hitPattern().trackerLayersWithMeasurement() > 8 "
                                          " && numberOfMatchedStations() > 1 && abs(eta)< 2.1 && "+ isolationCutString)
                     
 else : process.tightMuons.cut = cms.string(" pt>20 && isGlobalMuon && isPFMuon && abs(eta)<2.4 && globalTrack().normalizedChi2<10"
                                   " && globalTrack().hitPattern().numberOfValidMuonHits>0 && globalTrack().hitPattern().numberOfValidPixelHits>0 && numberOfMatchedStations>1"
                                   " && globalTrack().hitPattern().trackerLayersWithMeasurement>5 && " + isolationCutString)


 ## tight mu filter
 process.tightMuonFilter = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("tightMuons")                     
 )

 process.tightLeptonStep = AllPassFilter.clone()

 process.WToMunu = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tightMuons patMetShiftCorrected"),
    cut = cms.string(' daughter(0).pt >20 && daughter(1).pt >20 && sqrt(2*daughter(0).pt*daughter(1).pt*(1-cos(daughter(0).phi-daughter(1).phi)))>0'), 
    checkCharge = cms.bool(False),
 )

 if isTransverseMassCut : process.WToMunu.cut = cms.string(' daughter(1).pt >20  && sqrt(2*daughter(0).pt*daughter(1).pt*(1-cos(daughter(0).phi-daughter(1).phi)))>30')


 process.bestWmunu = cms.EDFilter("LargestPtCandViewSelector",
    maxNumber = cms.uint32(10),
    src = cms.InputTag("WToMunu")
 )

 process.bestWToLepnuStep = AllPassFilter.clone()
 
 ## --------- Loose Lepton Filters ----------

 LooseLeptonVetoPAT(process,isQCD, isHEEPID, isMuonAnalyzer)

 process.WSequence = cms.Sequence(process.tightMuons *
                         process.tightMuonFilter *
                         process.tightLeptonStep *
                         process.WToMunu *
                         process.bestWmunu *
                         process.bestWToLepnuStep
                         )

 process.WPath = cms.Sequence(process.WSequence*process.VetoSequence)
