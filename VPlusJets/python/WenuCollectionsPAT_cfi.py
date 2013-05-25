import FWCore.ParameterSet.Config as cms

from ElectroWeakAnalysis.VPlusJets.AllPassFilter_cfi import AllPassFilter
from ElectroWeakAnalysis.VPlusJets.LooseLeptonVetoPAT_cfi import *
from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *
 
isMuonAnalyzer = False


def WenuCollectionsPAT(process,isQCD,isHEEPID,isTransverseMassCut) :

#WP80 electrons, only track iso, remove H/E cut

 tightEleIdLabel = "tight"
 looseEleIdLabel = "loose"

 if isQCD:
  tightEleIdLabel = "qcd"
  looseEleIdLabel = "qcd"


 ## modified WP70
 if isHEEPID:

      process.heepPatElectrons = cms.EDProducer("HEEPAttStatusToPAT",
                                               eleLabel = cms.InputTag("selectedPatElectronsPFlow"),
                                               barrelCuts = cms.PSet(heepBarrelCuts),
                                               endcapCuts = cms.PSet(heepEndcapCuts),
                                               applyRhoCorrToEleIsol = cms.bool(True),
                                               eleIsolEffectiveAreas = cms.PSet (heepEffectiveAreas),
                                               eleRhoCorrLabel = cms.InputTag("kt6PFJetsPFlow","rho"),
                                               verticesLabel = cms.InputTag("offlinePrimaryVerticesWithBS")
                                               )  

      process.tightElectrons  = cms.EDProducer("HEEPElectronProducer",
                                               electronCollection = cms.InputTag("heepPatElectrons"),
                                               eleIdType = cms.string("TightID"),
                                              )

      process.tightElectronFilter = cms.EDFilter("PATCandViewCountFilter",
                                                  minNumber = cms.uint32(1),
                                                  maxNumber = cms.uint32(999999),
                                                  src = cms.InputTag("tightElectrons")
                                                )
          
      process.tightLeptonStep = AllPassFilter.clone()

      process.WToEnu = cms.EDProducer("CandViewShallowCloneCombiner",
                                  decay = cms.string("tightElectrons patMetShiftCorrected"),
                                  cut = cms.string('daughter(0).pt >20 && daughter(1).pt >20  && sqrt(2*daughter(0).pt*daughter(1).pt*(1-cos(daughter(0).phi-daughter(1).phi)))>0'),
                                  checkCharge = cms.bool(False),
      )

      if isTransverseMassCut : process.WToEnu.cut = cms.string('daughter(1).pt >20  && sqrt(2*daughter(0).pt*daughter(1).pt*(1-cos(daughter(0).phi-daughter(1).phi)))>30')


      process.bestWToEnu =cms.EDFilter("LargestPtCandViewSelector",
                                   maxNumber = cms.uint32(10),
                                   src = cms.InputTag("WToEnu")                 
      )
      process.bestWToLepnuStep = AllPassFilter.clone()


      process.WSequence = cms.Sequence(process.heepPatElectrons*
                          process.tightElectrons *
                          process.tightElectronFilter *
                          process.tightLeptonStep *
                          process.WToEnu *
                          process.bestWToEnu *
                          process.bestWToLepnuStep
      )

      LooseLeptonVetoPAT(process,isQCD, isHEEPID, isMuonAnalyzer, looseEleIdLabel)

      process.WPath = cms.Sequence(process.WSequence*process.VetoSequence)
      


 else:
      process.tightElectrons = cms.EDProducer("PATElectronIdSelector",
                                               src = cms.InputTag( "selectedPatElectronsPFlow" ),
                                               idLabel = cms.string(tightEleIdLabel),  # refers to Higgs Cut Based or MVA WP
                                               useMVAbasedID_   = cms.bool(True),
                                               idType = cms.string("HiggsID")
      )

      process.tightElectronFilter = cms.EDFilter("PATCandViewCountFilter",
        minNumber = cms.uint32(1),
        maxNumber = cms.uint32(999999),
        src = cms.InputTag("tightElectrons")                     
      )

      ## tight ele filter

      process.tightLeptonStep = AllPassFilter.clone()


      process.WToEnu = cms.EDProducer("CandViewShallowCloneCombiner",
                                  decay = cms.string("tightElectrons patMetShiftCorrected"),
                                  cut = cms.string('daughter(1).pt >20  && sqrt(2*daughter(0).pt*daughter(1).pt*(1-cos(daughter(0).phi-daughter(1).phi)))>0'),
                                  checkCharge = cms.bool(False),
      )

      if isTransverseMassCut : process.WToEnu.cut = cms.string('daughter(1).pt >20  && sqrt(2*daughter(0).pt*daughter(1).pt*(1-cos(daughter(0).phi-daughter(1).phi)))>30')


      process.bestWToEnu =cms.EDFilter("LargestPtCandViewSelector",
                                   maxNumber = cms.uint32(10),
                                   src = cms.InputTag("WToEnu")                 
      )
      process.bestWToLepnuStep = AllPassFilter.clone()


      process.WSequence = cms.Sequence(process.tightElectrons *
                         process.tightElectronFilter *
                         process.tightLeptonStep *
                         process.WToEnu *
                         process.bestWToEnu *
                         process.bestWToLepnuStep
      )

      LooseLeptonVetoPAT(process,isQCD, isHEEPID, isMuonAnalyzer, looseEleIdLabel)

      process.WPath = cms.Sequence(process.WSequence*process.VetoSequence)
