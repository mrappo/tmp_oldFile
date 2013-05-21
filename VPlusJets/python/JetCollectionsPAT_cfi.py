import FWCore.ParameterSet.Config as cms
from ElectroWeakAnalysis.VPlusJets.AllPassFilter_cfi import AllPassFilter
from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
from RecoJets.Configuration.GenJetParticles_cff import *
from RecoJets.JetProducers.ak5GenJets_cfi import *
from RecoMET.Configuration.GenMETParticles_cff import *
from RecoMET.METProducers.genMetTrue_cfi import *

##########################################################################

def JetCollectionsPAT(process,isHEEPID):
 
  if isHEEPID:
             # Apply loose PF jet ID
             process.ak5PFGoodJets = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                                   filterParams = pfJetIDSelector.clone(),
                                                   src = cms.InputTag("selectedPatJetsPFlow"),
                                                   filter = cms.bool(True))


  else : 

       # Apply loose PileUp PF jet ID
       process.ak5PFnoPUJets = cms.EDProducer("PATPuJetIdSelector",
                                               src = cms.InputTag( "selectedPatJetsPFlow" ),
                                               idLabel = cms.string("loose"),
                                               valueMapLabel = cms.string("puJetMvaChs"))
       # Apply loose PF jet ID
       process.ak5PFGoodJets = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                             filterParams = pfJetIDSelector.clone(),
                                             src = cms.InputTag("ak5PFnoPUJets"),
                                             filter = cms.bool(True))


  ##-------- Remove electrons and muons from jet collection ----------------------
  process.ak5PFJetsClean = cms.EDProducer("PFPATJetCleaner",
                                           srcJets = cms.InputTag("ak5PFGoodJets"),
                                           module_label = cms.string(""),
                                           srcObjects = cms.VInputTag(cms.InputTag("looseElectrons"),cms.InputTag("looseMuons")),
                                           deltaRMin = cms.double(0.3))


  process.ak5PFJetsLooseIdAll = cms.EDFilter("PATJetRefSelector",
                                              src = cms.InputTag("ak5PFJetsClean"),
                                              cut = cms.string('pt > 30.0'))


  process.ak5cleaningPTCut = AllPassFilter.clone()

  process.ak5PFJetsLooseId = cms.EDFilter("PATJetRefSelector",
                                           src = cms.InputTag("ak5PFJetsClean"),
                                           cut = cms.string('pt > 30.0 && abs(eta) < 2.4'))

  process.ak5PFJetsLooseIdVBFTag = cms.EDFilter("PATJetRefSelector",
                                                src = cms.InputTag("ak5PFJetsClean"),
                                                cut = cms.string('pt > 30.0 && abs(eta) > 2.4 && abs(eta) < 9.9'))

  ##########################################
  ## Filter to require at least two jets in the event

  process.RequireTwoJets = cms.EDFilter("PATCandViewCountFilter",
                                         minNumber = cms.uint32(2),
                                         maxNumber = cms.uint32(100),
                                         src = cms.InputTag("ak5PFJetsLooseId"),)


  ############################################
  if isHEEPID: process.PFJetPath = cms.Sequence( process.ak5PFGoodJets + process.ak5PFJetsClean + process.ak5PFJetsLooseId +  process.ak5PFJetsLooseIdAll + process.ak5cleaningPTCut +
                                    process.ak5PFJetsLooseIdVBFTag )

  else : process.PFJetPath = cms.Sequence( process.ak5PFnoPUJets +process.ak5PFGoodJets + process.ak5PFJetsClean + process.ak5PFJetsLooseId +  process.ak5PFJetsLooseIdAll +
                                           process.ak5cleaningPTCut + process.ak5PFJetsLooseIdVBFTag )

  ########################################################################
  #############  Jets in Monte Carlo  #############
  ##########################################################################


  ##################### Tag jets: Needed for MC flavor matching
  process.myPartons = cms.EDProducer("PartonSelector",
    src = cms.InputTag("genParticles"),
    withLeptons = cms.bool(False))

  process.ak5flavourByRef = cms.EDProducer("JetPartonMatcher",
    jets = cms.InputTag("selectedPatJetsPFlow"),
    coneSizeToAssociate = cms.double(0.3),
    partons = cms.InputTag("myPartons"))

  process.ak5tagJet = cms.EDProducer("JetFlavourIdentifier",
    srcByReference = cms.InputTag("ak5flavourByRef"),
    physicsDefinition = cms.bool(False))

  process.TagJetPath = cms.Sequence( process.myPartons+process.ak5flavourByRef*process.ak5tagJet) 

  process.GenJetPath = cms.Sequence( process.genParticlesForJets + process.ak5GenJets + process.genParticlesForMETAllVisible + process.genMetTrue)

#############################################

