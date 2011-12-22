import FWCore.ParameterSet.Config as cms

from PhysicsTools.MiBiCommonPAT.makeMiBiCommonNT_cff import *

process = cms.Process("MiBiCommonNT")
process.source = cms.Source("PoolSource",
     fileNames = cms.untracked.vstring('file:/tmp/govoni/reducedSample.root')
)


process.load("PhysicsTools.MiBiCommonPAT.SimpleNtupleGen_cfi") 

process.MiBiCommonPAT = cms.Sequence(process.SimpleNtupleGen)

process.MiBiPathPFlow = cms.Path(process.MiBiCommonPAT) 
			 
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.TFileService = cms.Service(
    "TFileService", 
    fileName = cms.string("MiBiCommonNT.root"),
    closeFileFast = cms.untracked.bool(True)
)
