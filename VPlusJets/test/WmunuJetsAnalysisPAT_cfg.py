import FWCore.ParameterSet.Config as cms
import pprint


### set some important flags
isMC = False
isQCD = False
isHEEPID = True
isTransverseMassCut = False

process = cms.Process("demo")

##---------  Load standard Reco modules ------------
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')



##----- this config frament brings you the generator information ----
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
process.load("Configuration.StandardSequences.Generator_cff")


##----- Detector geometry : some of these needed for b-tag -------
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")


##----- B-tags --------------
process.load("RecoBTag.Configuration.RecoBTag_cff")


##----- Global tag: conditions database ------------
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

##----- Counter module ------------
process.load("ElectroWeakAnalysis.VPlusJets.AllPassFilter_cfi")

## import skeleton process
#from PhysicsTools.PatAlgos.patTemplate_cfg import *


############################################
if not isMC:
    process.GlobalTag.globaltag = 'GR_P_V39_AN3::All'
else:
    process.GlobalTag.globaltag = 'START53_V15::All'

############################################
########################################################################################
########################################################################################

##---------  W-->munu Collection ------------
from ElectroWeakAnalysis.VPlusJets.WmunuCollectionsPAT_cfi import*

WmunuCollectionsPAT(process,isQCD,isHEEPID,isTransverseMassCut)

##---------  Jet Collection ----------------
process.load("ElectroWeakAnalysis.VPlusJets.JetCollectionsPAT_cfi")

##---------  Vertex and track Collections -----------
process.load("ElectroWeakAnalysis.VPlusJets.TrackCollections_cfi")
#

numEventsToRun = -1

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(numEventsToRun)
)

process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

OutputFileName = "WmunuJetAnalysisntuple_gr2_207231_920_1260547040.root"


process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP1_190703_11_8519680.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP1_191834_99_120189151.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP1_194199_258_229967752.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP1_194699_80_97577634.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP1_194712_559_405404310.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP1_194778_190_248720484.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP1_194912_390_657997761.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP1_195013_303_44790264.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP1_195099_171_245228003.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP1_195398_899_752941172.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP1_195399_159_130828536.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP1_195930_83_62880579.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP1_196438_586_479961054.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP1_196452_193_232654780.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP1_198522_104_75971527.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP1_199428_457_557329806.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP1_199428_86_80090879.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP1_199961_177_187360166.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP1_200188_104_142020231.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP1_203002_1517_1719001770.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP1_204544_102_152417542.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP1_204564_359_400306553.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP1_207372_380_547252863.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP1_207490_80_87321589.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP1_208487_118_209515269.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP2_194199_190_149891910.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP2_194897_50_91029885.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP2_195397_558_772053348.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP2_199608_1123_1247460423.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP2_200091_1678_1762413892.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP2_202299_249_336447785.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP2_202314_128_132710041.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP2_206476_128_148892985.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP2_206859_764_1094226833.root',
        #'file:/data2/rgerosa/SideBandClosureTest/PAT53X_Sideband/pat_53X_ofile_mu_GROUP2_207231_920_1260547040.root',
        #'/store/user/lnujj/PatTuples_8TeV_53X/custodio/BulkG_WW_lvjj_c0p2_M1500-JHU-v1/SQWaT_PAT_53X_Summer12_v1/829f288d768dd564418efaaf3a8ab9aa/pat_53x_test_v03_4_1_jvg.root'
) )




##-------- Electron events of interest --------
process.HLTMu =cms.EDFilter("HLTHighLevel",
     TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
     HLTPaths = cms.vstring('HLT_IsoMu24_*','HLT_IsoMu30_*'),
     eventSetupPathsKey = cms.string(''),
     andOr = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
     throw = cms.bool(False) # throw exception on unknown path names
)

process.primaryVertex.src = cms.InputTag("goodOfflinePrimaryVertices");
process.primaryVertex.cut = cms.string(" ");

##########################################
## Filter to require at least two jets in the event
process.RequireTwoJetsORboostedV = cms.EDFilter("JetsORboostedV",
    minNumber = cms.untracked.int32(2),
    maxNumber = cms.untracked.int32(100),
    srcJets = cms.InputTag("ak5PFJetsLooseIdAll"),
    srcVectorBoson = cms.InputTag("bestWmunu"),
    minVpt = cms.untracked.double(100.),
    minNumberPhotons = cms.untracked.int32(0)
)
process.RequireTwoJetsORboostedVStep = process.AllPassFilter.clone()


##-------- Save V+jets trees --------
process.VplusJets = cms.EDAnalyzer("VplusJetsAnalysis", 
    jetType = cms.string("PF"),
    srcPFCor = cms.InputTag("ak5PFJetsLooseId"),
    srcPhoton = cms.InputTag("photons"),
    IsoValPhoton = cms.VInputTag(cms.InputTag('phoPFIso:chIsoForGsfEle'),
                                 cms.InputTag('phoPFIso:phIsoForGsfEle'),
                                 cms.InputTag('phoPFIso:nhIsoForGsfEle'),
                                                           ),
    srcPFCorVBFTag = cms.InputTag("ak5PFJetsLooseIdVBFTag"), 
    srcVectorBoson = cms.InputTag("bestWmunu"),
    VBosonType     = cms.string('W'),
    LeptonType     = cms.string('muon'),                          
    TreeName    = cms.string('WJet'),
    srcPrimaryVertex = cms.InputTag("goodOfflinePrimaryVertices"),                               
    runningOverMC = cms.bool(isMC),			
    runningOverAOD = cms.bool(False),			
    srcMet = cms.InputTag("patMetShiftCorrected"), # type1 + shift corrections
    srcMetMVA = cms.InputTag("pfMEtMVA"),
    srcGen  = cms.InputTag("ak5GenJets"),
    srcMuons  = cms.InputTag("selectedPatMuonsPFlow"),
    srcBeamSpot  = cms.InputTag("offlineBeamSpot"),
    srcCaloMet  = cms.InputTag("patMETs"),
    srcgenMet  = cms.InputTag("genMetTrue"),
    srcGenParticles  = cms.InputTag("genParticles"),
    srcTcMet    = cms.InputTag("patMETsAK5TC"),
    srcJetsforRho = cms.string("kt6PFJetsPFlow"),                               
    srcJetsforRho_lepIso = cms.string("kt6PFJetsForIsolation"),       
    srcJetsforRhoCHS = cms.string("kt6PFJetsChsPFlow"),
    srcJetsforRho_lepIsoCHS = cms.string("kt6PFJetsChsForIsolationPFlow"),
    srcFlavorByValue = cms.InputTag("ak5tagJet"),
    bTagger=cms.string("simpleSecondaryVertexHighEffBJetTags"),

    applyJECToGroomedJets=cms.bool(True),
    doGroomedAK5 = cms.bool(True),
    doGroomedAK7 = cms.bool(True),
    doGroomedAK8 = cms.bool(False),
    doGroomedCA8 = cms.bool(True),
    doGroomedCA12 = cms.bool(False)
)

if isMC:
    process.VplusJets.JEC_GlobalTag_forGroomedJet = cms.string("START53_V15")
else:
    process.VplusJets.JEC_GlobalTag_forGroomedJet = cms.string("GR_P_V39_AN3")



process.TFileService = cms.Service(
    "TFileService", 
    fileName = cms.string( OutputFileName ),
    closeFileFast = cms.untracked.bool(False)
)


##
## MET shift correction on the fly
##
#from PhysicsTools.PatAlgos.tools.metTools import addPfMET

process.load("JetMETCorrections/Type1MET/pfMETsysShiftCorrections_cfi")
process.pfMEtSysShiftCorr.src = cms.InputTag('patMETsPFlow')
process.pfMEtSysShiftCorr.srcMEt = cms.InputTag('patMETsPFlow')
process.pfMEtSysShiftCorr.srcJets = cms.InputTag('selectedPatJetsPFlow')

if isMC:
# pfMEtSysShiftCorrParameters_2012runAplusBvsNvtx_mc
  process.pfMEtSysShiftCorr.parameter = cms.PSet(
    numJetsMin = cms.int32(-1),
    numJetsMax = cms.int32(-1),
    px = cms.string("+2.22335e-02 - 6.59183e-02*Nvtx"),
    py = cms.string("+1.52720e-01 - 1.28052e-01*Nvtx")
  )
else:
# pfMEtSysShiftCorrParameters_2012runAplusBvsNvtx_data
  process.pfMEtSysShiftCorr.parameter = cms.PSet(
    numJetsMin = cms.int32(-1),
    numJetsMax = cms.int32(-1),
    px = cms.string("+1.68804e-01 + 3.37139e-01*Nvtx"),
    py = cms.string("-1.72555e-01 - 1.79594e-01*Nvtx")
  )

process.patMetShiftCorrected = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag('patMETsPFlow'),
    applyType1Corrections = cms.bool(True),
    srcType1Corrections = cms.VInputTag(
        cms.InputTag('pfMEtSysShiftCorr')
    ),
    applyType2Corrections = cms.bool(False)
)   

process.myseq = cms.Sequence(
    process.pfMEtSysShiftCorrSequence *
    process.patMetShiftCorrected *
    process.TrackVtxPath *
    process.HLTMu *
    process.WPath *
    process.GenJetPath *
##    process.btagging * 
    process.TagJetPath *
    process.PFJetPath *
    process.RequireTwoJetsORboostedV *
    process.RequireTwoJetsORboostedVStep
    )

if isMC:
    process.myseq.remove ( process.noscraping)
    process.myseq.remove ( process.HBHENoiseFilter)
    process.myseq.remove ( process.HLTMu)
else:
    process.myseq.remove ( process.noscraping)
    process.myseq.remove ( process.HBHENoiseFilter)
    process.myseq.remove ( process.GenJetPath)
    process.myseq.remove ( process.TagJetPath)


##---- if do not want to require >= 2 jets then disable that filter ---
##process.myseq.remove ( process.RequireTwoJets)  

#process.outpath.remove(process.out)
process.p = cms.Path( process.myseq * process.VplusJets)
