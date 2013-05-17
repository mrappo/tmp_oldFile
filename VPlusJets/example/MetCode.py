process.load("JetMETCorrections/Type1MET/pfMETsysShiftCorrections_cfi")
process.pfMEtSysShiftCorr.src = cms.InputTag('patMETsPFlow')
process.pfMEtSysShiftCorr.srcMEt = cms.InputTag('patMETsPFlow')
process.pfMEtSysShiftCorr.srcJets = cms.InputTag('selectedPatJetsPFlow')

if isMC:
 # pfMEtSysShiftCorrParameters_2012runAplusBvsNvtx_mc
 process.pfMEtSysShiftCorr.parameter = cms.PSet(
      numJetsMin = cms.int32(-1),
      numJetsMax = cms.int32(-1),
      px = cms.string("+0.1166 + 0.0200*Nvtx"),
      py = cms.string("+0.2764 - 0.1280*Nvtx")
 #    px = cms.string("+2.22335e-02 - 6.59183e-02*Nvtx"),
 #    py = cms.string("+1.52720e-01 - 1.28052e-01*Nvtx")
 )
 else:
 # pfMEtSysShiftCorrParameters_2012runAplusBvsNvtx_data
 process.pfMEtSysShiftCorr.parameter = cms.PSet(
      numJetsMin = cms.int32(-1),
      numJetsMax = cms.int32(-1),
      px = cms.string("+0.2661 + 0.3217*Nvtx"),
      py = cms.string("-0.2251 - 0.1747*Nvtx")
 #    px = cms.string("+1.68804e-01 + 3.37139e-01*Nvtx"),
 #    py = cms.string("-1.72555e-01 - 1.79594e-01*Nvtx")
 )
            
process.patMetShiftCorrected = cms.EDProducer("CorrectedPATMETProducer",
                                               src = cms.InputTag('patMETsPFlow'),
                                               applyType1Corrections = cms.bool(True),
                                               srcType1Corrections = cms.VInputTag(
                                               cms.InputTag('pfMEtSysShiftCorr')),
                                               applyType2Corrections = cms.bool(False)
                                              )
#### This has to be run in this way
process.myseq = cms.Sequence(
        process.pfMEtSysShiftCorrSequence *
        process.patMetShiftCorrected
        )
