import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:pat.root'
    )
)
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('top.root')
                                   )

process.demo = cms.EDAnalyzer('FullAndWorking',
                              top = cms.untracked.bool(True),
                              neutralino = cms.untracked.bool(False),
                              RFijo = cms.untracked.bool(True),
                              LHCOWithMC = cms.untracked.bool(True),
                              LHCOWithData = cms.untracked.bool(True),
                              MC = cms.untracked.bool(True),
                              bdisc_name = cms.string('combinedSecondaryVertexBJetTags'),
                              selector = cms.untracked.double(9),
                              c1 = cms.untracked.double(2),
                              c2 = cms.untracked.double(5),
                              c3 = cms.untracked.double(10000),
                              c4 = cms.untracked.double(10),
                              c5 = cms.untracked.double(0.1),
                              c6 = cms.untracked.double(0.1),
                              c7 = cms.untracked.double(100),
                              c8 = cms.untracked.double(6),
                              c9 = cms.untracked.double(6),
                              c10 = cms.untracked.double(0.1),
                              c11 = cms.untracked.double(0.1),
                              c12 =cms.untracked.double(10000),
                              c13 =cms.untracked.double(4.5),
                              c14 =cms.untracked.double(4.5)
                              
                              )


process.p = cms.Path(process.demo)
