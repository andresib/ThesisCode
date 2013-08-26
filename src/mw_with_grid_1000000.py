import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:pat.root'
    )
)

process.demo = cms.EDAnalyzer('PruebaMWWithGrid',
                              s_lhco_source_address = cms.string("LHCO_data2/"),
                              b_runningLocally = cms.bool(False),
                              maxEvents = cms.int32(1000000)
)


process.p = cms.Path(process.demo)
