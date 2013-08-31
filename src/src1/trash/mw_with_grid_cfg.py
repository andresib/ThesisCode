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

process.demo = cms.EDAnalyzer('MW_With_Grid',
                              s_lhco_source_address = cms.string("LHCO/gen/background/semileptonic_tt/"),
                              b_runningLocally = cms.bool(False)
)


process.p = cms.Path(process.demo)
