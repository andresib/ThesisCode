import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:-input-'
    )
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('-output-/preselectionStudy_-postfix-')
                                   )


process.demo = cms.EDAnalyzer('PreselectionStudy',
                              bdisc_name = cms.string('combinedSecondaryVertexBJetTags'),
                              outputLHCO = cms.string('-outputLHCO-'),
                              outputDiscriminants = cms.string('-output-'),
                              postfix = cms.string('-postfix-'),
                              RFijo = cms.untracked.bool(True)
                              )


process.p = cms.Path(process.demo)
