import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:pat.root'
    )
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('preselectionStudy.root')
                                   )


process.demo = cms.EDAnalyzer('PreselectionStudy',
                              bdisc_name = cms.string('combinedSecondaryVertexBJetTags'),
                              RFijo = cms.untracked.bool(True),
                              outputLHCO = cms.string('pruebaPreselectionStudy'),
                              outputDiscriminants= cms.string('pruebaPreselectionStudy'),
                              postfix= cms.string('1')
                              )


process.p = cms.Path(process.demo)
