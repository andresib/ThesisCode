#!/bin/tcsh

pushd $CMSSW_BASE/src

cvs co -r V06-05-01    DataFormats/PatCandidates
cvs co -r V08-09-11-02 PhysicsTools/PatAlgos
cvs co -r V03-09-22    PhysicsTools/PatUtils
cvs co -r V00-03-14    CommonTools/ParticleFlow
cvs co -r V00-00-09    CommonTools/RecoUtils
cvs co -r V04-06-09    JetMETCorrections/Type1MET
cvs co -r V00-00-08    RecoMET/METAnalyzers
cvs co -r V00-00-07    RecoMET/METFilters
cvs co -r V00-00-13 -d EGamma/EGammaAnalysisTools UserCode/EGamma/EGammaAnalysisTools
cd EGamma/EGammaAnalysisTools/data
cat download.url | xargs wget
cd -

# checkdeps -a ; rm -r ElectroWeakAnalysis TauAnalysis AnalysisDataFormats TopQuarkAnalysis
# just results in these three:
addpkg PhysicsTools/PatExamples
addpkg PhysicsTools/SelectorUtils
addpkg PhysicsTools/TagAndProbe

scram b -j 4
popd
