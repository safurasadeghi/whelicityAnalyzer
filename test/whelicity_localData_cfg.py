import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
## configure process options
process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary      = cms.untracked.bool(True)
)

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
'file:../../../../../DoubleMu/Run2016B/02477A4E-C586-E611-BC6F-02163E013D1C.root'
] );
secFiles.extend( [
               ] )


#genParticleCollection = 'prunedGenParticles'
#process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")
#process.initSubset.src = genParticleCollection
#process.decaySubset.src = genParticleCollection
#process.decaySubset.runMode = "Run2"
process.load("whelicity1.MiniAnalyzer.whelicity_cff")
process.Whelicity.isData = cms.bool(True)
process.TFileService = cms.Service("TFileService",
 fileName = cms.string("localDataHistos.root")
)


process.p = cms.Path(process.Whelicity)
