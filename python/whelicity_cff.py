import FWCore.ParameterSet.Config as cms


Whelicity = cms.EDAnalyzer("MiniAnalyzer",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    taus = cms.InputTag("slimmedTaus"),
    photons = cms.InputTag("slimmedPhotons"),
    jets = cms.InputTag("slimmedJets"),
    fatjets = cms.InputTag("slimmedJetsAK8"),
    mets = cms.InputTag("slimmedMETs"),
	packed = cms.InputTag("packedGenParticles"),
	pruned = cms.InputTag("prunedGenParticles"),
	genEvtInfo = cms.InputTag("generator"),
	lheEvtInfo = cms.InputTag("externalLHEProducer"),
	rho = cms.InputTag("fixedGridRhoFastjetAll"),
    effAreasConfigFile = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Spring15/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_25ns.txt"),
	ttgen = cms.InputTag("genEvt"),
	isData = cms.bool(False)
)


