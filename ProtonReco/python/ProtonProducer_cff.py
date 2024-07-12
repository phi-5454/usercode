import FWCore.ParameterSet.Config as cms

# Energy Loss
protonProducer   = cms.EDProducer("ProtonProducer",
  totemRpClusters    = cms.InputTag('totemRPClusterProducer', '',  'RECO'),
  localTracks        = cms.InputTag('totemRPLocalTrackFitter','','reRECO'),
  #outFile = cms.string('out.dat.gz')
  outFile = cms.string('output.root')
)

# Paths
produceProtons = cms.Path(protonProducer)

