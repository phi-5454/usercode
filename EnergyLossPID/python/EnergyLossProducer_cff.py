import FWCore.ParameterSet.Config as cms

# Tracker local reco
from RecoLocalTracker.Configuration.RecoLocalTracker_cff import *

# Refitter
from RecoTracker.TrackProducer.TrackRefitters_cff import *

# Track refitter
from RecoTracker.TrackProducer.TrackRefitter_cfi import *
refitterForEnergyLoss     = TrackRefitter.clone()
refitterForEnergyLoss.src = 'generalTracks'

#refitterForEnergyLoss.MeasurementTrackerEvent = cms.InputTag('MeasurementTrackerEventPreSplitting')
refitterForEnergyLoss.MeasurementTrackerEvent = cms.InputTag('MeasurementTrackerEvent')

# Energy Loss
energyLossProducer   = cms.EDProducer("EnergyLossProducer",
  trackProducer      = cms.InputTag('refitterForEnergyLoss'),
  trajectoryProducer = cms.InputTag('refitterForEnergyLoss'),
  clusterShapes      = cms.InputTag('siPixelClusterShapeCache'),
  dedxProducer       = cms.InputTag('dedxHarmonic2'), # 'dedxTruncated40'
  tag = cms.string('totem')
)

# Micro Dst
microDstProducer  = cms.EDAnalyzer("MicroDstProducer",
  trackProducer      = cms.InputTag('refitterForEnergyLoss'),
  trajectoryProducer = cms.InputTag('refitterForEnergyLoss'),
  clusterShapes      = cms.InputTag('siPixelClusterShapeCache'),
  outFile            = cms.untracked.string('')
)


# Paths
produceEnergyLoss = cms.Path(refitterForEnergyLoss
                           * energyLossProducer)

produceMicroDst   = cms.Path(refitterForEnergyLoss
                           * microDstProducer)
