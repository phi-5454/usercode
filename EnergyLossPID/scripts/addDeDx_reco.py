import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process("addDeDx_reco", eras.Run2_2018_highBetaStar)

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("UserCode.EnergyLossPID.EnergyLossProducer_cff")

###############################################################################
# Message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.MessageLogger = cms.Service("MessageLogger",
  cerr = cms.untracked.PSet(
    threshold = cms.untracked.string('DEBUG'),
    DEBUG     = cms.untracked.PSet(
      limit = cms.untracked.int32(0)
    ),
    FwkReport = cms.untracked.PSet(
      optionalPSet = cms.untracked.bool(True),
      reportEvery = cms.untracked.int32(50),
      limit = cms.untracked.int32(1000000)
    )
  ),
  destinations = cms.untracked.vstring('cerr'),
)
 
###############################################################################
# Source
process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

process.options = cms.untracked.PSet(
)

###############################################################################
# Global tag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '101X_dataRun2_Prompt_v11', '')

process.energyLossProducer.tag = cms.string('totem')

process.reco = cms.Path(process.MeasurementTrackerEvent
                      * process.siPixelClusterShapeCache)

###############################################################################
# Schedule
process.schedule = cms.Schedule(process.reco,
                                process.produceEnergyLoss)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)

