import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process("reRECO", eras.Run2_2018_highBetaStar)

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("UserCode.ProtonReco.ProtonProducer_cff")

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
      reportEvery = cms.untracked.int32(  10000),
      limit       = cms.untracked.int32(1000000)
    )
  ),
  destinations = cms.untracked.vstring('cerr'),
)
 
###############################################################################
# Source
process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    skipBadFiles = cms.untracked.bool(True),
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
 'root://eostotem//eos/totem/data/cmstotem/2018/90m/RECO_copy/TOTEM23/120000/34B6109F-8450-E911-B40A-98039B3B01D2.root', 'root://eostotem//eos/totem/data/cmstotem/2018/90m/RECO_copy/TOTEM23/120000/34F9C73F-CB50-E911-AB8B-3417EBE706C3.root', 'root://eostotem//eos/totem/data/cmstotem/2018/90m/RECO_copy/TOTEM23/120000/3600B3CD-CB50-E911-8FAA-002590DE6E52.root', 'root://eostotem//eos/totem/data/cmstotem/2018/90m/RECO_copy/TOTEM23/120000/362959DA-0D51-E911-9736-3417EBE700A2.root', 'root://eostotem//eos/totem/data/cmstotem/2018/90m/RECO_copy/TOTEM23/120000/364AEE76-9350-E911-BC50-3417EBE74303.root', 'root://eostotem//eos/totem/data/cmstotem/2018/90m/RECO_copy/TOTEM23/120000/36579D40-AE50-E911-8392-98039B3B004E.root', 'root://eostotem//eos/totem/data/cmstotem/2018/90m/RECO_copy/TOTEM23/120000/366FD4D0-C450-E911-B8D4-A4BF012881D0.root', 'root://eostotem//eos/totem/data/cmstotem/2018/90m/RECO_copy/TOTEM23/120000/369083ED-9650-E911-9EF7-002590DE6DE4.root', 'root://eostotem//eos/totem/data/cmstotem/2018/90m/RECO_copy/TOTEM23/120000/36E3A2A2-C950-E911-9763-3417EBE7063F.root', 'root://eostotem//eos/totem/data/cmstotem/2018/90m/RECO_copy/TOTEM23/120000/3802F4C5-C050-E911-B76B-3417EBE705CD.root',
    )
)

process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(100)
# input = cms.untracked.int32(-1)
)

process.options = cms.untracked.PSet(
)

###############################################################################
# Global tag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '101X_dataRun2_Prompt_v11', '')

###############################################################################
# Re-reco
process.load("RecoCTPPS.TotemRPLocal.totemRPLocalReconstruction_cff")

process.totemRPUVPatternFinder.tagRecHit = cms.InputTag(
 'totemRPRecHitProducer', '','reRECO')
process.totemRPUVPatternFinder.max_a_toFit = cms.double(0.05)

process.totemRPLocalTrackFitter.tagUVPattern = cms.InputTag(
 'totemRPUVPatternFinder','','reRECO')

process.protonProducer.outFile = cms.string(
 'out.dat.gz'
)

#
process.rereco = cms.Path(process.totemRPRecHitProducer
                        * process.totemRPUVPatternFinder
                        * process.totemRPLocalTrackFitter)

###############################################################################
# Schedule
process.schedule = cms.Schedule(process.rereco,
                                process.produceProtons)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)

