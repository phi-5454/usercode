from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'addDeDx_reco'

config.section_('JobType')
config.JobType.psetName    = 'addDeDx_reco.py'
config.JobType.pluginName  = 'Analysis'
# config.JobType.outputFiles = ['result.dat']

config.section_('Data')
config.Data.inputDataset = '/TOTEM40/Run2018B-22Feb2019-v1/RECO'
config.Data.lumiMask = '90m.json'
config.Data.splitting   = 'LumiBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits  = 3
config.Data.publication = False

config.section_('User')

config.section_('Site')
config.Site.storageSite = 'T2_HU_Budapest'
