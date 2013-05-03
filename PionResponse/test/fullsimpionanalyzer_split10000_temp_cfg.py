import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load('Configuration.StandardSequences.Services_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag = autoCond['startup']

#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'MC_53_V6::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# this inputs the input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
		'file:pion_ENERGYIN_500k_part1.root',
		'file:pion_ENERGYIN_500k_part2.root',
		'file:pion_ENERGYIN_500k_part3.root',
		'file:pion_ENERGYIN_500k_part4.root',
		'file:pion_ENERGYIN_500k_part5.root',
		'file:pion_ENERGYIN_500k_part6.root',
		'file:pion_ENERGYIN_500k_part7.root',
		'file:pion_ENERGYIN_500k_part8.root',
		'file:pion_ENERGYIN_500k_part9.root',
		'file:pion_ENERGYIN_500k_part10.root',
		'file:pion_ENERGYIN_500k_part11.root',
		'file:pion_ENERGYIN_500k_part12.root',
		'file:pion_ENERGYIN_500k_part13.root',
		'file:pion_ENERGYIN_500k_part14.root',
		'file:pion_ENERGYIN_500k_part15.root',
		'file:pion_ENERGYIN_500k_part16.root',
		'file:pion_ENERGYIN_500k_part17.root',
		'file:pion_ENERGYIN_500k_part18.root',
		'file:pion_ENERGYIN_500k_part19.root',
		'file:pion_ENERGYIN_500k_part20.root',
		'file:pion_ENERGYIN_500k_part21.root',
		'file:pion_ENERGYIN_500k_part22.root',
		'file:pion_ENERGYIN_500k_part23.root',
		'file:pion_ENERGYIN_500k_part24.root',
		'file:pion_ENERGYIN_500k_part25.root',
		'file:pion_ENERGYIN_500k_part26.root',
		'file:pion_ENERGYIN_500k_part27.root',
		'file:pion_ENERGYIN_500k_part28.root',
		'file:pion_ENERGYIN_500k_part29.root',
		'file:pion_ENERGYIN_500k_part30.root',
		'file:pion_ENERGYIN_500k_part31.root',
		'file:pion_ENERGYIN_500k_part32.root',
		'file:pion_ENERGYIN_500k_part33.root',
		'file:pion_ENERGYIN_500k_part34.root',
		'file:pion_ENERGYIN_500k_part35.root',
		'file:pion_ENERGYIN_500k_part36.root',
		'file:pion_ENERGYIN_500k_part37.root',
		'file:pion_ENERGYIN_500k_part38.root',
		'file:pion_ENERGYIN_500k_part39.root',
		'file:pion_ENERGYIN_500k_part40.root',
		'file:pion_ENERGYIN_500k_part41.root',
		'file:pion_ENERGYIN_500k_part42.root',
		'file:pion_ENERGYIN_500k_part43.root',
		'file:pion_ENERGYIN_500k_part44.root',
		'file:pion_ENERGYIN_500k_part45.root',
		'file:pion_ENERGYIN_500k_part46.root',
		'file:pion_ENERGYIN_500k_part47.root',
		'file:pion_ENERGYIN_500k_part48.root',
		'file:pion_ENERGYIN_500k_part49.root',
		'file:pion_ENERGYIN_500k_part50.root',
   ),
	duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)

process.demo = cms.EDAnalyzer('FullSimPionAnalyzer',
    fileName = cms.string("tree_pion_ENERGYIN_500k.root"),
	dRcut = cms.double(0.5),
	samplingHBHE = cms.vdouble(125.44, 125.54, 125.32, 125.13, 124.46,
                               125.01, 125.22, 125.48, 124.45, 125.90,
                               125.83, 127.01, 126.82, 129.73, 131.83,
                               143.52, # HB
                               210.55, 197.93, 186.12, 189.64, 189.63,
                               190.28, 189.61, 189.60, 190.12, 191.22,
                               190.90, 193.06, 188.42, 188.42), #HE
    samplingHF   = cms.vdouble(0.383, 0.368),
	simHitToPhotoelectrons = cms.vdouble(6.0,6.0),
    samplingHO   = cms.vdouble(231.0, 231.0, 231.0, 231.0, 360.0, 
                               360.0, 360.0, 360.0, 360.0, 360.0,
                               360.0, 360.0, 360.0, 360.0, 360.0)
)

#random number seed
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    demo = cms.PSet(
        initialSeed = cms.untracked.uint32(83),
        engineName = cms.untracked.string('TRandom3')
    )
)

process.p = cms.Path(process.demo)

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#SOURCE process.load('#SRC')
