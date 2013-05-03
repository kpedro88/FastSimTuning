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
		'file:pion_ENERGYIN_500k_part51.root',
		'file:pion_ENERGYIN_500k_part52.root',
		'file:pion_ENERGYIN_500k_part53.root',
		'file:pion_ENERGYIN_500k_part54.root',
		'file:pion_ENERGYIN_500k_part55.root',
		'file:pion_ENERGYIN_500k_part56.root',
		'file:pion_ENERGYIN_500k_part57.root',
		'file:pion_ENERGYIN_500k_part58.root',
		'file:pion_ENERGYIN_500k_part59.root',
		'file:pion_ENERGYIN_500k_part60.root',
		'file:pion_ENERGYIN_500k_part61.root',
		'file:pion_ENERGYIN_500k_part62.root',
		'file:pion_ENERGYIN_500k_part63.root',
		'file:pion_ENERGYIN_500k_part64.root',
		'file:pion_ENERGYIN_500k_part65.root',
		'file:pion_ENERGYIN_500k_part66.root',
		'file:pion_ENERGYIN_500k_part67.root',
		'file:pion_ENERGYIN_500k_part68.root',
		'file:pion_ENERGYIN_500k_part69.root',
		'file:pion_ENERGYIN_500k_part70.root',
		'file:pion_ENERGYIN_500k_part71.root',
		'file:pion_ENERGYIN_500k_part72.root',
		'file:pion_ENERGYIN_500k_part73.root',
		'file:pion_ENERGYIN_500k_part74.root',
		'file:pion_ENERGYIN_500k_part75.root',
		'file:pion_ENERGYIN_500k_part76.root',
		'file:pion_ENERGYIN_500k_part77.root',
		'file:pion_ENERGYIN_500k_part78.root',
		'file:pion_ENERGYIN_500k_part79.root',
		'file:pion_ENERGYIN_500k_part80.root',
		'file:pion_ENERGYIN_500k_part81.root',
		'file:pion_ENERGYIN_500k_part82.root',
		'file:pion_ENERGYIN_500k_part83.root',
		'file:pion_ENERGYIN_500k_part84.root',
		'file:pion_ENERGYIN_500k_part85.root',
		'file:pion_ENERGYIN_500k_part86.root',
		'file:pion_ENERGYIN_500k_part87.root',
		'file:pion_ENERGYIN_500k_part88.root',
		'file:pion_ENERGYIN_500k_part89.root',
		'file:pion_ENERGYIN_500k_part90.root',
		'file:pion_ENERGYIN_500k_part91.root',
		'file:pion_ENERGYIN_500k_part92.root',
		'file:pion_ENERGYIN_500k_part93.root',
		'file:pion_ENERGYIN_500k_part94.root',
		'file:pion_ENERGYIN_500k_part95.root',
		'file:pion_ENERGYIN_500k_part96.root',
		'file:pion_ENERGYIN_500k_part97.root',
		'file:pion_ENERGYIN_500k_part98.root',
		'file:pion_ENERGYIN_500k_part99.root',
		'file:pion_ENERGYIN_500k_part100.root',
		'file:pion_ENERGYIN_500k_part101.root',
		'file:pion_ENERGYIN_500k_part102.root',
		'file:pion_ENERGYIN_500k_part103.root',
		'file:pion_ENERGYIN_500k_part104.root',
		'file:pion_ENERGYIN_500k_part105.root',
		'file:pion_ENERGYIN_500k_part106.root',
		'file:pion_ENERGYIN_500k_part107.root',
		'file:pion_ENERGYIN_500k_part108.root',
		'file:pion_ENERGYIN_500k_part109.root',
		'file:pion_ENERGYIN_500k_part110.root',
		'file:pion_ENERGYIN_500k_part111.root',
		'file:pion_ENERGYIN_500k_part112.root',
		'file:pion_ENERGYIN_500k_part113.root',
		'file:pion_ENERGYIN_500k_part114.root',
		'file:pion_ENERGYIN_500k_part115.root',
		'file:pion_ENERGYIN_500k_part116.root',
		'file:pion_ENERGYIN_500k_part117.root',
		'file:pion_ENERGYIN_500k_part118.root',
		'file:pion_ENERGYIN_500k_part119.root',
		'file:pion_ENERGYIN_500k_part120.root',
		'file:pion_ENERGYIN_500k_part121.root',
		'file:pion_ENERGYIN_500k_part122.root',
		'file:pion_ENERGYIN_500k_part123.root',
		'file:pion_ENERGYIN_500k_part124.root',
		'file:pion_ENERGYIN_500k_part125.root',
		'file:pion_ENERGYIN_500k_part126.root',
		'file:pion_ENERGYIN_500k_part127.root',
		'file:pion_ENERGYIN_500k_part128.root',
		'file:pion_ENERGYIN_500k_part129.root',
		'file:pion_ENERGYIN_500k_part130.root',
		'file:pion_ENERGYIN_500k_part131.root',
		'file:pion_ENERGYIN_500k_part132.root',
		'file:pion_ENERGYIN_500k_part133.root',
		'file:pion_ENERGYIN_500k_part134.root',
		'file:pion_ENERGYIN_500k_part135.root',
		'file:pion_ENERGYIN_500k_part136.root',
		'file:pion_ENERGYIN_500k_part137.root',
		'file:pion_ENERGYIN_500k_part138.root',
		'file:pion_ENERGYIN_500k_part139.root',
		'file:pion_ENERGYIN_500k_part140.root',
		'file:pion_ENERGYIN_500k_part141.root',
		'file:pion_ENERGYIN_500k_part142.root',
		'file:pion_ENERGYIN_500k_part143.root',
		'file:pion_ENERGYIN_500k_part144.root',
		'file:pion_ENERGYIN_500k_part145.root',
		'file:pion_ENERGYIN_500k_part146.root',
		'file:pion_ENERGYIN_500k_part147.root',
		'file:pion_ENERGYIN_500k_part148.root',
		'file:pion_ENERGYIN_500k_part149.root',
		'file:pion_ENERGYIN_500k_part150.root',
		'file:pion_ENERGYIN_500k_part151.root',
		'file:pion_ENERGYIN_500k_part152.root',
		'file:pion_ENERGYIN_500k_part153.root',
		'file:pion_ENERGYIN_500k_part154.root',
		'file:pion_ENERGYIN_500k_part155.root',
		'file:pion_ENERGYIN_500k_part156.root',
		'file:pion_ENERGYIN_500k_part157.root',
		'file:pion_ENERGYIN_500k_part158.root',
		'file:pion_ENERGYIN_500k_part159.root',
		'file:pion_ENERGYIN_500k_part160.root',
		'file:pion_ENERGYIN_500k_part161.root',
		'file:pion_ENERGYIN_500k_part162.root',
		'file:pion_ENERGYIN_500k_part163.root',
		'file:pion_ENERGYIN_500k_part164.root',
		'file:pion_ENERGYIN_500k_part165.root',
		'file:pion_ENERGYIN_500k_part166.root',
		'file:pion_ENERGYIN_500k_part167.root',
		'file:pion_ENERGYIN_500k_part168.root',
		'file:pion_ENERGYIN_500k_part169.root',
		'file:pion_ENERGYIN_500k_part170.root',
		'file:pion_ENERGYIN_500k_part171.root',
		'file:pion_ENERGYIN_500k_part172.root',
		'file:pion_ENERGYIN_500k_part173.root',
		'file:pion_ENERGYIN_500k_part174.root',
		'file:pion_ENERGYIN_500k_part175.root',
		'file:pion_ENERGYIN_500k_part176.root',
		'file:pion_ENERGYIN_500k_part177.root',
		'file:pion_ENERGYIN_500k_part178.root',
		'file:pion_ENERGYIN_500k_part179.root',
		'file:pion_ENERGYIN_500k_part180.root',
		'file:pion_ENERGYIN_500k_part181.root',
		'file:pion_ENERGYIN_500k_part182.root',
		'file:pion_ENERGYIN_500k_part183.root',
		'file:pion_ENERGYIN_500k_part184.root',
		'file:pion_ENERGYIN_500k_part185.root',
		'file:pion_ENERGYIN_500k_part186.root',
		'file:pion_ENERGYIN_500k_part187.root',
		'file:pion_ENERGYIN_500k_part188.root',
		'file:pion_ENERGYIN_500k_part189.root',
		'file:pion_ENERGYIN_500k_part190.root',
		'file:pion_ENERGYIN_500k_part191.root',
		'file:pion_ENERGYIN_500k_part192.root',
		'file:pion_ENERGYIN_500k_part193.root',
		'file:pion_ENERGYIN_500k_part194.root',
		'file:pion_ENERGYIN_500k_part195.root',
		'file:pion_ENERGYIN_500k_part196.root',
		'file:pion_ENERGYIN_500k_part197.root',
		'file:pion_ENERGYIN_500k_part198.root',
		'file:pion_ENERGYIN_500k_part199.root',
		'file:pion_ENERGYIN_500k_part200.root',
		'file:pion_ENERGYIN_500k_part201.root',
		'file:pion_ENERGYIN_500k_part202.root',
		'file:pion_ENERGYIN_500k_part203.root',
		'file:pion_ENERGYIN_500k_part204.root',
		'file:pion_ENERGYIN_500k_part205.root',
		'file:pion_ENERGYIN_500k_part206.root',
		'file:pion_ENERGYIN_500k_part207.root',
		'file:pion_ENERGYIN_500k_part208.root',
		'file:pion_ENERGYIN_500k_part209.root',
		'file:pion_ENERGYIN_500k_part210.root',
		'file:pion_ENERGYIN_500k_part211.root',
		'file:pion_ENERGYIN_500k_part212.root',
		'file:pion_ENERGYIN_500k_part213.root',
		'file:pion_ENERGYIN_500k_part214.root',
		'file:pion_ENERGYIN_500k_part215.root',
		'file:pion_ENERGYIN_500k_part216.root',
		'file:pion_ENERGYIN_500k_part217.root',
		'file:pion_ENERGYIN_500k_part218.root',
		'file:pion_ENERGYIN_500k_part219.root',
		'file:pion_ENERGYIN_500k_part220.root',
		'file:pion_ENERGYIN_500k_part221.root',
		'file:pion_ENERGYIN_500k_part222.root',
		'file:pion_ENERGYIN_500k_part223.root',
		'file:pion_ENERGYIN_500k_part224.root',
		'file:pion_ENERGYIN_500k_part225.root',
		'file:pion_ENERGYIN_500k_part226.root',
		'file:pion_ENERGYIN_500k_part227.root',
		'file:pion_ENERGYIN_500k_part228.root',
		'file:pion_ENERGYIN_500k_part229.root',
		'file:pion_ENERGYIN_500k_part230.root',
		'file:pion_ENERGYIN_500k_part231.root',
		'file:pion_ENERGYIN_500k_part232.root',
		'file:pion_ENERGYIN_500k_part233.root',
		'file:pion_ENERGYIN_500k_part234.root',
		'file:pion_ENERGYIN_500k_part235.root',
		'file:pion_ENERGYIN_500k_part236.root',
		'file:pion_ENERGYIN_500k_part237.root',
		'file:pion_ENERGYIN_500k_part238.root',
		'file:pion_ENERGYIN_500k_part239.root',
		'file:pion_ENERGYIN_500k_part240.root',
		'file:pion_ENERGYIN_500k_part241.root',
		'file:pion_ENERGYIN_500k_part242.root',
		'file:pion_ENERGYIN_500k_part243.root',
		'file:pion_ENERGYIN_500k_part244.root',
		'file:pion_ENERGYIN_500k_part245.root',
		'file:pion_ENERGYIN_500k_part246.root',
		'file:pion_ENERGYIN_500k_part247.root',
		'file:pion_ENERGYIN_500k_part248.root',
		'file:pion_ENERGYIN_500k_part249.root',
		'file:pion_ENERGYIN_500k_part250.root',
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