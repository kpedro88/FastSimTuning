# Auto generated configuration file
# using: 
# Revision: 1.386 
# Source: /local/reps/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: Configuration/Generator/python/SinglePiE50HCAL_cfi.py -s GEN,SIM --conditions auto:mc --datatier GEN-SIM-RAW --eventcontent RAWSIM -n 10 --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('SIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load("Configuration.StandardSequences.GeometryDB_cff")
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load('Configuration.StandardSequences.GeometrySimDB_cff')
#process.load('Configuration.Geometry.GeometryECALHCAL_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
#process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
#process.load('IOMC.EventVertexGenerators.VtxSmearedGaus_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(SPLITNUMBER)
)

#configuration changes
process.g4SimHits.UseMagneticField = cms.bool(False)

#turn off vertex smearing
process.VtxSmeared = cms.EDProducer("GaussEvtVtxGenerator",
    src = cms.InputTag("generator"),
    readDB = cms.bool(False),
    MeanX = cms.double(0.),
    MeanY = cms.double(0.),
    MeanZ = cms.double(0.),
    SigmaX = cms.double(0),
    SigmaY = cms.double(0),
    SigmaZ = cms.double(0),
	TimeOffset = cms.double(0.0)
   )

#random number seed
process.RandomNumberGeneratorService.generator.initialSeed = cms.untracked.uint32(ENERGYINPNUMBER)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.386 $'),
    annotation = cms.untracked.string('Configuration/Generator/python/SinglePiE50HCAL_cfi.py nevts:10'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    #outputCommands = process.RAWSIMEventContent.outputCommands,
	outputCommands = cms.untracked.vstring('drop *', 
                                           'keep PCaloHits_g4SimHits_EcalHitsEB_*',
										   'keep PCaloHits_g4SimHits_EcalHitsEE_*',
										   'keep PCaloHits_g4SimHits_EcalHitsES_*',
										   'keep PCaloHits_g4SimHits_HcalHits_*',
										   'keep recoGenParticles_genParticles_*_*',
										   ),
    fileName = cms.untracked.string('/data/users/pedrok/FullSim/response/pion_ENERGYIN_500k_partPNUMBER.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RAW')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag = autoCond['startup']

process.generator = cms.EDProducer("FlatRandomEGunProducer",
    PGunParameters = cms.PSet(
        MinE = cms.double(ENERGYIN),
        MaxE = cms.double(ENERGYIN),
		# Pi+ = 211, Pi0 = 111, nu_e = 12, e = 11, K_long = 130
        PartID = cms.vint32(-211),
		# Eta limits: ECAL end-cap = [1.6,2.5], HCAL end-cap = [1.39, 3.0], HCAL forward = [3.0, 5.0]
        MinEta = cms.double(-5.0),
        MaxEta = cms.double(5.0),
        MinPhi = cms.double(-3.14159265359),
        MaxPhi = cms.double(3.14159265359)
    ),
    Verbosity = cms.untracked.int32(0),
    psethack = cms.string('single pi E 50 HCAL'),
    AddAntiParticle = cms.bool(False),
    firstRun = cms.untracked.uint32(1)
)


# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.RAWSIMoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 

