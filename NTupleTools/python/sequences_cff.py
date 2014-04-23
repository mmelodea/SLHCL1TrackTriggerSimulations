import FWCore.ParameterSet.Config as cms

from SLHCL1TrackTriggerSimulations.NTupleTools.ntupleGen_cfi import *
from SLHCL1TrackTriggerSimulations.NTupleTools.ntupleSim_cfi import *
from SLHCL1TrackTriggerSimulations.NTupleTools.ntupleL1Track_cfi import *
from SLHCL1TrackTriggerSimulations.NTupleTools.ntupleMaker_cfi import *

ntupleSequence = cms.Sequence(ntupleGen * ntupleSim * ntupleL1Track * (ntupleEventInfo * ntupler))
