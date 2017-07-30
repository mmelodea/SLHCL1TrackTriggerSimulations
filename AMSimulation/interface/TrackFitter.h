#ifndef AMSimulation_TrackFitter_h_
#define AMSimulation_TrackFitter_h_

#include "SLHCL1TrackTriggerSimulations/AMSimulationDataFormats/interface/TTTrack2.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTTrackReader.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTTrackWriter.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/Helper.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/ProgramOption.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/TrackFitterAlgoPCA.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/TrackFitterAlgoATF.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/TrackFitterAlgoLTF.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/CombinationFactory.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/CombinationBuilderFactory.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/PairCombinationFactory.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/DuplicateRemoval.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/ParameterDuplicateRemoval.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/MCTruthAssociator.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/HTMemory.h"
using namespace slhcl1tt;


class TrackFitter {
  public:
    typedef TTTrackReaderT<kTrackFitter> TTRoadReader;
    typedef TTTrackWriterT<kTrackFitter> TTTrackWriter;

    // Constructor
    TrackFitter(const ProgramOption& po) :
      po_(po),
      nEvents_(po.maxEvents), verbose_(po.verbose) {

        // Decide the track fitter to use
        fitter_ = 0;
        if (po.algo == "PCA4" || po.algo == "PCA5") {
            fitter_ = new TrackFitterAlgoPCA(po);
        } else if (po.algo == "ATF4") {
            fitter_ = new TrackFitterAlgoATF(false);
        } else if (po.algo == "ATF5") {
            fitter_ = new TrackFitterAlgoATF(true);
        } else if (po.algo == "LTF") {
            fitter_ = new TrackFitterAlgoLTF(po);
        } else {
            throw std::invalid_argument("unknown track fitter algo.");
        }

        combinationBuilderFactory_ = std::make_shared<CombinationBuilderFactory>(po_.FiveOfSix);

        try {
            htm_.reset(new HTMemory("htm", po.htmfile));
        } catch (...) {
            throw;
        }
    }

    // Destructor
    ~TrackFitter() {
        if (fitter_)  delete fitter_;
    }

    // Main driver
    int run();


  private:
    // Member functions
    int makeTracks(TString src, TString out);

    // Program options
    const ProgramOption po_;
    long long nEvents_;
    int verbose_;

    // Track fitter
    TrackFitterAlgoBase * fitter_;

    // Combination factory
    CombinationFactory combinationFactory_;

    // pair combination factory
    PairCombinationFactory pairCombinationFactory_;

    // SCB and ACB combination factory
    std::shared_ptr<CombinationBuilderFactory> combinationBuilderFactory_;

    // MC truth associator
    MCTruthAssociator truthAssociator_;

    // Hough Transform filter
    std::unique_ptr<HTMemory> htm_;
};

#endif
