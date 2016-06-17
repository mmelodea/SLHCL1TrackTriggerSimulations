#ifndef AMSimulation_LUTGenerator_h_
#define AMSimulation_LUTGenerator_h_

#include "SLHCL1TrackTriggerSimulations/AMSimulationDataFormats/interface/Pattern.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/Helper.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/ProgramOption.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/TriggerTowerMap.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/LocalToGlobalMap.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/SuperstripArbiter.h"
using namespace slhcl1tt;


class LUTGenerator {
  public:
    // Constructor
    LUTGenerator(const ProgramOption& po)
    : po_(po),
      nEvents_(po.maxEvents), verbose_(po.verbose),
      prefix_(""), suffix_("") {

        // Initialize
        ttmap_ = new TriggerTowerMap();
        ttmap_->read(po_.datadir);

        l2gmap_ = new LocalToGlobalMap();
        l2gmap_->read(po_.datadir);

        arbiter_ = new SuperstripArbiter();
        arbiter_->setDefinition(po_.superstrip, po_.tower, ttmap_);
    }

    // Destructor
    ~LUTGenerator() {
        if (ttmap_)     delete ttmap_;
        if (l2gmap_)    delete l2gmap_;
        if (arbiter_)   delete arbiter_;
    }

    // Main driver
    int run();


  private:
    // Member functions

    // Make local-to-global LUT for superstrip conversion
    int makeLocalToGlobal();

    // Make local-to-global LUT for global r,phi,z conversion
    int makeLocalToGlobalStar();

    // Program options
    const ProgramOption po_;
    long long nEvents_;
    int verbose_;

    // Configurations
    const TString prefix_;
    const TString suffix_;

    // Operators
    TriggerTowerMap   * ttmap_;
    LocalToGlobalMap  * l2gmap_;
    SuperstripArbiter * arbiter_;
};

#endif
