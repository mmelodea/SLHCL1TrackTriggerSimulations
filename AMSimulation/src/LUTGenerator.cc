#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/LUTGenerator.h"

#include <bitset>


// _____________________________________________________________________________
int LUTGenerator::makeLocalToGlobal() {
    if (verbose_)  std::cout << Info() << "Making local-to-global LUT for superstrip conversion" << std::endl;

    float conv_r = 0., conv_phi = 0., conv_z = 0.;
    LocalToGlobal conv_l2g;
    LocalToGlobalInt conv_l2g_int;

    const std::vector<unsigned>& moduleIds = ttmap_->getTriggerTowerModules(po_.tower);
    for (unsigned i=0; i<moduleIds.size(); ++i) {

        for (unsigned j=0; j<8; ++j) {  // loop over 8 chips
            unsigned moduleId = moduleIds.at(i);
            unsigned chipId   = j;
            float    strip    = 0.5 * (chipId * 256);  // in full-strip unit
            float    segment  = 0;

            // Do local-to-global conversion
            l2gmap_ -> convert(moduleId, strip, segment, conv_r, conv_phi, conv_z, conv_l2g);

            // Find superstrip ID
            unsigned ssId = arbiter_ -> superstripLocal(moduleId, strip, segment, conv_l2g, conv_l2g_int);

            std::cout << moduleId << ", " << chipId << ", "
                      << std::bitset<64>(conv_l2g_int.i_phi0) << ", " << std::bitset<64>(conv_l2g_int.i_phi) << ", "
                      << std::bitset<64>(conv_l2g_int.i_z0  ) << ", " << std::bitset<64>(conv_l2g_int.i_z  ) << ", "
                      << std::bitset<64>(conv_l2g_int.i_r0  ) << ", " << std::bitset<64>(conv_l2g_int.i_r  )
                      << std::endl;
        }
    }

    return 0;
}

// _____________________________________________________________________________
int LUTGenerator::makeLocalToGlobalStar() {
    if (verbose_)  std::cout << Info() << "Making local-to-global LUT for global r,phi,z conversion" << std::endl;

    float conv_r = 0., conv_phi = 0., conv_z = 0.;
    int64_t conv_r_int = 0., conv_phi_int = 0., conv_z_int = 0.;
    LocalToGlobal conv_l2g;
    LocalToGlobalInt conv_l2g_int;

    const std::vector<unsigned>& moduleIds = ttmap_->getTriggerTowerModules(po_.tower);
    for (unsigned i=0; i<moduleIds.size(); ++i) {

        for (unsigned j=0; j<8; ++j) {  // loop over 8 chips
            unsigned moduleId = moduleIds.at(i);
            unsigned chipId   = j;
            float    strip    = 0.5 * (chipId * 256);  // in full-strip unit
            float    segment  = 0;

            // Do local-to-global conversion
            l2gmap_ -> convert(moduleId, strip, segment, conv_r, conv_phi, conv_z, conv_l2g);

            l2gmap_ -> convertInt(moduleId, strip, segment, po_.tower, conv_l2g, conv_r_int, conv_phi_int, conv_z_int, conv_l2g_int);

            std::cout << moduleId << ", " << chipId << ", "
                      << std::bitset<64>(conv_l2g_int.i_phi0) << ", " << std::bitset<64>(conv_l2g_int.i_phi) << ", "
                      << std::bitset<64>(conv_l2g_int.i_z0  ) << ", " << std::bitset<64>(conv_l2g_int.i_z  ) << ", "
                      << std::bitset<64>(conv_l2g_int.i_r0  ) << ", " << std::bitset<64>(conv_l2g_int.i_r  )
                      << std::endl;
        }
    }

    return 0;
}


// _____________________________________________________________________________
// Main driver
int LUTGenerator::run() {
    int exitcode = 0;
    Timing(1);

    exitcode = makeLocalToGlobal();
    if (exitcode)  return exitcode;
    Timing();

    exitcode = makeLocalToGlobalStar();
    if (exitcode)  return exitcode;
    Timing();

    return exitcode;
}
