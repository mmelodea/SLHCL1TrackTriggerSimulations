#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/LUTGenerator.h"

#include <bitset>
#include <fstream>
#include <iostream>


// _____________________________________________________________________________
unsigned LUTGenerator::localModuleId(unsigned moduleId) {
    unsigned lay16 = compressLayer(decodeLayer(moduleId));

    if (po_.tower == 27) {
        unsigned offsets[16] = {50230, 60326, 70426, 80611, 90711, 100811, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        moduleId -= offsets[lay16];
    }

    unsigned lad = decodeLadder(moduleId);
    unsigned mod = decodeModule(moduleId);
    unsigned newModuleId = ((lad & 0xf) << 4) | (mod & 0xf);

    return newModuleId;
}

// _____________________________________________________________________________
int LUTGenerator::makeLocalToGlobal() {
    if (verbose_)  std::cout << Info() << "Making local-to-global LUT for superstrip conversion" << std::endl;

    float conv_r = 0., conv_phi = 0., conv_z = 0.;
    LocalToGlobal conv_l2g;
    LocalToGlobalInt conv_l2g_int;

    std::vector<std::string> bitStrings;
    std::vector<unsigned> bitStrings_breaks(6, 0);

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

            if (verbose_ > 1) {
                std::cout << moduleId << ", " << chipId << ", "
                          << std::bitset<64>(conv_l2g_int.i_phi0) << ", " << std::bitset<64>(conv_l2g_int.i_phi) << ", "
                          << std::bitset<64>(conv_l2g_int.i_z0  ) << ", " << std::bitset<64>(conv_l2g_int.i_z  ) << ", "
                          << std::bitset<64>(conv_l2g_int.i_r0  ) << ", " << std::bitset<64>(conv_l2g_int.i_r  )
                          << std::endl;
            }

            unsigned newModuleId = localModuleId(moduleId);
            std::string bitString = "";
            bitString += std::bitset<18>(conv_l2g_int.i_phi0).to_string();
            bitString += std::bitset<18>(conv_l2g_int.i_phi).to_string();
            bitString += std::bitset<18>(conv_l2g_int.i_z0).to_string();
            bitString += std::bitset<18>(conv_l2g_int.i_z).to_string();
            bitString += std::bitset<8>(newModuleId).to_string();
            bitString += std::bitset<3>(chipId).to_string();
            bitStrings.push_back(bitString);

            unsigned lay16 = compressLayer(decodeLayer(moduleId));
            bitStrings_breaks.at(lay16) = bitStrings.size();
        }
    }

    for (unsigned i=0; i<bitStrings_breaks.size(); ++i) {
        // Open text file
        TString txtname = Form("lut_L2G_l%i.txt", i+5);
        ofstream txtfile(txtname.Data());

        // Write text file
        std::vector<std::string>::const_iterator it=bitStrings.begin()+(i == 0 ? 0 : bitStrings_breaks.at(i-1));
        for (; it!=bitStrings.begin() + bitStrings_breaks.at(i); ++it) {
            txtfile << *it << std::endl;
            if (verbose_ > 2) {
                std::cout << *it << std::endl;
            }
        }

        // Close text file
        txtfile.close();
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

    std::vector<std::string> bitStrings;
    std::vector<unsigned> bitStrings_breaks(6, 0);

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

            if (verbose_ > 1) {
                std::cout << moduleId << ", " << chipId << ", "
                          << std::bitset<64>(conv_l2g_int.i_phi0) << ", " << std::bitset<64>(conv_l2g_int.i_phi) << ", "
                          << std::bitset<64>(conv_l2g_int.i_z0  ) << ", " << std::bitset<64>(conv_l2g_int.i_z  ) << ", "
                          << std::bitset<64>(conv_l2g_int.i_r0  ) << ", " << std::bitset<64>(conv_l2g_int.i_r  )
                          << std::endl;
            }

            unsigned newModuleId = localModuleId(moduleId);
            std::string bitString = "";
            bitString += std::bitset<18>(conv_l2g_int.i_phi0).to_string();
            bitString += std::bitset<18>(conv_l2g_int.i_phi).to_string();
            bitString += std::bitset<18>(conv_l2g_int.i_z0).to_string();
            bitString += std::bitset<18>(conv_l2g_int.i_z).to_string();
            bitString += std::bitset<18>(conv_l2g_int.i_r0).to_string();
            bitString += std::bitset<18>(conv_l2g_int.i_r).to_string();
            bitString += std::bitset<8>(newModuleId).to_string();
            bitString += std::bitset<3>(chipId).to_string();
            bitStrings.push_back(bitString);

            unsigned lay16 = compressLayer(decodeLayer(moduleId));
            bitStrings_breaks.at(lay16) = bitStrings.size();
        }
    }

    for (unsigned i=0; i<bitStrings_breaks.size(); ++i) {
        // Open text file
        TString txtname = Form("lut_L2GStar_l%i.txt", i+5);
        ofstream txtfile(txtname.Data());

        // Write text file
        std::vector<std::string>::const_iterator it=bitStrings.begin()+(i == 0 ? 0 : bitStrings_breaks.at(i-1));
        for (; it!=bitStrings.begin()+bitStrings_breaks.at(i); ++it) {
            txtfile << *it << std::endl;
            if (verbose_ > 2) {
                std::cout << *it << std::endl;
            }
        }

        // Close text file
        txtfile.close();
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
