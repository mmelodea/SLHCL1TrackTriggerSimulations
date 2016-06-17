#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/LocalToGlobalMap.h"

#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/Helper.h"
using namespace slhcl1tt;

#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>


// _____________________________________________________________________________
LocalToGlobalMap::LocalToGlobalMap() {
    // Intentionally left blank
}

// _____________________________________________________________________________
void LocalToGlobalMap::read(TString datadir) {
    TString csvfile = datadir + "module_localtoglobal_map.csv";

    readLocalToGlobalMap(csvfile);
}

// _____________________________________________________________________________
void LocalToGlobalMap::readLocalToGlobalMap(TString csvfile) {

    if (!csvfile.EndsWith(".csv"))
        throw std::invalid_argument("Incorrect filename.");

    // Read trigger tower map
    std::string line, line2;
    std::ifstream ifs(csvfile.Data());  // open file

    l2gmap_.clear();
    unsigned i = 0;

    while (std::getline(ifs, line)) {  // split by line break
        std::istringstream iss(line);

        if (i != 0) {  // skip the first line
            unsigned moduleId = 0;
            unsigned chipId = 0;
            std::vector<float> values;

            unsigned j = 0;
            while (std::getline(iss, line2, ',')) {  // split by comma
                if (j == 0) {
                    moduleId = std::stoi(line2);
                } else if (j == 1) {
                    chipId = std::stoi(line2);
                } else {
                    values.push_back(std::stof(line2));
                }
                ++j;
            }
            if (ifs.eof())
                break;

            assert(values.size() == 6 || iss.eof());

            LocalToGlobal l2g;
            l2g.set(values);
            l2gmap_.insert(std::make_pair(std::make_pair(moduleId, chipId), l2g));
        }
        ++i;
    }
}

// _____________________________________________________________________________
void LocalToGlobalMap::convert(const unsigned moduleId, const float strip, const float segment,
    float& conv_r, float& conv_phi, float& conv_z, LocalToGlobal& conv_l2g) {

    unsigned istrip = halfStripRound(strip);
    unsigned isegment = segmentRound(segment);
    assert(istrip < (1<<11));   // 11-bit number
    assert(isegment < (1<<5));  // 5-bit number

    const unsigned chipId = (istrip >> 8);
    istrip = istrip & 0xff;
    //const unsigned cicId = isPSModule(moduleId) ? (isegment >> 4) : isegment;
    //isegment = isegment & 0xf;

    conv_l2g = l2gmap_.at(std::make_pair(moduleId, chipId));
    if (isBarrelModule(moduleId)) {
        conv_r   = conv_l2g.x_r0   + conv_l2g.x_r   * istrip;
        conv_phi = conv_l2g.x_phi0 + conv_l2g.x_phi * istrip;
        conv_z   = conv_l2g.x_z0   + conv_l2g.x_z   * isegment;
    } else {
        conv_r   = conv_l2g.x_r0   + conv_l2g.x_r   * isegment;
        conv_phi = conv_l2g.x_phi0 + conv_l2g.x_phi * istrip;
        conv_z   = conv_l2g.x_z0   + conv_l2g.x_z   * istrip;
    }

    return;
}

// _____________________________________________________________________________
void LocalToGlobalMap::convertInt(const unsigned moduleId, const float strip, const float segment, const unsigned tt, const LocalToGlobal& conv_l2g,
    int64_t& conv_r, int64_t& conv_phi, int64_t& conv_z, LocalToGlobalInt& conv_l2g_int) {

    unsigned istrip = halfStripRound(strip);
    unsigned isegment = segmentRound(segment);
    assert(istrip < (1<<11));   // 11-bit number
    assert(isegment < (1<<5));  // 5-bit number

    float deltaPhi = 1. * 2;     // [-1, 1] rad
    float deltaZ   = 1024. * 2;  // [-1024, 1024] cm
    float deltaR   = 1024. * 2;  // [-1024, 1024] cm

    unsigned ttphi = tt%8;
    float phi0 = -M_PI/2. + (2.*M_PI/8.) * (0.5+ttphi);  // center of trigger tower

    conv_l2g_int.i_phi  = std::round((conv_l2g.x_phi  - 0.  )/deltaPhi * std::pow(2,26));  // 18 + 8 = 26
    conv_l2g_int.i_phi0 = std::round((conv_l2g.x_phi0 - phi0)/deltaPhi * std::pow(2,18));  // 18
    if (isPSModule(moduleId)) {  // is PS
        if (isBarrelModule(moduleId)) {  // is barrel
            conv_l2g_int.i_z    = std::round((conv_l2g.x_z    - 0.  )/deltaZ * std::pow(2,23));  // 18 + 5 = 23
            conv_l2g_int.i_z0   = std::round((conv_l2g.x_z0   - 0.  )/deltaZ * std::pow(2,18));  // 18
            conv_l2g_int.i_r    = std::round((conv_l2g.x_r    - 0.  )/deltaR * std::pow(2,26));
            conv_l2g_int.i_r0   = std::round((conv_l2g.x_r0   - 0.  )/deltaR * std::pow(2,18));
        } else {  // is endcap
            conv_l2g_int.i_z    = std::round((conv_l2g.x_z    - 0.  )/deltaZ * std::pow(2,26));
            conv_l2g_int.i_z0   = std::round((conv_l2g.x_z0   - 0.  )/deltaZ * std::pow(2,18));
            conv_l2g_int.i_r    = std::round((conv_l2g.x_r    - 0.  )/deltaR * std::pow(2,23));  // 18 + 5 = 23
            conv_l2g_int.i_r0   = std::round((conv_l2g.x_r0   - 0.  )/deltaR * std::pow(2,18));  // 18
        }
    } else {  // is 2S
        if (isBarrelModule(moduleId)) {  // is barrel
            conv_l2g_int.i_z    = std::round((conv_l2g.x_z    - 0.  )/deltaZ * std::pow(2,19));  // 18 + 1 = 19
            conv_l2g_int.i_z0   = std::round((conv_l2g.x_z0   - 0.  )/deltaZ * std::pow(2,18));  // 18
            conv_l2g_int.i_r    = std::round((conv_l2g.x_r    - 0.  )/deltaR * std::pow(2,26));
            conv_l2g_int.i_r0   = std::round((conv_l2g.x_r0   - 0.  )/deltaR * std::pow(2,18));
        } else {  // is endcap
            conv_l2g_int.i_z    = std::round((conv_l2g.x_z    - 0.  )/deltaZ * std::pow(2,26));
            conv_l2g_int.i_z0   = std::round((conv_l2g.x_z0   - 0.  )/deltaZ * std::pow(2,18));
            conv_l2g_int.i_r    = std::round((conv_l2g.x_r    - 0.  )/deltaR * std::pow(2,19));  // 18 + 1 = 19
            conv_l2g_int.i_r0   = std::round((conv_l2g.x_r0   - 0.  )/deltaR * std::pow(2,18));  // 18
        }
    }

    if (isPSModule(moduleId)) {  // is PS
        if (isBarrelModule(moduleId)) {  // is barrel
            conv_phi = (conv_l2g_int.i_phi * istrip  ) + (conv_l2g_int.i_phi0 << 8);
            conv_z   = (conv_l2g_int.i_z   * isegment) + (conv_l2g_int.i_z0   << 5);
            conv_r   = (conv_l2g_int.i_r   * istrip  ) + (conv_l2g_int.i_r0   << 8);
            conv_phi >>= 8;
            conv_z   >>= 5;
            conv_r   >>= 8;
        } else {  // is endcap
            conv_phi = (conv_l2g_int.i_phi * istrip  ) + (conv_l2g_int.i_phi0 << 8);
            conv_z   = (conv_l2g_int.i_z   * istrip  ) + (conv_l2g_int.i_z0   << 8);
            conv_r   = (conv_l2g_int.i_r   * isegment) + (conv_l2g_int.i_r0   << 5);
            conv_phi >>= 8;
            conv_z   >>= 8;
            conv_r   >>= 5;
        }
    } else {  // is 2S
        if (isBarrelModule(moduleId)) {  // is barrel
            conv_phi = (conv_l2g_int.i_phi * istrip  ) + (conv_l2g_int.i_phi0 << 8);
            conv_z   = (conv_l2g_int.i_z   * isegment) + (conv_l2g_int.i_z0   << 1);
            conv_r   = (conv_l2g_int.i_r   * istrip  ) + (conv_l2g_int.i_r0   << 8);
            conv_phi >>= 8;
            conv_z   >>= 1;
            conv_r   >>= 8;
        } else {  // is endcap
            conv_phi = (conv_l2g_int.i_phi * istrip  ) + (conv_l2g_int.i_phi0 << 8);
            conv_z   = (conv_l2g_int.i_z   * istrip  ) + (conv_l2g_int.i_z0   << 8);
            conv_r   = (conv_l2g_int.i_r   * isegment) + (conv_l2g_int.i_r0   << 1);
            conv_phi >>= 8;
            conv_z   >>= 8;
            conv_r   >>= 1;
        }
    }

    return;
}

// _____________________________________________________________________________
void LocalToGlobalMap::print() {

}
