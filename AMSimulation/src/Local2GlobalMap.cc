#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/Local2GlobalMap.h"

#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/Helper.h"
using namespace slhcl1tt;

#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>


// _____________________________________________________________________________
Local2GlobalMap::Local2GlobalMap() {
    // Intentionally left blank
}

// _____________________________________________________________________________
void Local2GlobalMap::read(TString datadir) {
    TString csvfile = datadir + "module_local2global_map.csv";

    readLocal2GlobalMap(csvfile);
}

// _____________________________________________________________________________
void Local2GlobalMap::readLocal2GlobalMap(TString csvfile) {

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

            Local2Global l2g;
            l2g.set(values);
            l2gmap_.insert(std::make_pair(std::make_pair(moduleId, chipId), l2g));
        }
        ++i;
    }
}

// _____________________________________________________________________________
void Local2GlobalMap::convert(const unsigned moduleId, const float strip, const float segment,
    float& conv_r, float& conv_phi, float& conv_z, Local2Global& conv_l2g) {

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
void Local2GlobalMap::print() {

}
