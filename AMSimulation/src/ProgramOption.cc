#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/ProgramOption.h"

#include <iostream>
#include <iterator>


namespace slhcl1tt {

// Print vector
template<class T>
std::ostream& operator<<(std::ostream& o, const std::vector<T>& v) {
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(o, " "));
    return o;
}

std::ostream& operator<<(std::ostream& o, const ProgramOption& po) {
    o << "Parsed program options --"
      << "  input: "        << po.input
      << "  output: "       << po.output
      << "  bankfile: "     << po.bankfile
      << "  matrixfile: "   << po.matrixfile
      << "  roadfile: "     << po.roadfile
      << "  trackfile: "    << po.trackfile

      << "  verbose: "      << po.verbose
      << "  speedup: "      << po.speedup
      << "  maxEvents: "    << po.maxEvents

      << "  nLayers: "      << po.nLayers
      << "  nFakers: "      << po.nFakers
      << "  nDCBits: "      << po.nDCBits

      << "  tower: "        << po.tower
      << "  superstrip: "   << po.superstrip
      << "  algo: "         << po.algo

      << "  minPt: "        << po.minPt
      << "  maxPt: "        << po.maxPt
      << "  minInvPt: "     << po.minInvPt
      << "  maxInvPt: "     << po.maxInvPt
      << "  minEta: "       << po.minEta
      << "  maxEta: "       << po.maxEta
      << "  minPhi: "       << po.minPhi
      << "  maxPhi: "       << po.maxPhi
      << "  minVz: "        << po.minVz
      << "  maxVz: "        << po.maxVz

      << "  picky: "        << po.picky
      << "  removeOverlap: "<< po.removeOverlap
      << "  minFrequency: " << po.minFrequency
      << "  maxPatterns: "  << po.maxPatterns
      << "  maxMisses: "    << po.maxMisses
      << "  maxStubs: "     << po.maxStubs
      << "  maxRoads: "     << po.maxRoads

      << "  view: "         << po.view
      << "  hitBits: "      << po.hitBits

      //<< "  maxChi2: "      << po.maxChi2
      << "  maxRedChi2_6out6: " << po.maxRedChi2_6out6
      << "  maxRedChi2_5out6: " << po.maxRedChi2_5out6
      << "  CutPrincipals: "<< po.CutPrincipals
      << "  minNdof: "      << po.minNdof
      << "  maxCombs: "     << po.maxCombs
      << "  maxCombsPreCB: "<< po.maxCombsPreCB
      << "  maxTracks: "    << po.maxTracks

      << "  rmDuplicate: "  << po.rmDuplicate
      << "  maxChi2Match:  "  << po.maxChi2Match
      << "  rmParDuplicate: " << po.rmParDuplicate

      << "  oldCB: "        << po.oldCB
      << "  FiveOfSix: "    << po.FiveOfSix
      << "  PDDS: "         << po.PDDS

      << "  no_trim: "      << po.no_trim

      << "  deltaS: "       << po.deltaS
      << "  deltaSM: "      << po.deltaSM

      << "  datadir: "      << po.datadir
      ;
    return o;
}

}
