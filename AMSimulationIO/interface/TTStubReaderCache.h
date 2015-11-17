#ifndef AMSimulationIO_TTStubReaderCache_h_
#define AMSimulationIO_TTStubReaderCache_h_

#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTStubReader.h"

namespace slhcl1tt {

class TTStubReaderCache {
  public:
    TTStubReaderCache(const TTStubReader& reader, int verbose=1);
    ~TTStubReaderCache();

    // Keep these branches synchronized with TTStubReader!

    // genParticle information
    std::vector<float> *          vp_pt;
    std::vector<float> *          vp_eta;
    std::vector<float> *          vp_phi;
    std::vector<float> *          vp_vx;
    std::vector<float> *          vp_vy;
    std::vector<float> *          vp_vz;
    std::vector<int> *            vp_charge;
    std::vector<int> *            vp_pdgId;
    std::vector<int> *            vp_status;

    // Stub information
    std::vector<float> *          vb_x;
    std::vector<float> *          vb_y;
    std::vector<float> *          vb_z;
    std::vector<float> *          vb_r;
    std::vector<float> *          vb_eta;
    std::vector<float> *          vb_phi;
    std::vector<float> *          vb_coordx;
    std::vector<float> *          vb_coordy;
    std::vector<float> *          vb_trigBend;
    std::vector<float> *          vb_roughPt;
    std::vector<float> *          vb_clusWidth0;
    std::vector<float> *          vb_clusWidth1;
    std::vector<unsigned> *       vb_modId;
    std::vector<int> *            vb_tpId;

  protected:
    const int verbose_;
};

}  // namespace slhcl1tt

#endif
