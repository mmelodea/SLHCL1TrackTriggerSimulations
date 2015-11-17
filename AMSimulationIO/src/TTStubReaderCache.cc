#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTStubReaderCache.h"

//#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/Helper.h"
using namespace slhcl1tt;

TTStubReaderCache::TTStubReaderCache(const TTStubReader& reader, int verbose)
: vp_pt               (new std::vector<float>   (*reader.vp_pt)),
  vp_eta              (new std::vector<float>   (*reader.vp_eta)),
  vp_phi              (new std::vector<float>   (*reader.vp_phi)),
  vp_vx               (new std::vector<float>   (*reader.vp_vx)),
  vp_vy               (new std::vector<float>   (*reader.vp_vy)),
  vp_vz               (new std::vector<float>   (*reader.vp_vz)),
  vp_charge           (new std::vector<int>     (*reader.vp_charge)),
  vp_pdgId            (new std::vector<int>     (*reader.vp_pdgId)),
  vp_status           (new std::vector<int>     (*reader.vp_status)),
  //
  //vb_x                (new std::vector<float>   (*reader.vb_x)),
  //vb_y                (new std::vector<float>   (*reader.vb_y)),
  vb_z                (new std::vector<float>   (*reader.vb_z)),
  vb_r                (new std::vector<float>   (*reader.vb_r)),
  vb_eta              (new std::vector<float>   (*reader.vb_eta)),
  vb_phi              (new std::vector<float>   (*reader.vb_phi)),
  vb_coordx           (new std::vector<float>   (*reader.vb_coordx)),
  vb_coordy           (new std::vector<float>   (*reader.vb_coordy)),
  vb_trigBend         (new std::vector<float>   (*reader.vb_trigBend)),
  //vb_roughPt          (new std::vector<float>   (*reader.vb_roughPt)),
  //vb_clusWidth0       (new std::vector<float>   (*reader.vb_clusWidth0)),
  //vb_clusWidth1       (new std::vector<float>   (*reader.vb_clusWidth1)),
  vb_modId            (new std::vector<unsigned>(*reader.vb_modId)),
  vb_tpId             (new std::vector<int>     (*reader.vb_tpId)),
  //
  verbose_(verbose) {}

TTStubReaderCache::~TTStubReaderCache() {
    delete vp_pt;
    delete vp_eta;
    delete vp_phi;
    delete vp_vx;
    delete vp_vy;
    delete vp_vz;
    delete vp_charge;
    delete vp_pdgId;
    delete vp_status;
    //
    //delete vb_x;
    //delete vb_y;
    delete vb_z;
    delete vb_r;
    delete vb_eta;
    delete vb_phi;
    delete vb_coordx;
    delete vb_coordy;
    delete vb_trigBend;
    //delete vb_roughPt;
    //delete vb_clusWidth0;
    //delete vb_clusWidth1;
    delete vb_modId;
    delete vb_tpId;
}
