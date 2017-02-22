#ifndef __TTTrackReader_h__
#define __TTTrackReader_h__

#include "TChain.h"
#include "TFile.h"
#include "TFileCollection.h"
#include "TString.h"
#include "TTree.h"
#include <memory>
#include <vector>

#include "BasicReader.h"

// _____________________________________________________________________________
// This is a simple wrapper around TChain. It sets the branch names, addresses
// etc. Its functions are essentially the same as the functions of TChain.

template<size_t T>
class TTTrackReaderT {
public:
    TTTrackReaderT(int verbose=1);
    ~TTTrackReaderT();

    void init(TString src);

    Long64_t loadTree(Long64_t entry) { return tchain_->LoadTree(entry); }

    Int_t getEntry(Long64_t entry) { return tchain_->GetEntry(entry); }

    Long64_t getEntries() const { return tchain_->GetEntries(); }

    TChain* getChain() { return tchain_; }

    template <typename T2>
    void nullVectorElements(std::vector<T2>* v, const std::vector<bool>& nulling) {
        assert(v->size() == nulling.size());
        for (unsigned i=0; i<nulling.size(); ++i)
            if (nulling[i])
                (*v)[i] = 0;
    }

    void nullStubs(const std::vector<bool>& nulling);

    void nullParticles(const std::vector<bool>& nulling);

    // genParticle information
    std::vector<float> *          vp_pt;        // particle pT [GeV]
    std::vector<float> *          vp_eta;       // particle eta
    std::vector<float> *          vp_phi;       // particle phi [rad]
    std::vector<float> *          vp_vx;        // particle vertex, x [cm]
    std::vector<float> *          vp_vy;        // particle vertex, y [cm]
    std::vector<float> *          vp_vz;        // particle vertex, z [cm]
    std::vector<int> *            vp_charge;    // particle charge

    // genJet information
    std::vector<float> *          vj_pt;        // jet pT [GeV]
    std::vector<float> *          vj_eta;       // jet eta
    std::vector<float> *          vj_phi;       // jet phi [rad]

    // trkParticle information
    std::vector<float> *          vp2_pt;
    std::vector<float> *          vp2_eta;
    std::vector<float> *          vp2_phi;
    std::vector<float> *          vp2_vx;
    std::vector<float> *          vp2_vy;
    std::vector<float> *          vp2_vz;
    std::vector<int> *            vp2_charge;
    std::vector<int> *            vp2_pdgId;    // particle pdgID
    std::vector<bool> *           vp2_signal;   // is particle from hard scattering vertex? (vertex 0)
    std::vector<bool> *           vp2_intime;   // is particle from BX 0?
    std::vector<bool> *           vp2_primary;  // is particle produced at a primary vertex?

    // Stub information
    std::vector<float> *          vb_z;             // stub global z [cm]
    std::vector<float> *          vb_r;             // stub global r [cm]
    std::vector<float> *          vb_eta;           // stub global eta
    std::vector<float> *          vb_phi;           // stub global phi [rad]
    std::vector<float> *          vb_coordx;        // stub local phi
    std::vector<float> *          vb_coordy;        // stub local z (barrel) or r (endcap)
    std::vector<float> *          vb_trigBend;      // stub local bend distance (a.k.a. delta-s)
    std::vector<unsigned> *       vb_modId;         // stub module ID
    std::vector<unsigned> *       vb_clusRef0;      // stub inner cluster index
    std::vector<unsigned> *       vb_clusRef1;      // stub outer cluster index
    std::vector<int> *            vb_tpId;          // stub tracking particle ID
    std::vector<std::string> *    vb_bitString;     // stub info encoded as in firmware
    std::vector<unsigned> *       vb_superstripId;  // stub superstrip ID

    // Road information
    std::vector<unsigned> *                             vr_patternRef;    // road pattern id
    std::vector<unsigned> *                             vr_tower;         // road tower
    std::vector<unsigned> *                             vr_nstubs;        // road total number of stubs
    std::vector<float> *                                vr_patternInvPt;  // road pattern q/pT [1/GeV]
    std::vector<unsigned> *                             vr_patternFreq;   // road pattern frequency
    std::vector<std::vector<unsigned> > *               vr_superstripIds; // road superstrip ids
    std::vector<std::vector<std::vector<unsigned> > > * vr_stubRefs;      // road stub ids

    // Track information
    //std::vector<float> *                  vt_px;          // track momentum [GeV]
    //std::vector<float> *                  vt_py;          // track momentum [GeV]
    //std::vector<float> *                  vt_pz;          // track momentum [GeV]
    std::vector<float> *                  vt_pt;          // track momentum [GeV]
    std::vector<float> *                  vt_eta;         // track momentum [GeV]
    //std::vector<float> *                  vt_phi;         // track momentum [GeV]
    //std::vector<float> *                  vt_vx;          // track vertex [cm]
    //std::vector<float> *                  vt_vy;          // track vertex [cm]
    //std::vector<float> *                  vt_vz;          // track vertex [cm]
    std::vector<float> *                  vt_rinv;        // = 0.003 * 3.8 * track fitted q/pT
    std::vector<float> *                  vt_invPt;       // track q/pT [1/GeV]
    std::vector<float> *                  vt_phi0;        // track phi at r=0
    std::vector<float> *                  vt_cottheta;    // track cot theta at r=0
    std::vector<float> *                  vt_z0;          // track vertex z [cm]
    std::vector<float> *                  vt_d0;          // track impact parameter [cm]
    std::vector<float> *                  vt_chi2;        // track fit chi-squared
    std::vector<int> *                    vt_ndof;        // track fit num of degrees of freedom
    std::vector<float> *                  vt_synMatchChi2;// track MC match chi-squared
    std::vector<int> *                    vt_synMatchCat; // track MC match category
    std::vector<int> *                    vt_tpId;        // associated tracking particle id
    std::vector<int> *                    vt_synTpId;     // associated tracking particle id using 'synthetic' approach
    std::vector<unsigned> *               vt_tower;       // track tower
    std::vector<unsigned> *               vt_hitBits;     // track stub composition by layers
    std::vector<unsigned> *               vt_ptSegment;   // obsolete
    std::vector<unsigned> *               vt_roadRef;     // track road id
    std::vector<unsigned> *               vt_combRef;     // track combination id
    std::vector<unsigned> *               vt_patternRef;  // track pattern id
    std::vector<std::vector<unsigned> > * vt_stubRefs;    // track stub ids
    std::vector<std::vector<float> > *    vt_principals;  // track principle components

protected:
    TChain* tchain_;
    int treenumber_;
    const int verbose_;
};

// _____________________________________________________________________________
// Implementation is included in the header file to simplify ROOT library generation

#include <cassert>
#include <stdexcept>

template<size_t T>
TTTrackReaderT<T>::TTTrackReaderT(int verbose) :
    vp_pt               (0),
    vp_eta              (0),
    vp_phi              (0),
    vp_vx               (0),
    vp_vy               (0),
    vp_vz               (0),
    vp_charge           (0),
    //
    vj_pt               (0),
    vj_eta              (0),
    vj_phi              (0),
    //
    vp2_pt              (0),
    vp2_eta             (0),
    vp2_phi             (0),
    vp2_vx              (0),
    vp2_vy              (0),
    vp2_vz              (0),
    vp2_charge          (0),
    vp2_pdgId           (0),
    vp2_signal          (0),
    vp2_intime          (0),
    vp2_primary         (0),
    //
    vb_z                (0),
    vb_r                (0),
    vb_eta              (0),
    vb_phi              (0),
    vb_coordx           (0),
    vb_coordy           (0),
    vb_trigBend         (0),
    vb_modId            (0),
    vb_clusRef0         (0),
    vb_clusRef1         (0),
    vb_tpId             (0),
    vb_bitString        (0),
    vb_superstripId     (0),
    //
    vr_patternRef       (0),
    vr_tower            (0),
    vr_nstubs           (0),
    vr_patternInvPt     (0),
    vr_patternFreq      (0),
    vr_superstripIds    (0),
    vr_stubRefs         (0),
    //
    //vt_px               (0),
    //vt_py               (0),
    //vt_pz               (0),
    vt_pt               (0),
    vt_eta              (0),
    //vt_phi              (0),
    //vt_vx               (0),
    //vt_vy               (0),
    //vt_vz               (0),
    vt_rinv             (0),
    vt_invPt            (0),
    vt_phi0             (0),
    vt_cottheta         (0),
    vt_z0               (0),
    vt_d0               (0),
    vt_chi2             (0),
    vt_ndof             (0),
    vt_synMatchChi2     (0),
    vt_synMatchCat      (0),
    vt_tpId             (0),
    vt_synTpId          (0),
    vt_tower            (0),
    vt_hitBits          (0),
    vt_ptSegment        (0),
    vt_roadRef          (0),
    vt_combRef          (0),
    vt_patternRef       (0),
    vt_stubRefs         (0),
    vt_principals       (0),
    //
    tchain_(0), treenumber_(0), verbose_(verbose) {}

template<size_t T>
TTTrackReaderT<T>::~TTTrackReaderT() {
    if (tchain_)  delete tchain_;
}

template<size_t T>
void TTTrackReaderT<T>::init(TString src) {
    if (!src.EndsWith(".root") && !src.EndsWith(".txt")) {
        TString msg = "Input source must be either .root or .txt";
        throw std::invalid_argument(msg.Data());
    }

    //if (verbose_)  std::cout << "Opening " << src << std::endl;
    tchain_ = new TChain("ntupler/tree");

    if (src.EndsWith(".root")) {
        if (tchain_->Add(src) ) {
            //if (verbose_)  std::cout << "Successfully read " << src << std::endl;
        } else {
            TString msg = "Failed to read " + src;
            throw std::invalid_argument(msg.Data());
        }

    } else if (src.EndsWith(".txt")) {
        TFileCollection fc("fileinfolist", "", src);
        if (tchain_->AddFileInfoList((TCollection*) fc.GetList()) ) {
            //if (verbose_)  std::cout << "Successfully read " << src << std::endl;
        } else {
            TString msg = "Failed to read " + src;
            throw std::invalid_argument(msg.Data());
        }
    }

    TString prefix = "", suffix = "";

    assert(tchain_ != 0);
    treenumber_ = tchain_->GetTreeNumber();

    if (T == kPatternMatcher) {
        tchain_->SetBranchStatus("*"                 , 0);
        tchain_->SetBranchStatus("genParts_pt"       , 1);
        tchain_->SetBranchStatus("genParts_eta"      , 1);
        tchain_->SetBranchStatus("genParts_phi"      , 1);
        tchain_->SetBranchStatus("genParts_vx"       , 1);
        tchain_->SetBranchStatus("genParts_vy"       , 1);
        tchain_->SetBranchStatus("genParts_vz"       , 1);
        tchain_->SetBranchStatus("genParts_charge"   , 1);
        //
        tchain_->SetBranchStatus("genJets_pt"        , 1);
        tchain_->SetBranchStatus("genJets_eta"       , 1);
        tchain_->SetBranchStatus("genJets_phi"       , 1);
        //
        tchain_->SetBranchStatus("trkParts_pt"       , 1);
        tchain_->SetBranchStatus("trkParts_eta"      , 1);
        tchain_->SetBranchStatus("trkParts_phi"      , 1);
        tchain_->SetBranchStatus("trkParts_vx"       , 1);
        tchain_->SetBranchStatus("trkParts_vy"       , 1);
        tchain_->SetBranchStatus("trkParts_vz"       , 1);
        tchain_->SetBranchStatus("trkParts_charge"   , 1);
        tchain_->SetBranchStatus("trkParts_pdgId"    , 1);
        tchain_->SetBranchStatus("trkParts_signal"   , 1);
        tchain_->SetBranchStatus("trkParts_intime"   , 1);
        tchain_->SetBranchStatus("trkParts_primary"  , 1);
        //
        tchain_->SetBranchStatus("TTStubs_z"         , 1);
        tchain_->SetBranchStatus("TTStubs_r"         , 1);
        tchain_->SetBranchStatus("TTStubs_eta"       , 1);
        tchain_->SetBranchStatus("TTStubs_phi"       , 1);
        tchain_->SetBranchStatus("TTStubs_coordx"    , 1);
        tchain_->SetBranchStatus("TTStubs_coordy"    , 1);
        tchain_->SetBranchStatus("TTStubs_trigBend"  , 1);
        tchain_->SetBranchStatus("TTStubs_modId"     , 1);
        tchain_->SetBranchStatus("TTStubs_clusRef0"  , 1);
        tchain_->SetBranchStatus("TTStubs_clusRef1"  , 1);
        tchain_->SetBranchStatus("TTStubs_tpId"      , 1);
    }

    if (T == kFuture) {
        tchain_->SetBranchAddress("genParts_pt"       , &(vp_pt));
        tchain_->SetBranchAddress("genParts_eta"      , &(vp_eta));
        tchain_->SetBranchAddress("genParts_phi"      , &(vp_phi));
        tchain_->SetBranchAddress("genParts_vx"       , &(vp_vx));
        tchain_->SetBranchAddress("genParts_vy"       , &(vp_vy));
        tchain_->SetBranchAddress("genParts_vz"       , &(vp_vz));
        tchain_->SetBranchAddress("genParts_charge"   , &(vp_charge));
    }
    //
    if (T == kFuture) {
        tchain_->SetBranchAddress("genJets_pt"        , &(vj_pt));
        tchain_->SetBranchAddress("genJets_eta"       , &(vj_eta));
        tchain_->SetBranchAddress("genJets_phi"       , &(vj_phi));
    }
    //
    if (T == kPatternMatcher || T == kTrackFitter || T == kFuture) {
        tchain_->SetBranchAddress("trkParts_pt"       , &(vp2_pt));
        tchain_->SetBranchAddress("trkParts_eta"      , &(vp2_eta));
        tchain_->SetBranchAddress("trkParts_phi"      , &(vp2_phi));
        tchain_->SetBranchAddress("trkParts_vx"       , &(vp2_vx));
        tchain_->SetBranchAddress("trkParts_vy"       , &(vp2_vy));
        tchain_->SetBranchAddress("trkParts_vz"       , &(vp2_vz));
        tchain_->SetBranchAddress("trkParts_charge"   , &(vp2_charge));
        tchain_->SetBranchAddress("trkParts_pdgId"    , &(vp2_pdgId));
        tchain_->SetBranchAddress("trkParts_signal"   , &(vp2_signal));
        tchain_->SetBranchAddress("trkParts_intime"   , &(vp2_intime));
        tchain_->SetBranchAddress("trkParts_primary"  , &(vp2_primary));
    }
    //
    if (T == kPatternMatcher || T == kTrackFitter || T == kFuture) {
        tchain_->SetBranchAddress("TTStubs_z"         , &(vb_z));
        tchain_->SetBranchAddress("TTStubs_r"         , &(vb_r));
        tchain_->SetBranchAddress("TTStubs_eta"       , &(vb_eta));
        tchain_->SetBranchAddress("TTStubs_phi"       , &(vb_phi));
        tchain_->SetBranchAddress("TTStubs_coordx"    , &(vb_coordx));
        tchain_->SetBranchAddress("TTStubs_coordy"    , &(vb_coordy));
        tchain_->SetBranchAddress("TTStubs_trigBend"  , &(vb_trigBend));
        tchain_->SetBranchAddress("TTStubs_modId"     , &(vb_modId));
        tchain_->SetBranchAddress("TTStubs_clusRef0"  , &(vb_clusRef0));
        tchain_->SetBranchAddress("TTStubs_clusRef1"  , &(vb_clusRef1));
        tchain_->SetBranchAddress("TTStubs_tpId"      , &(vb_tpId));
    }
    //
    if (T == kTrackFitter || T == kFuture) {
        tchain_->SetBranchAddress("TTStubs_bitString"   , &(vb_bitString));
        tchain_->SetBranchAddress("TTStubs_superstripId", &(vb_superstripId));
    }
    //
    if (T == kRoadMerging || T == kTrackFitter || T == kFuture) {
        prefix = "AMTTRoads_";
        tchain_->SetBranchAddress(prefix + "patternRef"    + suffix, &(vr_patternRef));
        tchain_->SetBranchAddress(prefix + "tower"         + suffix, &(vr_tower));
        tchain_->SetBranchAddress(prefix + "nstubs"        + suffix, &(vr_nstubs));
        tchain_->SetBranchAddress(prefix + "patternInvPt"  + suffix, &(vr_patternInvPt));
        tchain_->SetBranchAddress(prefix + "patternFreq"   + suffix, &(vr_patternFreq));
        tchain_->SetBranchAddress(prefix + "superstripIds" + suffix, &(vr_superstripIds));
        tchain_->SetBranchAddress(prefix + "stubRefs"      + suffix, &(vr_stubRefs));
    }
    //
    if (T == kFuture) {
        prefix = "AMTTTracks_";
        //tchain_->SetBranchAddress(prefix + "px"           + suffix, &(vt_px));
        //tchain_->SetBranchAddress(prefix + "py"           + suffix, &(vt_py));
        //tchain_->SetBranchAddress(prefix + "pz"           + suffix, &(vt_pz));
        tchain_->SetBranchAddress(prefix + "pt"           + suffix, &(vt_pt));
        tchain_->SetBranchAddress(prefix + "eta"          + suffix, &(vt_eta));
        //tchain_->SetBranchAddress(prefix + "phi"          + suffix, &(vt_phi));
        //tchain_->SetBranchAddress(prefix + "vx"           + suffix, &(vt_vx));
        //tchain_->SetBranchAddress(prefix + "vy"           + suffix, &(vt_vy));
        //tchain_->SetBranchAddress(prefix + "vz"           + suffix, &(vt_vz));
        tchain_->SetBranchAddress(prefix + "rinv"         + suffix, &(vt_rinv));
        tchain_->SetBranchAddress(prefix + "invPt"        + suffix, &(vt_invPt));
        tchain_->SetBranchAddress(prefix + "phi0"         + suffix, &(vt_phi0));
        tchain_->SetBranchAddress(prefix + "cottheta"     + suffix, &(vt_cottheta));
        tchain_->SetBranchAddress(prefix + "z0"           + suffix, &(vt_z0));
        tchain_->SetBranchAddress(prefix + "d0"           + suffix, &(vt_d0));
        tchain_->SetBranchAddress(prefix + "chi2"         + suffix, &(vt_chi2));
        tchain_->SetBranchAddress(prefix + "ndof"         + suffix, &(vt_ndof));
        tchain_->SetBranchAddress(prefix + "synMatchChi2" + suffix, &(vt_synMatchChi2));
        tchain_->SetBranchAddress(prefix + "synMatchCat"  + suffix, &(vt_synMatchCat));
        tchain_->SetBranchAddress(prefix + "tpId"         + suffix, &(vt_tpId));
        tchain_->SetBranchAddress(prefix + "synTpId"      + suffix, &(vt_synTpId));
        tchain_->SetBranchAddress(prefix + "tower"        + suffix, &(vt_tower));
        tchain_->SetBranchAddress(prefix + "hitBits"      + suffix, &(vt_hitBits));
        tchain_->SetBranchAddress(prefix + "ptSegment"    + suffix, &(vt_ptSegment));
        tchain_->SetBranchAddress(prefix + "roadRef"      + suffix, &(vt_roadRef));
        tchain_->SetBranchAddress(prefix + "combRef"      + suffix, &(vt_combRef));
        tchain_->SetBranchAddress(prefix + "patternRef"   + suffix, &(vt_patternRef));
        tchain_->SetBranchAddress(prefix + "stubRefs"     + suffix, &(vt_stubRefs));
        tchain_->SetBranchAddress(prefix + "principals"   + suffix, &(vt_principals));
    }
}

template<size_t T>
void TTTrackReaderT<T>::nullStubs(const std::vector<bool>& nulling) {
    nullVectorElements(vb_z         , nulling);
    nullVectorElements(vb_r         , nulling);
    nullVectorElements(vb_eta       , nulling);
    nullVectorElements(vb_phi       , nulling);
    nullVectorElements(vb_coordx    , nulling);
    nullVectorElements(vb_coordy    , nulling);
    nullVectorElements(vb_trigBend  , nulling);
    //nullVectorElements(vb_modId     , nulling);  // don't null this guy
    nullVectorElements(vb_clusRef0  , nulling);
    nullVectorElements(vb_clusRef1  , nulling);
    nullVectorElements(vb_tpId      , nulling);
}

template<size_t T>
void TTTrackReaderT<T>::nullParticles(const std::vector<bool>& nulling) {
    nullVectorElements(vp2_pt        , nulling);
    nullVectorElements(vp2_eta       , nulling);
    nullVectorElements(vp2_phi       , nulling);
    nullVectorElements(vp2_vx        , nulling);
    nullVectorElements(vp2_vy        , nulling);
    nullVectorElements(vp2_vz        , nulling);
    //nullVectorElements(vp2_charge    , nulling);  // don't null this guy
    nullVectorElements(vp2_pdgId     , nulling);
    nullVectorElements(vp2_signal    , nulling);
    nullVectorElements(vp2_intime    , nulling);
    nullVectorElements(vp2_primary   , nulling);
}

#endif  // __TTTrackReader_h__
