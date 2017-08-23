#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/MCTruthAssociator.h"
using namespace slhcl1tt;

#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>

#include <TMath.h>

namespace {
// Comparator
bool sortByPt(const TrackingParticle& lhs, const TrackingParticle& rhs) {
    return (1.0/std::abs(lhs.invPt)) > (1.0/std::abs(rhs.invPt));
}

bool sortByLogicPt(const TTTrack2& ltk, const TTTrack2& rtk) {
    return (ltk.ndof() > rtk.ndof()) || (ltk.ndof() == rtk.ndof() && ltk.pt() > rtk.pt());
}

bool sortByMatchQuality(const std::pair<unsigned, float>& lhs, const std::pair<unsigned, float>& rhs) {
    return lhs.second < rhs.second;  // smaller is better
}

float absDiff(float lhs, float rhs) {
    return std::abs(lhs - rhs);
}

float squaredNormDiff(float lhs, float rhs, float scale) {
    return std::pow((lhs - rhs)/scale, 2.0);
}

float resolution(float qbpT, std::string trk_param){
  float Const=1, G_Const=1, G_Media=1, G_Sigma=1, r=1;

  if(trk_param == "invPt"){
    Const   =  0.001702;
    G_Const = -0.001517;
    G_Media =  0.000256;
    G_Sigma =  0.148572;
  }
  else if(trk_param == "phi0"){
    Const   =  0.000663;
    G_Const = -0.000559;
    G_Media = -0.000678;
    G_Sigma =  0.155241;
  }
  else if(trk_param == "z0"){
    Const   =  0.178339;
    G_Const = -0.099997;
    G_Media =  0.009426;
    G_Sigma =  0.838881;
  }
  else if(trk_param == "cottheta"){
    Const   =  0.003157;
    G_Const = -0.001072;
    G_Media =  0.005015;
    G_Sigma =  0.268481;
  }
  else
    std::cout<<"Parameter not recognized!"<<std::endl;

  r = Const + G_Const*TMath::Gaus(qbpT, G_Media, G_Sigma);
  return r;
}

  static const float degrees_of_freedom = 4.0;
  //static const float match_chi2_cut   = 16.0;
  //static const float match_chi2_cut = po_.maxSynChi2;
}

// _____________________________________________________________________________
//MCTruthAssociator::MCTruthAssociator() {
    //rms_invPt_    = 0.000165559;
    //rms_phi0_     = 7.29527e-05;
    //rms_cottheta_ = 0.00215212 ;
    //rms_z0_       = 0.0819977  ;
    //rms_d0_       = 1.00000    ;

    //rms_invPt_    = 0.002759;
    //rms_phi0_     = 0.000646;
    //rms_cottheta_ = 0.002695;
    //rms_z0_       = 0.08586 ;
    //rms_d0_       = 1.00000 ;

    //rms_invPt_    = 0.001283;
    //rms_phi0_     = 0.000453;
    //rms_cottheta_ = 0.002403;
    //rms_z0_       = 0.083406;
    //rms_d0_       = 1.000000;
//}

// _____________________________________________________________________________
void MCTruthAssociator::associate(std::vector<TrackingParticle>& trkParts, std::vector<TTTrack2>& tracks, float match_chi2_cut) {
    

    // Sort tracking particles by pT
    std::sort(trkParts.begin(), trkParts.end(), sortByPt);

    // Sort by Logic and pT combined
    std::sort(tracks.begin(), tracks.end(), sortByLogicPt);

    const unsigned nparts  = trkParts.size();
    const unsigned ntracks = tracks.size();

    //// Create the map for matching
    //std::map<unsigned, std::vector<std::pair<unsigned, float> > > matches;  // key=ipart, value=vector of (itrac
    //
    //// Loop on all tracking particles
    //for (unsigned ipart=0; ipart<nparts; ++ipart) {
    //    const TrackingParticle& trkPart = trkParts.at(ipart);
    //    assert(trkParts.at(ipart).tpId >= 0);
    //
    //    std::vector<std::pair<unsigned, float> >& ipart_matches = matches[ipart];
    //    assert(ipart_matches.size() == 0);
    //
    //    // Loop on all reconstructed tracks
    //    for (unsigned itrack=0; itrack<ntracks; ++itrack) {
    //        const TTTrack2& track = tracks.at(itrack);
    //        float quality = 999.;  // smaller is better
    //
    //        if (accept(trkPart, track, quality)) {
    //            ipart_matches.push_back(std::make_pair(itrack, quality));
    //        }
    //    }  // end loop on all reconstructed tracks
    //}  // end loop on all tracking particles
    //
    //// Debug
    //for (unsigned ipart=0; ipart<nparts; ++ipart) {
    //    const std::vector<std::pair<unsigned, float> >& ipart_matches = matches.at(ipart);
    //    for (unsigned imatch=0; imatch<ipart_matches.size(); ++imatch) {
    //        unsigned itrack = ipart_matches.at(imatch).first;
    //        float quality = ipart_matches.at(imatch).second;
    //        std::cout << "ipart: " << ipart << " itrack: " << itrack << " quality: " << quality << std::endl;
    //    }
    //}

    // Initiate each tracking particle as not found
    std::vector<int> mcCategories(nparts, ParticleCategory::NOTFOUND);

    // Loop on all tracking particles
    for (unsigned ipart=0; ipart<nparts; ++ipart) {

        bool foundTheBest = false;
        float best_quality = match_chi2_cut, quality = -999999.;
        int best_match_id = -1;

        // Loop on all reconstructed tracks
        for (unsigned itrack=0; itrack<ntracks; ++itrack) {

            // Compute our match chi2 definition for each combination and
            // checks if it is below our defined threshold (12.8)
	    bool accept_quality = accept(trkParts.at(ipart), tracks.at(itrack), quality, match_chi2_cut);
	    //if(abs(trkParts.at(ipart).pdgId) == 13){
	      //std::cout<<"trkParts_pdgId: "<<trkParts.at(ipart).pdgId<<", AMtrack: "<<itrack<<", roadRef: "<<tracks.at(itrack).roadRef()<<", combRef: "<<tracks.at(itrack).combRef()<<", MatchChi2: "<<quality<<std::endl;
	      //std::cout<<"trkPart: "<<trkParts.at(ipart)<<" -- AMtrack: "<<tracks.at(itrack)<<std::endl;
	    //}

            // If matchChi2 > 12.8 keeps the FAKE flag for AM track
            if (!accept_quality) continue;


            // If AM track is not marked as good track, then categorize the track
            if (tracks.at(itrack).synMatchCat() != TrackCategory::GOOD) {

                // Any AM track at this stage is below the match chi2 cut
                // then it will be good or duplicate
                tracks.at(itrack).setSynMatchChi2(quality);
                tracks.at(itrack).setSynMatchCat(TrackCategory::DUPLICATE);
                tracks.at(itrack).setSynTpId(trkParts.at(ipart).tpId);

                // Find the AM track that best matches with the real track
                if (quality < best_quality) {
                   best_match_id = itrack;
                   best_quality  = quality;
                   foundTheBest  = true;
                }
            }
        }  // end loop on all reconstructed tracks

        if (foundTheBest) {
           mcCategories.at(ipart) = ParticleCategory::FOUND;
           tracks.at(best_match_id).setSynMatchChi2(best_quality);
           tracks.at(best_match_id).setSynMatchCat(TrackCategory::GOOD);
           tracks.at(best_match_id).setSynTpId(trkParts.at(ipart).tpId);
        }

    }  // end loop on all tracking particles


    // Sanity check
    unsigned ngoods = 0, nduplicates = 0, nfakes = 0, nfounds = 0, nnotfounds = 0;

    for (unsigned ipart=0; ipart<nparts; ++ipart) {
        const int cat = mcCategories.at(ipart);
        if (cat == ParticleCategory::NOTFOUND)
            ++nnotfounds;
        else
            ++nfounds;
    }

    for (unsigned itrack=0; itrack<ntracks; ++itrack) {
        const int cat = tracks.at(itrack).synMatchCat();
        if (cat == TrackCategory::FAKE)
            ++nfakes;
        else if (cat == TrackCategory::DUPLICATE)
            ++nduplicates;
        else
            ++ngoods;
    }

    assert(nfounds + nnotfounds == nparts);
    assert(ngoods + nduplicates + nfakes == ntracks);
    assert(nfounds == ngoods);
}

// _____________________________________________________________________________
bool MCTruthAssociator::accept(const TrackingParticle& trkPart, const TTTrack2& track, float& quality, float match_chi2_cut) {
    quality = squaredNormDiff(trkPart.invPt   , track.invPt()   , resolution(track.invPt(),"invPt")    ) +
              squaredNormDiff(trkPart.phi0    , track.phi0()    , resolution(track.invPt(),"phi0")     ) +
              squaredNormDiff(trkPart.cottheta, track.cottheta(), resolution(track.invPt(),"cottheta") ) +
              squaredNormDiff(trkPart.z0      , track.z0()      , resolution(track.invPt(),"z0")       );
    quality /= degrees_of_freedom;
    bool acc = quality < match_chi2_cut;
    return acc;
}

// _____________________________________________________________________________
void MCTruthAssociator::print() {
    //std::cout << "rms invPt: " << rms_invPt_ << " phi0: " << rms_phi0_ << " cottheta: " << rms_cottheta_ << " z0: " << rms_z0_ << " d0: " << rms_d0_ << std::endl;
}
