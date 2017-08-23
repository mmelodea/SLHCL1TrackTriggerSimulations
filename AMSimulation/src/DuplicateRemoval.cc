#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/DuplicateRemoval.h"
using namespace slhcl1tt;

#include <vector>
#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <iostream>

namespace {
  bool sortByLogicPt(const TTTrack2& ltk, const TTTrack2& rtk) {
    return (ltk.ndof() > rtk.ndof()) || (ltk.ndof() == rtk.ndof() && ltk.pt() > rtk.pt());
  }

  bool sortByChi2(const TTTrack2& ltk, const TTTrack2& rtk) {
    return (ltk.chi2Red() < rtk.chi2Red());
    //int lchi2 = std::floor(ltk.chi2Red()/2.0);  // round to 2.0 step
    //int rchi2 = std::floor(rtk.chi2Red()/2.0);  // round to 2.0 step
    //return (lchi2 < rchi2) || (lchi2 == rchi2 && ltk.ndof() > rtk.ndof());
  }
}


void DuplicateRemoval::checkTracks(std::vector<TTTrack2>& all_tracks, int dupRm) {

  // Prevents wrong checking
  if (dupRm > 6) {
    throw std::invalid_argument("Number of stubs above maximum (6)");
  }

  // If setted, duplicate removal is done
  if (dupRm != -1) {
    // Sort AM tracks by logic and pt (decreasing)
    //std::sort(all_tracks.begin(), all_tracks.end(), sortByLogicPt);
    //std::sort(all_tracks.begin(), all_tracks.end(), sortByChi2);
    //std::stable_sort(all_tracks.begin(), all_tracks.end(), sortByChi2);

    // The duplicate removal itself
    std::vector<TTTrack2> unique_tracks;

    for (unsigned int itrack = 0; itrack < all_tracks.size(); itrack++) {
      //if (all_tracks.at(itrack).chi2()/all_tracks.at(itrack).ndof() > 3) std::cout<<"Violates Chi2/ndof<3 cut"<<std::endl;

      bool duplicate_found = false;
      //int duplicateTpId = -1;
      for (unsigned int jtrack = 0; jtrack < unique_tracks.size(); jtrack++) {
        // Sanity check
        const unsigned int nstubs = all_tracks.at(itrack).stubRefs().size();
        assert(nstubs == unique_tracks.at(jtrack).stubRefs().size());

        int sharedStubs = 0;
        for (unsigned int istub = 0; istub < nstubs; istub++) {
          // When comparing two 5/6's tracks we don't consider the pedestal value as a stub reference
          if (all_tracks.at(itrack).stubRefs()[istub] == 999999999 && unique_tracks.at(jtrack).stubRefs()[istub] == 999999999)
            continue;

          // Compares the stubs in each layer
          if (all_tracks.at(itrack).stubRefs()[istub] == unique_tracks.at(jtrack).stubRefs()[istub])
            sharedStubs++;
        }

        if (sharedStubs > dupRm) {
          duplicate_found = true;
	  //duplicateTpId = jtrack;
          break;
        }
      }  // loop over non-duplicate tracks (unique_tracks vector)

      // If track is not sharing more than allowed number of stubs, store it
      if (!duplicate_found)
        unique_tracks.push_back(all_tracks.at(itrack));


      
      //std::cout<<"Road: "<<all_tracks.at(itrack).roadRef()<<", Comb: "<<all_tracks.at(itrack).combRef();
      //for(int ist=0; ist<(int)all_tracks.at(itrack).stubRefs().size(); ++ist){
	//if( all_tracks.at(itrack).stubRefs()[ist] != 999999999 )
	// std::cout<<reader.vb_tpId->at( all_tracks.at(itrack).stubRefs()[ist] );
	//else std::cout<<"-1";

	//if(ist < 5) std::cout<<", ";
	//}
      //std::cout<<", FitterChi2/ndof: "<<all_tracks.at(itrack).chi2Red()<<", PattRef: "<<all_tracks.at(itrack).patternRef();
      //std::cout<<", DuplicateTpId: "<<duplicateTpId<<std::endl;
      

    }  // loop over all tracks

    // Remove the duplicate tracks from the original list
    // Keeps only the unique tracks
    std::swap(all_tracks, unique_tracks);

  } else {
    // The standard situation... no duplicate removal
    // So, nothing is done even calling the duplication removal member
  }

  return;
}
