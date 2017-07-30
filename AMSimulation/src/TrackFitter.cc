#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/TrackFitter.h"

namespace {
unsigned getPtSegment(float invPt) {  // for PCA
    return (invPt - PCA_MIN_INVPT) / (PCA_MAX_INVPT - PCA_MIN_INVPT) * PCA_NSEGMENTS;
}

unsigned getHitBits(const std::vector<bool>& stubs_bool) {
    unsigned bitset = 0;

    for (unsigned i=0; i<stubs_bool.size(); ++i) {
        bitset |= (stubs_bool.at(i) << i);
    }

    switch (bitset) {
    case 0b111111:  return 0;
    case 0b111110:  return 1;
    case 0b111101:  return 2;
    case 0b111011:  return 3;
    case 0b110111:  return 4;
    case 0b101111:  return 5;
    case 0b011111:  return 6;
    default      :  return 7;
    }
}

// Comparator
bool sortByPt(const TTTrack2& lhs, const TTTrack2& rhs) {
    return lhs.pt() > rhs.pt();
}

// principal component cut function
bool PrincipalCuts(std::vector<float> princes){
  bool pass=true;
  const float Cuts6[12]={58.,33.,22.,12.,-1.,-1.,9.,3.,3.,3.,-1.,-1.};
  const float Cuts5[10]={33.,22.,12.,-1.,-1.,3.,3.,3.,-1.,-1.};
  if(princes.size()==12) {
    for(unsigned i=0; i<princes.size(); ++i){
      if(Cuts6[i]==-1) continue;
      if(fabs(princes[i])>Cuts6[i]){
        pass=false;
        break;
      }
    }
  } else {
    for(unsigned i=0; i<princes.size(); ++i){
      if(Cuts5[i]==-1) continue;
      if(fabs(princes[i])>Cuts5[i]){
        pass=false;
        break;
      }
    }
  }
  return pass;
}
}


// _____________________________________________________________________________
// Do track fitting
int TrackFitter::makeTracks(TString src, TString out) {
    if (verbose_)  std::cout << Info() << "Reading " << nEvents_ << " events and fitting tracks." << std::endl;

    // _________________________________________________________________________
    // For reading
    TTRoadReader reader(verbose_);
    reader.init(src);

    // _________________________________________________________________________
    // For writing
    TTTrackWriter writer(verbose_);
    writer.init(reader.getChain(), out);

    // _________________________________________________________________________
    // Loop over all events

    // Containers
    std::vector<TTTrack2> tracks;
    tracks.reserve(300);

    // Bookkeepers
    long int nRead = 0, nKept = 0;

    for (long long ievt=0; ievt<nEvents_; ++ievt) {
        if (reader.loadTree(ievt) < 0)  break;
        reader.getEntry(ievt);

        const unsigned nroads = reader.vr_patternRef->size();
        if (verbose_>1 && ievt%100==0)  std::cout << Debug() << Form("... Processing event: %7lld, fitting: %7ld", ievt, nKept) << std::endl;
        if (verbose_>2)  std::cout << Debug() << "... evt: " << ievt << " # roads: " << nroads << std::endl;

        if (!nroads) {  // skip if no road
            writer.fillTracks(std::vector<TTTrack2>());
            ++nRead;
            continue;
        }

        tracks.clear();
        int fitstatus = 0;


        // _____________________________________________________________________
        // Track fitters taking fit combinations

        // Loop over the roads
        for (unsigned iroad=0; iroad<nroads; ++iroad) {
            if (iroad >= (unsigned) po_.maxRoads)  break;

            const unsigned patternRef = reader.vr_patternRef->at(iroad);
            if (patternRef >= (unsigned) po_.maxPatterns)  continue;

            // Apply Hough Transform filter
            std::vector<std::vector<unsigned> > stubRefs;
            const HTMemory::Candidate& htm_cand = htm_->getBestCandidate_amsim(po_.htmconf, reader.vr_stubRefs->at(iroad), *reader.vb_z, stubRefs);
            if (htm_cand.nTotalHits == 0)  continue;

            // Get combinations of stubRefs
            //std::vector<std::vector<unsigned> > stubRefs = reader.vr_stubRefs->at(iroad);
            std::vector<std::vector<float> > stubDeltaS(stubRefs.size(), std::vector<float>());
            for (unsigned ilayer=0; ilayer<stubRefs.size(); ++ilayer) {
                for(unsigned istub=0; istub<stubRefs.at(ilayer).size(); ++istub) {
                    const unsigned stubRef = stubRefs.at(ilayer).at(istub);
                    stubDeltaS.at(ilayer).push_back(reader.vb_trigBend->at(stubRef));
                }
            }

            // Choose either the normal combination building or the 5/6 permutations per 6/6 road in addition
            // and/or pairwise Delta Delta S cleaning (PDDS)
            // Quote from Marco: The compiler will likely do RVO so the move might not be needed, we prefer to be explicit about it.
            std::vector<std::vector<unsigned> > combinations;
            if (po_.oldCB)
                combinations = std::move(combinationFactory_.combine(stubRefs, po_.FiveOfSix));
            else if (po_.PDDS)
                combinations = std::move(pairCombinationFactory_.combine(stubRefs, stubDeltaS, po_.FiveOfSix));
            else
                combinations = std::move(combinationBuilderFactory_->combine(stubRefs));

            assert(combinations.size() > 0);
            for (unsigned icomb=0; icomb<combinations.size(); ++icomb)
                assert(combinations.at(icomb).size() == reader.vr_stubRefs->at(iroad).size());

            if (verbose_>2) {
                std::cout << Debug() << "... ... road: " << iroad << " # combinations: " << combinations.size() << std::endl;
            }

            // Loop over the combinations
            for (unsigned icomb=0; icomb<combinations.size(); ++icomb) {
                if (icomb >= (unsigned) po_.maxCombs)  break;

                // Create and set TTRoadComb
                TTRoadComb acomb;
                acomb.roadRef    = iroad;
                acomb.combRef    = icomb;
                acomb.patternRef = patternRef;
                acomb.ptSegment  = getPtSegment(reader.vr_patternInvPt->at(iroad));
                acomb.stubRefs   = combinations.at(icomb);

                acomb.stubs_r   .clear();
                acomb.stubs_phi .clear();
                acomb.stubs_z   .clear();
                acomb.stubs_bool.clear();

                for (unsigned istub=0; istub<acomb.stubRefs.size(); ++istub) {
                    const unsigned stubRef = acomb.stubRefs.at(istub);
                    if (stubRef != CombinationFactory::BAD) {
                        acomb.stubs_r   .push_back(reader.vb_r   ->at(stubRef));
                        acomb.stubs_phi .push_back(reader.vb_phi ->at(stubRef));
                        acomb.stubs_z   .push_back(reader.vb_z   ->at(stubRef));
                        acomb.stubs_bool.push_back(true);
                    } else {
                        acomb.stubs_r   .push_back(0.);
                        acomb.stubs_phi .push_back(0.);
                        acomb.stubs_z   .push_back(0.);
                        acomb.stubs_bool.push_back(false);
                    }
                }

                acomb.hitBits = getHitBits(acomb.stubs_bool);

                if (verbose_>2) {
                    std::cout << Debug() << "... ... ... comb: " << icomb << " " << acomb;
                    std::cout << std::endl;
                }

                // _____________________________________________________________
                // Fit
                TTTrack2 atrack;
                fitstatus = fitter_->fit(acomb, atrack);

                atrack.setTower     (po_.tower);
                atrack.setRoadRef   (acomb.roadRef);
                atrack.setCombRef   (acomb.combRef);
                atrack.setPatternRef(acomb.patternRef);
                atrack.setPtSegment (acomb.ptSegment);
                atrack.setHitBits   (acomb.hitBits);
                atrack.setStubRefs  (acomb.stubRefs);

                if (!po_.CutPrincipals) {
                    if (atrack.chi2Red() < po_.maxChi2)  // reduced chi^2 = chi^2 / ndof
                        tracks.push_back(atrack);
                } else {
                    if (PrincipalCuts(atrack.principals()))
                        tracks.push_back(atrack);
                }

                if (verbose_>2)  std::cout << Debug() << "... ... ... track: " << icomb << " status: " << fitstatus << " reduced chi2: " << atrack.chi2Red() << " invPt: " << atrack.invPt() << " phi0: " << atrack.phi0() << " cottheta: " << atrack.cottheta() << " z0: " << atrack.z0() << std::endl;
            }
        }  // loop over the roads

        std::sort(tracks.begin(), tracks.end(), sortByPt);


        if (verbose_>2)  std::cout << Debug() << "... evt: " << ievt << " # tracks: " << tracks.size() << std::endl;
        if (verbose_>3) {
            for (unsigned itrack=0; itrack!=tracks.size(); ++itrack) {
                std::cout << "... ... track: " << itrack << " " << tracks.at(itrack) << std::endl;
            }
        }

        if (! tracks.empty())
            ++nKept;

        if (tracks.size() > (unsigned) po_.maxTracks)
            tracks.resize(po_.maxTracks);


        // ---------------------------------------------------------------------
        // Classify tracks as duplicates or not (for duplicate removal)
        // In the algorithm, AM tracks are sorted by chi2 (used to be matching logic and pT)
        // ---------------------------------------------------------------------
        DuplicateRemoval flagDuplicates;
        if (po_.rmDuplicate != -1) flagDuplicates.checkTracks(tracks, po_.rmDuplicate);

        //----------------------------------------------------------------------
        // Identify and flag duplicates by defining a track-parameter space
        // inside of which anything is considered to be a single track
        //----------------------------------------------------------------------
        ParameterDuplicateRemoval removeParameterDuplicates;
        if (po_.rmParDuplicate) removeParameterDuplicates.ReduceTracks(tracks);


        // _____________________________________________________________________
        // Track categorization

        if (po_.speedup<1) {
            const unsigned nparts = reader.vp2_primary->size();
            if (verbose_>2)  std::cout << Debug() << "... evt: " << ievt << " # particles: " << nparts << std::endl;

            std::vector<TrackingParticle> trkParts;
            for (unsigned ipart=0; ipart<nparts; ++ipart) {
                bool  primary         = reader.vp2_primary->at(ipart);
                bool  intime          = reader.vp2_intime->at(ipart);
                int   simCharge       = reader.vp2_charge->at(ipart);
                float simPt           = reader.vp2_pt->at(ipart);

                if (simCharge!=0 && primary && intime && simPt>=1) {
                    float simEta          = reader.vp2_eta->at(ipart);
                    float simPhi          = reader.vp2_phi->at(ipart);
                    //float simVx           = reader.vp2_vx->at(ipart);
                    //float simVy           = reader.vp2_vy->at(ipart);
                    float simVz           = reader.vp2_vz->at(ipart);
                    int   simCharge       = reader.vp2_charge->at(ipart);
                    int   simPdgId        = reader.vp2_pdgId->at(ipart);

                    float simCotTheta     = std::sinh(simEta);
                    float simChargeOverPt = float(simCharge)/simPt;

                    trkParts.emplace_back(TrackingParticle{  // using POD type constructor
                        (int) ipart,
                        simPdgId,
                        simChargeOverPt,
                        simPhi,
                        simCotTheta,
                        simVz,
                        0.
                    });

                    if (verbose_>3)  std::cout << Debug() << "... ... part: " << ipart << " primary: " << primary << " " << trkParts.back();
                }
            }
            truthAssociator_.associate(trkParts, tracks);
        }

        writer.fillTracks(tracks);
        ++nRead;
    }

    if (nRead == 0) {
        std::cout << Error() << "Failed to read any event." << std::endl;
        return 1;
    }

    if (verbose_)  std::cout << Info() << Form("Read: %7ld, triggered: %7ld", nRead, nKept) << std::endl;


    // _________________________________________________________________________
    // Write histograms

    for (std::map<TString, TH1F *>::const_iterator it=fitter_->histograms_.begin();
         it!=fitter_->histograms_.end(); ++it) {
        if (it->second)  it->second->SetDirectory(gDirectory);
    }

    writer.write();

    return 0;
}


// _____________________________________________________________________________
// Main driver
int TrackFitter::run() {
    int exitcode = 0;
    Timing(1);

    exitcode = makeTracks(po_.input, po_.output);
    if (exitcode)  return exitcode;
    Timing();

    return exitcode;
}
