#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/StubCleanerPU.h"

#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTStubReader.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTStubReaderCache.h"

static const unsigned MIN_NGOODSTUBS = 3;
static const unsigned MAX_NGOODSTUBS = 8;

namespace {
// Comparator
bool sortByFloat(const std::pair<unsigned, float>& lhs, const std::pair<unsigned, float>& rhs) {
    return lhs.second < rhs.second;
}

// Comparator
bool sortByUnsignedThenFloat(const std::pair<unsigned, std::pair<unsigned, float> >& lhs,
                             const std::pair<unsigned, std::pair<unsigned, float> >& rhs) {
    // Primary condition
    if (lhs.second.first < rhs.second.first)  return true;
    if (lhs.second.first > rhs.second.first)  return false;

    // Secondary condition
    return lhs.second.second < rhs.second.second;
}

// Insert an element into a specific position in a sorted list
// The algorithm moves every element behind that specific location from i to
// i+1, creating a hole at that location to accomodate the new element.
// 'first' points to the first element in the list
// 'len' is the length of the list
// 'pos' is where the new element is going to be inserted
// 'value' is the value of the new element
template<typename RandomAccessIterator, typename Size, typename T>
void insertSorted(RandomAccessIterator first, Size len, Size pos, const T value) {
    first += len;
    len -= pos;
    while (len>0) {
        *first = std::move(*(first-1));
        --first;
        --len;
    }
    *first = std::move(value);
}

float calcIdealPhi(float simPhi, float simChargeOverPt, float r) {
    static const float mPtFactor = 0.3*3.8*1e-2/2.0;
    //return simPhi - mPtFactor * r * simChargeOverPt;
    return simPhi - std::asin(mPtFactor * r * simChargeOverPt);
}

float calcIdealZ(float simVz, float simCotTheta, float simChargeOverPt, float r) {
    static const float mPtFactor = 0.3*3.8*1e-2/2.0;
    //return simVz + r * simCotTheta;
    return simVz + (1.0 / (mPtFactor * simChargeOverPt) * std::asin(mPtFactor * r * simChargeOverPt)) * simCotTheta;
}
}


// _____________________________________________________________________________
int StubCleanerPU::cleanStubs(TString src, TString out) {
    if (verbose_)  std::cout << Info() << "Reading " << nEvents_ << " events and cleaning them." << std::endl;

    // _________________________________________________________________________
    // For reading
    TTStubReader reader(verbose_);
    if (reader.init(src, false)) {
        std::cout << Error() << "Failed to initialize TTStubReader." << std::endl;
        return 1;
    }

    // For writing
    TTStubWriter writer(verbose_);
    if (writer.init(reader.getChain(), out)) {
        std::cout << Error() << "Failed to initialize TTStubWriter." << std::endl;
        return 1;
    }


    // _________________________________________________________________________
    // Loop over all events

    // Bookkeepers
    long int nRead2 = 0, nKept2 = 0;  // count each event
    long int nRead1 = 0, nKept1 = 0;  // count each particle

    for (long long ievt=0; ievt<nEvents_; ++ievt) {
        if (reader.loadTree(ievt) < 0)  break;
        reader.getEntry(ievt);

        // Make a read-only cache of the reader
        const TTStubReaderCache readerCache(reader);

        const unsigned nstubs = readerCache.vb_modId->size();
        const unsigned nparts = readerCache.vp_pt->size();
        if (verbose_>1 && ievt%50000==0)  std::cout << Debug() << Form("... Processing event   : %7lld, keeping: %7ld (%7ld particles)", ievt, nKept2, nKept1) << std::endl;
        if (verbose_>2)  std::cout << Debug() << "... evt: " << ievt << " # stubs: " << nstubs << " # particles: " << nparts << std::endl;

        bool keepevent = false;

        // Loop over genParticles
        for (unsigned ipart=0; ipart<nparts; ++ipart) {

            // _________________________________________________________________
            // Check sim info

            float simPt           = readerCache.vp_pt->at(ipart);
            float simEta          = readerCache.vp_eta->at(ipart);
            float simPhi          = readerCache.vp_phi->at(ipart);
            //float simVx           = readerCache.vp_vx->at(ipart);
            //float simVy           = readerCache.vp_vy->at(ipart);
            float simVz           = readerCache.vp_vz->at(ipart);
            int   simCharge       = readerCache.vp_charge->at(ipart);
            int   simPdgId        = readerCache.vp_pdgId->at(ipart);
            int   simStatus       = readerCache.vp_status->at(ipart);

            float simCotTheta     = std::sinh(simEta);
            float simChargeOverPt = float(simCharge)/simPt;

            // Apply charge and status requirements
            bool sim = (simCharge != 0 && simStatus == 1);
            if (!sim)
                continue;

            // Apply pt, eta, phi requirements
            sim = (po_.minPt  <= simPt  && simPt  <= po_.maxPt  &&
                   po_.minEta <= simEta && simEta <= po_.maxEta &&
                   po_.minPhi <= simPhi && simPhi <= po_.maxPhi &&
                   po_.minVz  <= simVz  && simVz  <= po_.maxVz);
            if (!sim)
                continue;

            if (verbose_>2)  std::cout << Debug() << "... ... part: " << ipart << " simPt: " << simPt << " simEta: " << simEta << " simPhi: " << simPhi << " simVz: " << simVz << " simChargeOverPt: " << simChargeOverPt << " simPdgId: " << simPdgId << " simStatus: " << simStatus << std::endl;

            // _________________________________________________________________
            // Start cleaning

            // Events that fail don't exit the loop immediately, so that event info
            // can still be printed when verbosity is turned on.
            bool keep = true;

            // Check min # of stubs
            bool require = (nstubs >= MIN_NGOODSTUBS);
            if (!require)
                keep = false;

            // _________________________________________________________________
            // Remove multiple stubs in one layer

            // Make a vector of pairs, each pair has an id and a 2D (R,D) value,
            // where R is rank based on radius or z coord, D is (dx**2 + dy**2 + dz**2)**(1/2)
            std::vector<std::pair<unsigned, std::pair<unsigned, float> > > vec_index_dist;

            for (unsigned istub=0; (istub<nstubs) && keep; ++istub) {
                //int tpId = readerCache.vb_tpId->at(istub);
                //if (tpId < 0)
                //    continue;

                unsigned moduleId = readerCache.vb_modId   ->at(istub);
                float    stub_r   = readerCache.vb_r       ->at(istub);
                float    stub_phi = readerCache.vb_phi     ->at(istub);
                float    stub_z   = readerCache.vb_z       ->at(istub);
                float    stub_ds  = readerCache.vb_trigBend->at(istub);

                unsigned lay16    = compressLayer(decodeLayer(moduleId));
                assert(lay16 < 16);

                // CUIDADO: simVx and simVy are currently not used in the calculation
                //          therefore d0 is assumed to be zero, and z0 is assumed to be equal to vz
                float idealPhi = calcIdealPhi(simPhi, simChargeOverPt, stub_r);
                float idealZ   = calcIdealZ(simVz, simCotTheta, simChargeOverPt, stub_r);
                float idealR   = stub_r;

                if (lay16 >= 6) {  // for endcap
                    idealR     = (stub_z - simVz) / simCotTheta;
                    //if (idealR <= 0) {
                    //    std::cout << Warning() << "Stub ideal r <= 0! moduleId: " << moduleId << " r: " << stub_r << " z: " << stub_z << " simVz: " << simVz << " simCotTheta: " << simCotTheta << std::endl;
                    //}
                    idealPhi   = calcIdealPhi(simPhi, simChargeOverPt, idealR);
                    idealZ     = stub_z;
                }

                float deltaPhi = stub_phi - idealPhi;
                float deltaZ   = stub_z - idealZ;
                float deltaR   = stub_r - idealR;

                if (verbose_>2)  std::cout << Debug() << "... ... ... stub: " << istub << " moduleId: " << moduleId << " r: " << stub_r << " phi: " << stub_phi << " z: " << stub_z << " ds: " << stub_ds << " lay16: " << lay16 << " deltaPhi: " << deltaPhi << " deltaR: " << deltaR << " deltaZ: " << deltaZ << std::endl;

                bool picked = picky_ -> applyCuts(lay16, deltaPhi, deltaR, deltaZ);
                if (!po_.picky || (po_.picky && picked) ) {
                    unsigned rank = picky_ -> findRank(lay16, stub_r, stub_z);

                    float deltaX = stub_r * (std::cos(stub_phi) - std::cos(idealPhi));
                    float deltaY = stub_r * (std::sin(stub_phi) - std::sin(idealPhi));
                    float dist   = std::sqrt(deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ);

                    if (lay16 >= 6) {  // for endcap
                        deltaX = stub_r * std::cos(stub_phi) - idealR * std::cos(idealPhi);
                        deltaY = stub_r * std::sin(stub_phi) - idealR * std::sin(idealPhi);
                        dist   = std::sqrt(deltaX*deltaX + deltaY*deltaY);
                    }

                    if (verbose_>2)  std::cout << Debug() << "... ... ... stub: " << istub << " rank: " << rank << " dist: " << dist << std::endl;

                    vec_index_dist.push_back(std::make_pair(istub, std::make_pair(rank, dist)));

                } else {
                    if (verbose_>2)  std::cout << Debug() << "... ... ... stub: " << istub << " fail cut!" << std::endl;
                }

            }  // end loop over stubs

            // Sort by rank, then by smallest dist to largest
            // For future: also include delta_s?
            std::sort(vec_index_dist.begin(), vec_index_dist.end(), sortByUnsignedThenFloat);

            // Select only one stub per layer
            std::vector<unsigned> goodIndices(16, 999999);
            if (vec_index_dist.size()) {
                for (unsigned iistub=0; iistub<vec_index_dist.size(); ++iistub) {
                    unsigned istub = vec_index_dist.at(iistub).first;
                    float    dist  = vec_index_dist.at(iistub).second.second;

                    unsigned moduleId = readerCache.vb_modId->at(istub);
                    unsigned lay16    = compressLayer(decodeLayer(moduleId));

                    // For each layer, takes the stub with min dist to simTrack
                    if (goodIndices.at(lay16) == 999999 && dist < 26.0) {  // gets rid of stubs due to loopers
                        goodIndices.at(lay16) = istub;
                    }
                }
            }

            //if (keep && goodIndices.at(0) == 999999)
            //    std::cout << Warning() << "... ... part: " << ipart << " no stub in the first layer of the barrel!" << std::endl;
            //if (keep && goodIndices.at(6) != 999999 && goodIndices.at(11) != 999999)
            //    std::cout << Warning() << "... ... part: " << ipart << " found stubs in the first layers of both positive and negative endcaps!" << std::endl;

            // _________________________________________________________________
            // Now make keep-or-ignore decision per stub
            unsigned ngoodstubs = 0;

            for (unsigned istub=0; (istub<nstubs) && keep; ++istub) {
                bool keepstub = true;

                unsigned moduleId = readerCache.vb_modId->at(istub);

                // Check whether istub was an index stored for a good stub
                const unsigned count = std::count(goodIndices.begin(), goodIndices.end(), istub);
                if (!count)
                    keepstub = false;

                if (keepstub) {
                    if (verbose_>2)  std::cout << Debug() << "... ... ... stub: " << istub << " moduleId: " << moduleId << " is selected!" << std::endl;

                    // Keep the stub and do something similar to insertion sort
                    // First, find the position to insert (determined by moduleId)
                    std::vector<unsigned>::const_iterator pos = std::upper_bound(readerCache.vb_modId->begin(), readerCache.vb_modId->begin()+ngoodstubs, moduleId);
                    unsigned ipos = pos - readerCache.vb_modId->begin();

                    // Insert while keeping only the 'ngoodstubs' elements
                    // (read from readerCache, write into reader)
                  //insertSorted(reader.vb_x->begin()         , ngoodstubs, ipos, readerCache.vb_x->at(istub));
                  //insertSorted(reader.vb_y->begin()         , ngoodstubs, ipos, readerCache.vb_y->at(istub));
                    insertSorted(reader.vb_z->begin()         , ngoodstubs, ipos, readerCache.vb_z->at(istub));
                    insertSorted(reader.vb_r->begin()         , ngoodstubs, ipos, readerCache.vb_r->at(istub));
                    insertSorted(reader.vb_eta->begin()       , ngoodstubs, ipos, readerCache.vb_eta->at(istub));
                    insertSorted(reader.vb_phi->begin()       , ngoodstubs, ipos, readerCache.vb_phi->at(istub));
                    insertSorted(reader.vb_coordx->begin()    , ngoodstubs, ipos, readerCache.vb_coordx->at(istub));
                    insertSorted(reader.vb_coordy->begin()    , ngoodstubs, ipos, readerCache.vb_coordy->at(istub));
                    insertSorted(reader.vb_trigBend->begin()  , ngoodstubs, ipos, readerCache.vb_trigBend->at(istub));
                  //insertSorted(reader.vb_roughPt->begin()   , ngoodstubs, ipos, readerCache.vb_roughPt->at(istub));
                  //insertSorted(reader.vb_clusWidth0->begin(), ngoodstubs, ipos, readerCache.vb_clusWidth0->at(istub));
                  //insertSorted(reader.vb_clusWidth1->begin(), ngoodstubs, ipos, readerCache.vb_clusWidth1->at(istub));
                    insertSorted(reader.vb_modId->begin()     , ngoodstubs, ipos, readerCache.vb_modId->at(istub));
                    insertSorted(reader.vb_tpId->begin()      , ngoodstubs, ipos, readerCache.vb_tpId->at(istub));

                    ++ngoodstubs;  // remember to increment
                }

            }  // end loop over stubs
            assert(ngoodstubs <= nstubs);

            unsigned ngoodparts = 0;
            {
                insertSorted(reader.vp_pt->begin()        , ngoodparts, 0u, readerCache.vp_pt->at(ipart));
                insertSorted(reader.vp_eta->begin()       , ngoodparts, 0u, readerCache.vp_eta->at(ipart));
                insertSorted(reader.vp_phi->begin()       , ngoodparts, 0u, readerCache.vp_phi->at(ipart));
                insertSorted(reader.vp_vx->begin()        , ngoodparts, 0u, readerCache.vp_vx->at(ipart));
                insertSorted(reader.vp_vy->begin()        , ngoodparts, 0u, readerCache.vp_vy->at(ipart));
                insertSorted(reader.vp_vz->begin()        , ngoodparts, 0u, readerCache.vp_vz->at(ipart));
                insertSorted(reader.vp_charge->begin()    , ngoodparts, 0u, readerCache.vp_charge->at(ipart));
                insertSorted(reader.vp_pdgId->begin()     , ngoodparts, 0u, readerCache.vp_pdgId->at(ipart));
                insertSorted(reader.vp_status->begin()    , ngoodparts, 0u, readerCache.vp_status->at(ipart));

                ++ngoodparts;
            }


            // _________________________________________________________________
            // Now make keep-or-ignore decision per particle

            // Check again min # of stubs
            require = (ngoodstubs >= MIN_NGOODSTUBS);
            if (!require)
                keep = false;

            if (keep) {
                ++nKept1;
                keepevent = true;
            } else {  // do not keep any stub
                ngoodstubs = 0;
            }

            if (keep && ngoodstubs > MAX_NGOODSTUBS) {
                std::cout << Warning() << "... ... part: " << ipart << " simPt: " << simPt << " simEta: " << simEta << " simPhi: " << simPhi <<  " ngoodstubs: " << ngoodstubs << std::endl;
            }

            if (verbose_>2)  std::cout << Debug() << "... ... part: " << ipart << " # good stubs: " << ngoodstubs << " keep? " << keep << std::endl;

          //reader.vb_x         ->resize(ngoodstubs);
          //reader.vb_y         ->resize(ngoodstubs);
            reader.vb_z         ->resize(ngoodstubs);
            reader.vb_r         ->resize(ngoodstubs);
            reader.vb_eta       ->resize(ngoodstubs);
            reader.vb_phi       ->resize(ngoodstubs);
            reader.vb_coordx    ->resize(ngoodstubs);
            reader.vb_coordy    ->resize(ngoodstubs);
            reader.vb_trigBend  ->resize(ngoodstubs);
          //reader.vb_roughPt   ->resize(ngoodstubs);
          //reader.vb_clusWidth0->resize(ngoodstubs);
          //reader.vb_clusWidth1->resize(ngoodstubs);
            reader.vb_modId     ->resize(ngoodstubs);
            reader.vb_tpId      ->resize(ngoodstubs);

            reader.vp_pt        ->resize(ngoodparts);
            reader.vp_eta       ->resize(ngoodparts);
            reader.vp_phi       ->resize(ngoodparts);
            reader.vp_vx        ->resize(ngoodparts);
            reader.vp_vy        ->resize(ngoodparts);
            reader.vp_vz        ->resize(ngoodparts);
            reader.vp_charge    ->resize(ngoodparts);
            reader.vp_pdgId     ->resize(ngoodparts);
            reader.vp_status    ->resize(ngoodparts);

            ++nRead1;
            writer.fill();

            // Resize to the original size
          //reader.vb_x         ->resize(nstubs);
          //reader.vb_y         ->resize(nstubs);
            reader.vb_z         ->resize(nstubs);
            reader.vb_r         ->resize(nstubs);
            reader.vb_eta       ->resize(nstubs);
            reader.vb_phi       ->resize(nstubs);
            reader.vb_coordx    ->resize(nstubs);
            reader.vb_coordy    ->resize(nstubs);
            reader.vb_trigBend  ->resize(nstubs);
          //reader.vb_roughPt   ->resize(nstubs);
          //reader.vb_clusWidth0->resize(nstubs);
          //reader.vb_clusWidth1->resize(nstubs);
            reader.vb_modId     ->resize(nstubs);
            reader.vb_tpId      ->resize(nstubs);

            reader.vp_pt        ->resize(nparts);
            reader.vp_eta       ->resize(nparts);
            reader.vp_phi       ->resize(nparts);
            reader.vp_vx        ->resize(nparts);
            reader.vp_vy        ->resize(nparts);
            reader.vp_vz        ->resize(nparts);
            reader.vp_charge    ->resize(nparts);
            reader.vp_pdgId     ->resize(nparts);
            reader.vp_status    ->resize(nparts);
        }  // end loop over genParticles

        if (keepevent)
            ++nKept2;

        ++nRead2;
    }  // end loop over events

    if (nRead2 == 0) {
        std::cout << Error() << "Failed to read any event." << std::endl;
        return 1;
    }

    if (verbose_)  std::cout << Info() << Form("Read event   : %7ld, kept: %7ld", nRead2, nKept2) << std::endl;
    if (verbose_)  std::cout << Info() << Form("Read particle: %7ld, kept: %7ld", nRead1, nKept1) << std::endl;

    long long nentries = writer.writeTree();
    assert(nentries == nRead1);

    return 0;
}


// _____________________________________________________________________________
// Main driver
int StubCleanerPU::run() {
    int exitcode = 0;
    Timing(1);

    exitcode = cleanStubs(po_.input, po_.output);
    if (exitcode)  return exitcode;
    Timing();

    return exitcode;
}
