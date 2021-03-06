#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/PatternAnalyzer.h"

#include <fstream>

namespace {
// Comparator
bool sortByInvPt(const std::pair<pattern_type, Attributes *>& lhs, const std::pair<pattern_type, Attributes *>& rhs) {
    return (lhs.second)->invPt.getMean() < (rhs.second)->invPt.getMean();
}
}


// _____________________________________________________________________________
int PatternAnalyzer::loadPatterns(TString bank) {
    if (verbose_)  std::cout << Info() << "Loading patterns from " << bank << std::endl;

    // _________________________________________________________________________
    // For reading pattern bank
    PatternBankReader pbreader(verbose_);
    pbreader.init(bank);

    long long npatterns = pbreader.getEntries();
    if (npatterns > po_.maxPatterns)
        npatterns = po_.maxPatterns;
    assert(npatterns > 0);

    // Allocate memory
    patternAttributes_.clear();
    patternAttributes_.resize(npatterns);

    // _________________________________________________________________________
    // Load the patterns

    pattern_type patt;
    patt.fill(0);
    std::pair<std::map<pattern_type, Attributes *>::iterator,bool> ret;

    for (long long ipatt=0; ipatt<npatterns; ++ipatt) {
        pbreader.getPattern(ipatt);
        if (pbreader.pb_frequency < po_.minFrequency)
            break;

        assert(pbreader.pb_superstripIds->size() == po_.nLayers);

        // Make the pattern
        std::copy(pbreader.pb_superstripIds->begin(), pbreader.pb_superstripIds->end(), patt.begin());

        // Insert the pattern
        ret = patternAttributes_map_.insert(std::make_pair(patt, &(patternAttributes_.at(ipatt))));

        if (!ret.second) {
            std::cout << Warning() << "Failed to insert: " << patt << std::endl;
        }
    }
    if (verbose_)  std::cout << Info() << "Successfully loaded " << npatterns << " patterns." << std::endl;

    return 0;
}

// _____________________________________________________________________________
// Make the patterns
int PatternAnalyzer::makePatterns(TString src) {
    if (verbose_)  std::cout << Info() << "Reading " << nEvents_ << " events and generating patterns." << std::endl;

    // _________________________________________________________________________
    // For reading
    TTStubReader reader(verbose_);
    reader.init(src);

    // _________________________________________________________________________
    // Get trigger tower reverse map
    const std::map<unsigned, bool>& ttrmap = ttmap_ -> getTriggerTowerReverseMap(po_.tower);

    // _________________________________________________________________________
    // Book histograms
    TH1::AddDirectory(kFALSE);
    histogram2Ds["diff_invPt_vs_invPt1"] = new TH2F("diff_invPt_vs_invPt1", "; signed 1/p_{T} [1/GeV]; Diff from mean value of road [1/GeV]", 1000, -0.5, 0.5, 1000, -0.05, 0.05);
    histogram2Ds["diff_invPt_vs_invPt2"] = new TH2F("diff_invPt_vs_invPt2", "; signed 1/p_{T} [1/GeV]; Diff from linear approx [1/GeV]", 1000, -0.5, 0.5, 1000, -0.05, 0.05);
    histogramPRs["diff_invPt_pr_invPt1"] = new TProfile("diff_invPt_pr_invPt1", "; signed 1/p_{T} [1/GeV]; Diff from mean value of road [1/GeV]", 1000, -0.5, 0.5, "s");
    histogramPRs["diff_invPt_pr_invPt2"] = new TProfile("diff_invPt_pr_invPt2", "; signed 1/p_{T} [1/GeV]; Diff from linear approx [1/GeV]", 1000, -0.5, 0.5, "s");


    // _________________________________________________________________________
    // Loop over all events (filter)

    if (verbose_)  std::cout << Info() << "Begin event filtering" << std::endl;

    // Event decisions
    std::vector<bool> keepEvents;

    // Bookkeepers
    long int nRead = 0, nKept = 0;

    for (long long ievt=0; ievt<nEvents_; ++ievt) {
        if (reader.loadTree(ievt) < 0)  break;
        reader.getEntry(ievt);

        const unsigned nstubs = reader.vb_modId->size();
        if (verbose_>1 && ievt%100000==0)  std::cout << Debug() << Form("... Processing event: %7lld, keeping: %7ld", ievt, nKept) << std::endl;
        if (verbose_>2)  std::cout << Debug() << "... evt: " << ievt << " # stubs: " << nstubs << std::endl;

        // Apply track pt requirement
        float simPt = reader.vp_pt->front();
        if (simPt < po_.minPt || po_.maxPt < simPt) {
            ++nRead;
            keepEvents.push_back(false);
            continue;
        }

        // Apply trigger tower acceptance
        unsigned ngoodstubs = 0;
        for (unsigned istub=0; istub<nstubs; ++istub) {
            unsigned moduleId = reader.vb_modId   ->at(istub);
            if (ttrmap.find(moduleId) != ttrmap.end()) {
                ++ngoodstubs;
            }
        }
        if (ngoodstubs != po_.nLayers) {
            ++nRead;
            keepEvents.push_back(false);
            continue;
        }
        assert(nstubs == po_.nLayers);

        ++nKept;
        ++nRead;
        keepEvents.push_back(true);
    }

    if (nRead == 0) {
        std::cout << Error() << "Failed to read any event." << std::endl;
        return 1;
    }

    if (verbose_)  std::cout << Info() << Form("Read: %7ld, kept: %7ld", nRead, nKept) << std::endl;


    // _________________________________________________________________________
    // Loop over all events

    if (verbose_)  std::cout << Info() << "Begin first loop on tracks" << std::endl;

    pattern_type patt;
    patt.fill(0);

    // Save pointer to the attribute for every valid track
    std::vector<std::pair<float, Attributes *> >  attrs;

    // Bookkeepers
    nRead = 0, nKept = 0;


    //Get the layers to be applied the deltaS ss definition                                                                                                 
    std::vector<unsigned int> dslayers;
    int deltaS_layers = po_.deltaS;
    while( deltaS_layers > 0 ){
      int ilayer = deltaS_layers % 10;
      dslayers.push_back(ilayer);
      deltaS_layers /= 10;
    }
    std::reverse(dslayers.begin(),dslayers.end());
    const int ndslayers = dslayers.size();
    std::cout<<"Configuration for SS ID: "<< po_.deltaSM <<" -- (layer/#bins): ";
    for(int idsl=0; idsl<ndslayers; ++idsl)
      std::cout<<" l"<<idsl<<"/bs"<<dslayers[idsl];
    std::cout<<std::endl;


    for (long long ievt=0; ievt<nEvents_; ++ievt) {
        if (reader.loadTree(ievt) < 0)  break;
        reader.getEntry(ievt);

        const unsigned nstubs = reader.vb_modId->size();
        if (verbose_>1 && ievt%100000==0)  std::cout << Debug() << Form("... Processing event: %7lld, keeping: %7ld", ievt, nKept) << std::endl;

        if (!keepEvents.at(ievt)) {
            ++nRead;
            continue;
        }

        // Get sim info
        assert(reader.vp_pt->size() == 1);
        float simPt           = reader.vp_pt->front();
        float simEta          = reader.vp_eta->front();
        float simPhi          = reader.vp_phi->front();
        //float simVx           = reader.vp_vx->front();
        //float simVy           = reader.vp_vy->front();
        float simVz           = reader.vp_vz->front();
        int   simCharge       = reader.vp_charge->front();

        float simCotTheta     = std::sinh(simEta);
        float simChargeOverPt = float(simCharge)/simPt;

        // _____________________________________________________________________
        // Start generating patterns
        patt.fill(0);

        // Loop over reconstructed stubs
        for (unsigned istub=0; istub<nstubs; ++istub) {
            unsigned moduleId = reader.vb_modId   ->at(istub);
            float    strip    = reader.vb_coordx  ->at(istub);  // in full-strip unit
            float    segment  = reader.vb_coordy  ->at(istub);  // in full-strip unit

            float    stub_r   = reader.vb_r       ->at(istub);
            float    stub_phi = reader.vb_phi     ->at(istub);
            float    stub_z   = reader.vb_z       ->at(istub);
            float    stub_ds  = reader.vb_trigBend->at(istub);  // in full-strip unit

            //Get layer & number of bins to split deltaS
            unsigned lay16 = compressLayer(decodeLayer(moduleId));
            unsigned nDSbins = dslayers[lay16];

            // Find superstrip ID
            unsigned ssId = 0;
            if (!arbiter_ -> useGlobalCoord()) {  // local coordinates
                ssId = arbiter_ -> superstripLocal(moduleId, strip, segment);

            } else {                              // global coordinates
	      ssId = arbiter_ -> superstripGlobal(moduleId, stub_r, stub_phi, stub_z, stub_ds, nDSbins, po_.deltaSM);
            }
            patt.at(istub) = ssId;

            if (verbose_>2) {
                std::cout << Debug() << "... ... stub: " << istub << " moduleId: " << moduleId << " strip: " << strip << " segment: " << segment << " r: " << stub_r << " phi: " << stub_phi << " z: " << stub_z << " ds: " << stub_ds << std::endl;
                std::cout << Debug() << "... ... stub: " << istub << " ssId: " << ssId << std::endl;
            }

        }

        // Find pattern in the bank, update the attributes
        std::map<pattern_type, Attributes *>::iterator found = patternAttributes_map_.find(patt);
        if (found != patternAttributes_map_.end()) {
            Attributes * attr = found->second;
            ++ attr->n;
            attr->invPt.fill(simChargeOverPt);
            attr->cotTheta.fill(simCotTheta);
            attr->phi.fill(simPhi);
            attr->z0.fill(simVz);

            attrs.push_back(std::make_pair(simChargeOverPt, attr));

        } else {
            //std::cout << Warning() << "Failed to find: " << patt << std::endl;

            keepEvents.at(ievt) = false;
        }

        if (verbose_>2)  std::cout << Debug() << "... evt: " << ievt << " patt: " << patt << std::endl;

        ++nKept;
        ++nRead;
    }


    // _________________________________________________________________________
    // Sort by mean invPt

    // Convert map to vector of pairs
    const unsigned origSize = patternAttributes_map_.size();
    //patternAttributes_pairs_.reserve(patternAttributes_map_.size());  // can cause bad_alloc
    //patternAttributes_pairs_.insert(patternAttributes_pairs_.end(), patternAttributes_map_.begin(), patternAttributes_map_.end());

    for (std::map<pattern_type, Attributes *>::const_iterator it = patternAttributes_map_.begin();
         it != patternAttributes_map_.end(); ) {  // should not cause bad_alloc
        patternAttributes_pairs_.push_back(*it);
        it = patternAttributes_map_.erase(it);
    }
    assert(patternAttributes_pairs_.size() == origSize);

    // Clear map and release memory
    std::map<pattern_type, Attributes *> mapEmpty;
    patternAttributes_map_.clear();
    patternAttributes_map_.swap(mapEmpty);

    // Sort by invPt
    std::stable_sort(patternAttributes_pairs_.begin(), patternAttributes_pairs_.end(), sortByInvPt);

    // Assign sorted pattern id
    for (unsigned i=0; i<patternAttributes_pairs_.size(); ++i) {
        patternAttributes_pairs_.at(i).second->id = i;
    }


    // _________________________________________________________________________
    // Fill histograms

    // CUIDADO: this assumes all the pattern attributes are populated well (ntracks > npatterns)
    const unsigned npatterns = patternAttributes_pairs_.size();
    for (std::vector<std::pair<float, Attributes *> >::const_iterator it=attrs.begin();
         it!=attrs.end(); ++it) {

        float approx = 1.0 / po_.minPt * (-1. + 2./(npatterns-1) * it->second->id);
        if (it->second->n <= 1)
            std::cout << Warning() << "Too few entries: " << it->second->n << " id: " << it->second->id << std::endl;
        if (verbose_>3)
            std::cout << Debug() << "... good evt: " << it-attrs.begin() << " invPt: " << it->first << " mean: " << it->second->invPt.getMean() << " id: " << it->second->id << " approx: " << approx << std::endl;

        // Mean is like a measured quantity. Usually difference is defined as "measured" - "true"
        histogram2Ds["diff_invPt_vs_invPt1"]->Fill(it->first, it->second->invPt.getMean() - it->first);
        histogram2Ds["diff_invPt_vs_invPt2"]->Fill(it->first, approx - it->first);
        histogramPRs["diff_invPt_pr_invPt1"]->Fill(it->first, it->second->invPt.getMean() - it->first);
        histogramPRs["diff_invPt_pr_invPt2"]->Fill(it->first, approx - it->first);
    }


    // _________________________________________________________________________
    // Dump event decisions

    TString out = "keepEvents.txt";
    std::ofstream outfile(out.Data());
    if (!outfile) {
        std::cout << Error() << "Unable to open " << out << std::endl;
        return 1;
    }
    for (unsigned i=0; i<keepEvents.size(); ++i)
        outfile << keepEvents.at(i) << std::endl;
    outfile.close();

    return 0;
}


// _____________________________________________________________________________
// Output patterns into a TTree
int PatternAnalyzer::writePatterns(TString out) {

    // _________________________________________________________________________
    // For writing
    PatternBankWriter writer(verbose_);
    writer.init(out);

    // _________________________________________________________________________
    // Save pattern bank statistics
    *(writer.pb_coverage)   = 0.;  // dummy
    *(writer.pb_count)      = 0;   // dummy
    *(writer.pb_tower)      = po_.tower;
    *(writer.pb_superstrip) = po_.superstrip;
    writer.fillPatternBankInfo();

    // _________________________________________________________________________
    // Save pattern bank
    const long long npatterns = patternAttributes_pairs_.size();

    // Bookkeepers
    long int nKept = 0;

    for (long long ipatt=0; ipatt<npatterns; ++ipatt) {
        if (verbose_>1 && ipatt%1000==0) {
            std::cout << Debug() << Form("... Writing event: %7lld, total freq: %7lu", ipatt, nKept) << std::endl;
        }

        const Attributes * attr = patternAttributes_pairs_.at(ipatt).second;

        writer.pb_superstripIds->clear();
        const pattern_type& patt = patternAttributes_pairs_.at(ipatt).first;
        for (unsigned ilayer=0; ilayer<po_.nLayers; ++ilayer) {
            writer.pb_superstripIds->push_back(patt.at(ilayer));
        }
        *(writer.pb_frequency)      = attr->n;

        *(writer.pb_invPt_mean)     = attr->invPt.getMean();
        *(writer.pb_invPt_sigma)    = attr->invPt.getSigma();
        *(writer.pb_cotTheta_mean)  = attr->cotTheta.getMean();
        *(writer.pb_cotTheta_sigma) = attr->cotTheta.getSigma();
        *(writer.pb_phi_mean)       = attr->phi.getMean();
        *(writer.pb_phi_sigma)      = attr->phi.getSigma();
        *(writer.pb_z0_mean)        = attr->z0.getMean();
        *(writer.pb_z0_sigma)       = attr->z0.getSigma();

        writer.fillPatternBank();
        writer.fillPatternAttributes();

        nKept += attr->n;
    }

    // _________________________________________________________________________
    // Write histograms
    for (std::map<TString, TH1F *>::const_iterator it=histograms.begin();
         it!=histograms.end(); ++it) {
        if (it->second)  it->second->SetDirectory(gDirectory);
    }

    for (std::map<TString, TH2F *>::const_iterator it=histogram2Ds.begin();
         it!=histogram2Ds.end(); ++it) {
        if (it->second)  it->second->SetDirectory(gDirectory);
    }

    for (std::map<TString, TProfile *>::const_iterator it=histogramPRs.begin();
         it!=histogramPRs.end(); ++it) {
        if (it->second)  it->second->SetDirectory(gDirectory);
    }

    writer.write();

    return 0;
}


// _____________________________________________________________________________
// Main driver
int PatternAnalyzer::run() {
    int exitcode = 0;
    Timing(1);

    exitcode = loadPatterns(po_.bankfile);
    if (exitcode)  return exitcode;
    Timing();

    exitcode = makePatterns(po_.input);
    if (exitcode)  return exitcode;
    Timing();

    exitcode = writePatterns(po_.output);
    if (exitcode)  return exitcode;
    Timing();

    return exitcode;
}
