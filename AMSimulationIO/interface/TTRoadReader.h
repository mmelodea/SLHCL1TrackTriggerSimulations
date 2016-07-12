#ifndef AMSimulationIO_TTRoadReader_h_
#define AMSimulationIO_TTRoadReader_h_

#include "SLHCL1TrackTriggerSimulations/AMSimulationDataFormats/interface/TTRoad.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTStubPlusTPReader.h"


namespace slhcl1tt {

// _____________________________________________________________________________
class TTRoadReader : public TTStubPlusTPReader {
  public:
    TTRoadReader(int verbose=1);
    ~TTRoadReader();

    int init(TString src, TString prefix, TString suffix);

    // Stubs
    // Adding stubs branches here is inconsistent and ugly, but changing
    // TTStubReader is problematic now and breaks backward compatibility
    std::vector<std::string> *                          vb_bitString;
    std::vector<unsigned> *                             vb_superstripId;

    // Roads
    std::vector<unsigned> *                             vr_patternRef;
    std::vector<unsigned> *                             vr_tower;
    std::vector<unsigned> *                             vr_nstubs;
    std::vector<float> *                                vr_patternInvPt;
    std::vector<std::vector<unsigned> > *               vr_superstripIds;
    std::vector<std::vector<std::vector<unsigned> > > * vr_stubRefs;
};


// _____________________________________________________________________________
class TTRoadWriter : public BasicWriter {
  public:
    TTRoadWriter(int verbose=1);
    ~TTRoadWriter();

    int init(TChain* tchain, TString out, TString prefix, TString suffix);

    void fill(const std::vector<TTRoad>& roads);

    void fill(const std::vector<TTRoad>& roads, const std::vector<std::string>& stubs_bitString, const std::vector<unsigned>& stubs_superstripId);

  protected:
    // Stubs
    // Adding stubs branches here is inconsistent and ugly, but changing
    // TTStubReader is problematic now and breaks backward compatibility
    std::auto_ptr<std::vector<std::string> >                          vb_bitString;
    std::auto_ptr<std::vector<unsigned> >                             vb_superstripId;

    // Roads
    std::auto_ptr<std::vector<unsigned> >                             vr_patternRef;
    std::auto_ptr<std::vector<unsigned> >                             vr_tower;
    std::auto_ptr<std::vector<unsigned> >                             vr_nstubs;
    std::auto_ptr<std::vector<float> >                                vr_patternInvPt;
    std::auto_ptr<std::vector<std::vector<unsigned> > >               vr_superstripIds;
    std::auto_ptr<std::vector<std::vector<std::vector<unsigned> > > > vr_stubRefs;
};

}  // namespace slhcl1tt

#endif
