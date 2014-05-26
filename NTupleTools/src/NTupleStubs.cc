#include "SLHCL1TrackTriggerSimulations/NTupleTools/interface/NTupleStubs.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/StackedTrackerDetId.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"
//#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
//#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"


unsigned NTupleStubs::findById(const std::vector<unsigned>& vec, unsigned id, bool throwError) {
    std::vector<unsigned>::const_iterator found;
    found = std::find(vec.begin(), vec.end(), id);
    if (found == vec.end() && throwError)
        throw cms::Exception("LogicError") << "Id not found: " << id << ".\n";
    return (found - vec.begin());
}

unsigned NTupleStubs::getModuleLayer(const DetId& id) {
    if (id.subdetId() == (int) PixelSubdetector::PixelBarrel) {
        PXBDetId pxbId(id);
        uint32_t layer  = pxbId.layer();
        return layer;
    } else if (id.subdetId() == (int) PixelSubdetector::PixelEndcap) {
        PXFDetId pxfId(id);
        uint32_t side   = pxfId.side();
        uint32_t disk   = pxfId.disk();
        uint32_t layer  = (side == 2) ? disk : disk+7;
        return layer;
    } else { // error
        return 999999;
    }
}

unsigned NTupleStubs::getModuleLadder(const DetId& id) {
    if (id.subdetId() == (int) PixelSubdetector::PixelBarrel) {
        PXBDetId pxbId(id);
        uint32_t ladder = pxbId.ladder();
        return ladder;
    } else if (id.subdetId() == (int) PixelSubdetector::PixelEndcap) {
        PXFDetId pxfId(id);
        uint32_t ring   = pxfId.ring();
        uint32_t blade  = pxfId.blade();
        uint32_t panel  = pxfId.panel();
        if(panel != 1 || ring != blade) {  // error: panel must be 1 and ring is the same as blade now
            return 999999;
        }
        return ring;
    } else {  // error
        return 999999;
    }
}

unsigned NTupleStubs::getModuleModule(const DetId& id) {
    if (id.subdetId() == (int) PixelSubdetector::PixelBarrel) {
        PXBDetId pxbId(id);
        uint32_t module = pxbId.module();
        return module;
    } else if (id.subdetId() == (int) PixelSubdetector::PixelEndcap) {
        PXFDetId pxfId(id);
        uint32_t module = pxfId.module();
        return module;
    } else { // error
        return 999999;
    }
}

unsigned NTupleStubs::getModuleId(const DetId& id) {
    uint32_t layer  = getModuleLayer(id);
    uint32_t ladder = getModuleLadder(id);
    uint32_t module = getModuleModule(id);
    return 10000*layer + 100*(ladder-1) + (module-1)/2;
}


NTupleStubs::NTupleStubs(const edm::ParameterSet& iConfig) :
  inputTagClus_        (iConfig.getParameter<edm::InputTag>("inputTagClus")),
  inputTagStub_        (iConfig.getParameter<edm::InputTag>("inputTagStub")),
  inputTagTrack_       (iConfig.getParameter<edm::InputTag>("inputTagTrack")),
  inputTagDigi_        (iConfig.getParameter<edm::InputTag>("inputTagDigi")),
  inputTagClusMCAssoc_ (iConfig.getParameter<edm::InputTag>("inputTagClusMCAssoc")),
  inputTagStubMCAssoc_ (iConfig.getParameter<edm::InputTag>("inputTagStubMCAssoc")),
  inputTagTrackMCAssoc_(iConfig.getParameter<edm::InputTag>("inputTagTrackMCAssoc")),
  prefixClus_          (iConfig.getParameter<std::string>("prefixClus")),
  prefixStub_          (iConfig.getParameter<std::string>("prefixStub")),
  prefixTrack_         (iConfig.getParameter<std::string>("prefixTrack")),
  prefixDigi_          (iConfig.getParameter<std::string>("prefixDigi")),
  suffix_              (iConfig.getParameter<std::string>("suffix")) {

    // Remember r is rho in cylindrical coordinate system
    produces<std::vector<float> >                   (prefixClus_ + "x"              + suffix_);
    produces<std::vector<float> >                   (prefixClus_ + "y"              + suffix_);
    produces<std::vector<float> >                   (prefixClus_ + "z"              + suffix_);
    produces<std::vector<float> >                   (prefixClus_ + "r"              + suffix_);
    produces<std::vector<float> >                   (prefixClus_ + "phi"            + suffix_);
    produces<std::vector<float> >                   (prefixClus_ + "localx"         + suffix_);
    produces<std::vector<float> >                   (prefixClus_ + "localy"         + suffix_);
    produces<std::vector<float> >                   (prefixClus_ + "coordx"         + suffix_);
    produces<std::vector<float> >                   (prefixClus_ + "coordy"         + suffix_);
    produces<std::vector<float> >                   (prefixClus_ + "surfx"          + suffix_);
    produces<std::vector<float> >                   (prefixClus_ + "surfy"          + suffix_);
    produces<std::vector<float> >                   (prefixClus_ + "surfz"          + suffix_);
    produces<std::vector<unsigned> >                (prefixClus_ + "stackId"        + suffix_);
    produces<std::vector<unsigned> >                (prefixClus_ + "iLayer"         + suffix_);
    produces<std::vector<unsigned> >                (prefixClus_ + "iRing"          + suffix_);
    produces<std::vector<unsigned> >                (prefixClus_ + "iSide"          + suffix_);
    produces<std::vector<unsigned> >                (prefixClus_ + "iPhi"           + suffix_);
    produces<std::vector<unsigned> >                (prefixClus_ + "iZ"             + suffix_);
    produces<std::vector<unsigned> >                (prefixClus_ + "iModLayer"      + suffix_);
    produces<std::vector<unsigned> >                (prefixClus_ + "iModLadder"     + suffix_);
    produces<std::vector<unsigned> >                (prefixClus_ + "iModModule"     + suffix_);
    produces<std::vector<unsigned> >                (prefixClus_ + "iModCols"       + suffix_);
    produces<std::vector<unsigned> >                (prefixClus_ + "iModRows"       + suffix_);
    produces<std::vector<float> >                   (prefixClus_ + "iModPitchX"     + suffix_);
    produces<std::vector<float> >                   (prefixClus_ + "iModPitchY"     + suffix_);
    produces<std::vector<bool> >                    (prefixClus_ + "barrel"         + suffix_);
    produces<std::vector<bool> >                    (prefixClus_ + "psmodule"       + suffix_);
    produces<std::vector<unsigned> >                (prefixClus_ + "modId"          + suffix_);
    produces<std::vector<unsigned> >                (prefixClus_ + "geoId"          + suffix_);
    produces<std::vector<unsigned> >                (prefixClus_ + "stack"          + suffix_);
    produces<std::vector<unsigned> >                (prefixClus_ + "width"          + suffix_);
    produces<std::vector<unsigned> >                (prefixClus_ + "nhits"          + suffix_);
    produces<std::vector<std::vector<int> > >       (prefixClus_ + "hitCols"        + suffix_);
    produces<std::vector<std::vector<int> > >       (prefixClus_ + "hitRows"        + suffix_);
    produces<std::vector<std::vector<int> > >       (prefixClus_ + "hitADCs"        + suffix_);
    produces<std::vector<std::vector<int> > >       (prefixClus_ + "hitChans"       + suffix_);
    produces<std::vector<std::vector<int> > >       (prefixClus_ + "hitTrkIds"      + suffix_);
    produces<std::vector<std::vector<unsigned> > >  (prefixClus_ + "hitEvtIds"      + suffix_);
    produces<std::vector<std::vector<float> > >     (prefixClus_ + "hitFracs"       + suffix_);
    produces<std::vector<std::vector<float> > >     (prefixClus_ + "hitXs"          + suffix_);
    produces<std::vector<std::vector<float> > >     (prefixClus_ + "hitYs"          + suffix_);
    produces<std::vector<std::vector<float> > >     (prefixClus_ + "hitZs"          + suffix_);
    produces<std::vector<bool> >                    (prefixClus_ + "isGenuine"      + suffix_);
    produces<std::vector<bool> >                    (prefixClus_ + "isUnknown"      + suffix_);
    produces<std::vector<bool> >                    (prefixClus_ + "isCombinatoric" + suffix_);
    produces<std::vector<int> >                     (prefixClus_ + "tpId"           + suffix_);
    produces<std::vector<int> >                     (prefixClus_ + "trkId"          + suffix_);
    produces<std::vector<int> >                     (prefixClus_ + "pdgId"          + suffix_);
    produces<std::vector<float> >                   (prefixClus_ + "simPt"          + suffix_);
    produces<std::vector<float> >                   (prefixClus_ + "simEta"         + suffix_);
    produces<std::vector<float> >                   (prefixClus_ + "simPhi"         + suffix_);
    produces<unsigned>                              (prefixClus_ + "size"           + suffix_);

    produces<std::vector<float> >                   (prefixStub_ + "x"              + suffix_);
    produces<std::vector<float> >                   (prefixStub_ + "y"              + suffix_);
    produces<std::vector<float> >                   (prefixStub_ + "z"              + suffix_);
    produces<std::vector<float> >                   (prefixStub_ + "r"              + suffix_);
    produces<std::vector<float> >                   (prefixStub_ + "phi"            + suffix_);
    produces<std::vector<float> >                   (prefixStub_ + "coordx"         + suffix_);
    produces<std::vector<float> >                   (prefixStub_ + "coordy"         + suffix_);
    produces<std::vector<float> >                   (prefixStub_ + "dirx"           + suffix_);
    produces<std::vector<float> >                   (prefixStub_ + "diry"           + suffix_);
    produces<std::vector<float> >                   (prefixStub_ + "dirz"           + suffix_);
    produces<std::vector<float> >                   (prefixStub_ + "roughPt"        + suffix_);
    produces<std::vector<unsigned> >                (prefixStub_ + "stackId"        + suffix_);
    produces<std::vector<unsigned> >                (prefixStub_ + "iLayer"         + suffix_);
    produces<std::vector<unsigned> >                (prefixStub_ + "iRing"          + suffix_);
    produces<std::vector<unsigned> >                (prefixStub_ + "iSide"          + suffix_);
    produces<std::vector<unsigned> >                (prefixStub_ + "iPhi"           + suffix_);
    produces<std::vector<unsigned> >                (prefixStub_ + "iZ"             + suffix_);
    produces<std::vector<unsigned> >                (prefixStub_ + "iModLayer"      + suffix_);
    produces<std::vector<unsigned> >                (prefixStub_ + "iModLadder"     + suffix_);
    produces<std::vector<unsigned> >                (prefixStub_ + "iModModule"     + suffix_);
    produces<std::vector<unsigned> >                (prefixStub_ + "iModCols"       + suffix_);
    produces<std::vector<unsigned> >                (prefixStub_ + "iModRows"       + suffix_);
    produces<std::vector<float> >                   (prefixStub_ + "iModPitchX"     + suffix_);
    produces<std::vector<float> >                   (prefixStub_ + "iModPitchY"     + suffix_);
    produces<std::vector<bool> >                    (prefixStub_ + "barrel"         + suffix_);
    produces<std::vector<bool> >                    (prefixStub_ + "psmodule"       + suffix_);
    produces<std::vector<unsigned> >                (prefixStub_ + "modId"          + suffix_);
    produces<std::vector<unsigned> >                (prefixStub_ + "geoId0"         + suffix_);
    produces<std::vector<unsigned> >                (prefixStub_ + "geoId1"         + suffix_);
    produces<std::vector<unsigned> >                (prefixStub_ + "clusId0"        + suffix_);
    produces<std::vector<unsigned> >                (prefixStub_ + "clusId1"        + suffix_);
    produces<std::vector<unsigned> >                (prefixStub_ + "nhits"          + suffix_);
    produces<std::vector<std::vector<int> > >       (prefixStub_ + "hitCols"        + suffix_);
    produces<std::vector<std::vector<int> > >       (prefixStub_ + "hitRows"        + suffix_);
    produces<std::vector<std::vector<int> > >       (prefixStub_ + "hitADCs"        + suffix_);
    produces<std::vector<std::vector<int> > >       (prefixStub_ + "hitChans"       + suffix_);
    produces<std::vector<std::vector<int> > >       (prefixStub_ + "hitTrkIds"      + suffix_);
    produces<std::vector<std::vector<unsigned> > >  (prefixStub_ + "hitEvtIds"      + suffix_);
    produces<std::vector<std::vector<float> > >     (prefixStub_ + "hitFracs"       + suffix_);
    produces<std::vector<std::vector<float> > >     (prefixStub_ + "hitXs"          + suffix_);
    produces<std::vector<std::vector<float> > >     (prefixStub_ + "hitYs"          + suffix_);
    produces<std::vector<std::vector<float> > >     (prefixStub_ + "hitZs"          + suffix_);
    produces<std::vector<float> >                   (prefixStub_ + "trigDist"       + suffix_);
    produces<std::vector<float> >                   (prefixStub_ + "trigOffset"     + suffix_);
    produces<std::vector<float> >                   (prefixStub_ + "trigPos"        + suffix_);
    produces<std::vector<float> >                   (prefixStub_ + "trigBend"       + suffix_);
    produces<std::vector<float> >                   (prefixStub_ + "separation"     + suffix_);
    produces<std::vector<int> >                     (prefixStub_ + "detWindow"      + suffix_);
    produces<std::vector<bool> >                    (prefixStub_ + "isGenuine"      + suffix_);
    produces<std::vector<bool> >                    (prefixStub_ + "isUnknown"      + suffix_);
    produces<std::vector<bool> >                    (prefixStub_ + "isCombinatoric" + suffix_);
    produces<std::vector<int> >                     (prefixStub_ + "tpId"           + suffix_);
    produces<std::vector<int> >                     (prefixStub_ + "trkId"          + suffix_);
    produces<std::vector<int> >                     (prefixStub_ + "pdgId"          + suffix_);
    produces<std::vector<float> >                   (prefixStub_ + "simPt"          + suffix_);
    produces<std::vector<float> >                   (prefixStub_ + "simEta"         + suffix_);
    produces<std::vector<float> >                   (prefixStub_ + "simPhi"         + suffix_);
    produces<unsigned>                              (prefixStub_ + "size"           + suffix_);

    produces<std::vector<float> >                   (prefixTrack_ + "px"             + suffix_);
    produces<std::vector<float> >                   (prefixTrack_ + "py"             + suffix_);
    produces<std::vector<float> >                   (prefixTrack_ + "pz"             + suffix_);
    produces<std::vector<float> >                   (prefixTrack_ + "pt"             + suffix_);
    produces<std::vector<float> >                   (prefixTrack_ + "eta"            + suffix_);
    produces<std::vector<float> >                   (prefixTrack_ + "phi"            + suffix_);
    produces<std::vector<float> >                   (prefixTrack_ + "vx"             + suffix_);
    produces<std::vector<float> >                   (prefixTrack_ + "vy"             + suffix_);
    produces<std::vector<float> >                   (prefixTrack_ + "vz"             + suffix_);
    produces<std::vector<float> >                   (prefixTrack_ + "rinv"           + suffix_);
    produces<std::vector<float> >                   (prefixTrack_ + "chi2"           + suffix_);
    produces<std::vector<float> >                   (prefixTrack_ + "ptconsistency"  + suffix_);
    produces<std::vector<unsigned> >                (prefixTrack_ + "nstubs"         + suffix_);
    produces<std::vector<bool> >                    (prefixTrack_ + "isGenuine"      + suffix_);
    produces<std::vector<bool> >                    (prefixTrack_ + "isUnknown"      + suffix_);
    produces<std::vector<bool> >                    (prefixTrack_ + "isCombinatoric" + suffix_);
    produces<std::vector<int> >                     (prefixTrack_ + "tpId"           + suffix_);
    produces<std::vector<int> >                     (prefixTrack_ + "pdgId"          + suffix_);
    produces<std::vector<float> >                   (prefixTrack_ + "simPt"          + suffix_);
    produces<std::vector<float> >                   (prefixTrack_ + "simEta"         + suffix_);
    produces<std::vector<float> >                   (prefixTrack_ + "simPhi"         + suffix_);
    produces<unsigned>                              (prefixTrack_ + "size"           + suffix_);
}

void NTupleStubs::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {
    /// Geometry setup
    /// Set pointers to Geometry
    edm::ESHandle<TrackerGeometry> geometryHandle;
    iSetup.get<TrackerDigiGeometryRecord>().get(geometryHandle);
    theGeometry = geometryHandle.product();

    /// Set pointers to Stacked Modules
    edm::ESHandle<StackedTrackerGeometry> stackedGeometryHandle;
    iSetup.get<StackedTrackerGeometryRecord>().get(stackedGeometryHandle);
    theStackedGeometry = stackedGeometryHandle.product();

    /// Magnetic field setup
    edm::ESHandle<MagneticField> magneticFieldHandle;
    iSetup.get<IdealMagneticFieldRecord>().get(magneticFieldHandle);
    theMagneticField = magneticFieldHandle.product();

    /// Make sure uint32_t == unsigned int
    unsigned checkunsigned = 4294967295;
    uint32_t checkuint32_t = 4294967295;
    assert(checkunsigned == checkuint32_t);

    /// Learn the geometry
    /*
    std::map<uint32_t, uint32_t> moduleIdMap;
    TrackingGeometry::DetUnitContainer::const_iterator trkIterator;
    for(trkIterator = theGeometry->detUnits().begin(); trkIterator != theGeometry->detUnits().end(); ++trkIterator) {
        DetId id = (**trkIterator).geographicalId();
        double r = (**trkIterator).position().perp();
        if ( (**trkIterator).type().isBarrel() &&
             (**trkIterator).type().isTrackerPixel() &&
             (r>20.0) &&
             id.subdetId() == PixelSubdetector::PixelBarrel) {
            //uint32_t lay = PXBDetId(id).layer();
            //uint32_t rod = PXBDetId(id).ladder();
            uint32_t mod = PXBDetId(id).module();
            uint32_t moduleId = getModuleId(id);
            //double zModule = (**trkIterator).position().z();
            //std::cout << "moduleId: " << moduleId << " module: " << mod << " zModule: " << zModule << std::endl;

            if (mod % 2 != 0) { // inner stack
                if (moduleIdMap.find(moduleId) == moduleIdMap.end() )
                    moduleIdMap.insert(std::make_pair(moduleId, 0) );
                else
                    edm::LogError("NTupleStubs") << "Found duplicate: " << moduleId;
            } else { // outer stack
                if (moduleIdMap.find(moduleId) == moduleIdMap.end() )
                    edm::LogError("NTupleStubs") << "Cannot find the inner stack paired with : " << moduleId;
            }

        } else if ( (**trkIterator).type().isEndcap() &&
                    (**trkIterator).type().isTrackerPixel() &&
                    (fabs(r)>20.0) &&
                    id.subdetId() == PixelSubdetector::PixelEndcap) {
            //uint32_t side = PXFDetId(id).side();
            //uint32_t disk = PXFDetId(id).disk();
            //uint32_t ring = PXFDetId(id).ring();
            uint32_t mod  = PXFDetId(id).module();
            uint32_t moduleId = getModuleId(id);
            //double phiModule = (**trkIterator).position().phi();
            //if (phiModule < 0) phiModule += 2*M_PI;
            //std::cout << "moduleId: " << moduleId << " module: " << mod << " phiModule: " << phiModule << std::endl;

            if (mod % 2 != 0) { // inner stack
                if (moduleIdMap.find(moduleId) == moduleIdMap.end() )
                    moduleIdMap.insert(std::make_pair(moduleId, 0) );
                else
                    edm::LogError("NTupleStubs") << "Found duplicate: " << moduleId;
            } else { // outer stack
                if (moduleIdMap.find(moduleId) == moduleIdMap.end() )
                    edm::LogError("NTupleStubs") << "Cannot find the inner stack paired with : " << moduleId;
            }
        }
    }

    std::map<uint32_t, uint32_t> moduleIdMap2;
    StackedTrackerGeometry::StackContainerIterator stkIterator;
    for(stkIterator = theStackedGeometry->stacks().begin(); stkIterator != theStackedGeometry->stacks().end(); ++stkIterator) {
        DetId id0 = (**stkIterator).stackMember(0);
        DetId id1 = (**stkIterator).stackMember(1);
        uint32_t moduleId0 = getModuleId(id0);
        uint32_t moduleId1 = getModuleId(id1);

        if (!(*stkIterator)->Id().isBarrel() && !(*stkIterator)->Id().isEndcap())
            edm::LogError("NTupleStubs") << "Neither barrel nor endcap: " << (*stkIterator)->Id();

        if (moduleIdMap.find(moduleId0) == moduleIdMap.end() )
            edm::LogError("NTupleStubs") << "stack id0 not found in theGeometry: " << moduleId0;

        if (moduleIdMap2.find(moduleId0) == moduleIdMap2.end() )
            moduleIdMap2.insert(std::make_pair(moduleId0, 0) );
        else
            edm::LogError("NTupleStubs") << "Found duplicate: " << moduleId0;

        if (moduleId0 != moduleId1)
            edm::LogError("NTupleStubs") << "ModuleId should be identical for: " << moduleId0 << ", " << moduleId1;

        std::cout << moduleId0 << std::endl;
    }
    std::cout << "moduleIdMap size: " << moduleIdMap.size() << " moduleIdMap2 size: " << moduleIdMap2.size() << std::endl;
    */

}

void NTupleStubs::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    std::auto_ptr<std::vector<float> >                  vc_x                (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vc_y                (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vc_z                (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vc_r                (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vc_phi              (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vc_localx           (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vc_localy           (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vc_coordx           (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vc_coordy           (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vc_surfx            (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vc_surfy            (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vc_surfz            (new std::vector<float>());
    std::auto_ptr<std::vector<unsigned> >               vc_stackId          (new std::vector<unsigned>());
    std::auto_ptr<std::vector<unsigned> >               vc_iLayer           (new std::vector<unsigned>());
    std::auto_ptr<std::vector<unsigned> >               vc_iRing            (new std::vector<unsigned>());
    std::auto_ptr<std::vector<unsigned> >               vc_iSide            (new std::vector<unsigned>());
    std::auto_ptr<std::vector<unsigned> >               vc_iPhi             (new std::vector<unsigned>());
    std::auto_ptr<std::vector<unsigned> >               vc_iZ               (new std::vector<unsigned>());
    std::auto_ptr<std::vector<unsigned> >               vc_iModLayer        (new std::vector<unsigned>());
    std::auto_ptr<std::vector<unsigned> >               vc_iModLadder       (new std::vector<unsigned>());
    std::auto_ptr<std::vector<unsigned> >               vc_iModModule       (new std::vector<unsigned>());
    std::auto_ptr<std::vector<unsigned> >               vc_iModCols         (new std::vector<unsigned>());
    std::auto_ptr<std::vector<unsigned> >               vc_iModRows         (new std::vector<unsigned>());
    std::auto_ptr<std::vector<float> >                  vc_iModPitchX       (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vc_iModPitchY       (new std::vector<float>());
    std::auto_ptr<std::vector<bool> >                   vc_barrel           (new std::vector<bool>());
    std::auto_ptr<std::vector<bool> >                   vc_psmodule         (new std::vector<bool>());
    std::auto_ptr<std::vector<unsigned> >               vc_modId            (new std::vector<unsigned>());
    std::auto_ptr<std::vector<unsigned> >               vc_geoId            (new std::vector<unsigned>());
    std::auto_ptr<std::vector<unsigned> >               vc_stack            (new std::vector<unsigned>());
    std::auto_ptr<std::vector<unsigned> >               vc_width            (new std::vector<unsigned>());
    std::auto_ptr<std::vector<unsigned> >               vc_nhits            (new std::vector<unsigned>());
    std::auto_ptr<std::vector<std::vector<int> > >      vc_hitCols          (new std::vector<std::vector<int> >());
    std::auto_ptr<std::vector<std::vector<int> > >      vc_hitRows          (new std::vector<std::vector<int> >());
    std::auto_ptr<std::vector<std::vector<int> > >      vc_hitADCs          (new std::vector<std::vector<int> >());
    std::auto_ptr<std::vector<std::vector<int> > >      vc_hitChans         (new std::vector<std::vector<int> >());
    std::auto_ptr<std::vector<std::vector<int> > >      vc_hitTrkIds        (new std::vector<std::vector<int> >());
    std::auto_ptr<std::vector<std::vector<unsigned> > > vc_hitEvtIds        (new std::vector<std::vector<unsigned> >());
    std::auto_ptr<std::vector<std::vector<float> > >    vc_hitFracs         (new std::vector<std::vector<float> >());
    std::auto_ptr<std::vector<std::vector<float> > >    vc_hitXs            (new std::vector<std::vector<float> >());
    std::auto_ptr<std::vector<std::vector<float> > >    vc_hitYs            (new std::vector<std::vector<float> >());
    std::auto_ptr<std::vector<std::vector<float> > >    vc_hitZs            (new std::vector<std::vector<float> >());
    std::auto_ptr<std::vector<bool> >                   vc_isGenuine        (new std::vector<bool>());
    std::auto_ptr<std::vector<bool> >                   vc_isUnknown        (new std::vector<bool>());
    std::auto_ptr<std::vector<bool> >                   vc_isCombinatoric   (new std::vector<bool>());
    std::auto_ptr<std::vector<int> >                    vc_tpId             (new std::vector<int>());
    std::auto_ptr<std::vector<int> >                    vc_trkId            (new std::vector<int>());
    std::auto_ptr<std::vector<int> >                    vc_pdgId            (new std::vector<int>());
    std::auto_ptr<std::vector<float> >                  vc_simPt            (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vc_simEta           (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vc_simPhi           (new std::vector<float>());
    std::auto_ptr<unsigned>                             vc_size             (new unsigned(0));

    std::auto_ptr<std::vector<float> >                  vb_x                (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vb_y                (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vb_z                (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vb_r                (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vb_phi              (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vb_coordx           (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vb_coordy           (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vb_dirx             (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vb_diry             (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vb_dirz             (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vb_roughPt          (new std::vector<float>());
    std::auto_ptr<std::vector<unsigned> >               vb_stackId          (new std::vector<unsigned>());
    std::auto_ptr<std::vector<unsigned> >               vb_iLayer           (new std::vector<unsigned>());
    std::auto_ptr<std::vector<unsigned> >               vb_iRing            (new std::vector<unsigned>());
    std::auto_ptr<std::vector<unsigned> >               vb_iSide            (new std::vector<unsigned>());
    std::auto_ptr<std::vector<unsigned> >               vb_iPhi             (new std::vector<unsigned>());
    std::auto_ptr<std::vector<unsigned> >               vb_iZ               (new std::vector<unsigned>());
    std::auto_ptr<std::vector<unsigned> >               vb_iModLayer        (new std::vector<unsigned>());
    std::auto_ptr<std::vector<unsigned> >               vb_iModLadder       (new std::vector<unsigned>());
    std::auto_ptr<std::vector<unsigned> >               vb_iModModule       (new std::vector<unsigned>());
    std::auto_ptr<std::vector<unsigned> >               vb_iModCols         (new std::vector<unsigned>());
    std::auto_ptr<std::vector<unsigned> >               vb_iModRows         (new std::vector<unsigned>());
    std::auto_ptr<std::vector<float> >                  vb_iModPitchX       (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vb_iModPitchY       (new std::vector<float>());
    std::auto_ptr<std::vector<bool> >                   vb_barrel           (new std::vector<bool>());
    std::auto_ptr<std::vector<bool> >                   vb_psmodule         (new std::vector<bool>());
    std::auto_ptr<std::vector<unsigned> >               vb_modId            (new std::vector<unsigned>());
    std::auto_ptr<std::vector<unsigned> >               vb_geoId0           (new std::vector<unsigned>());
    std::auto_ptr<std::vector<unsigned> >               vb_geoId1           (new std::vector<unsigned>());
    std::auto_ptr<std::vector<unsigned> >               vb_clusId0          (new std::vector<unsigned>());
    std::auto_ptr<std::vector<unsigned> >               vb_clusId1          (new std::vector<unsigned>());
    std::auto_ptr<std::vector<unsigned> >               vb_nhits            (new std::vector<unsigned>());
    std::auto_ptr<std::vector<std::vector<int> > >      vb_hitCols          (new std::vector<std::vector<int> >());
    std::auto_ptr<std::vector<std::vector<int> > >      vb_hitRows          (new std::vector<std::vector<int> >());
    std::auto_ptr<std::vector<std::vector<int> > >      vb_hitADCs          (new std::vector<std::vector<int> >());
    std::auto_ptr<std::vector<std::vector<int> > >      vb_hitChans         (new std::vector<std::vector<int> >());
    std::auto_ptr<std::vector<std::vector<int> > >      vb_hitTrkIds        (new std::vector<std::vector<int> >());
    std::auto_ptr<std::vector<std::vector<unsigned> > > vb_hitEvtIds        (new std::vector<std::vector<unsigned> >());
    std::auto_ptr<std::vector<std::vector<float> > >    vb_hitFracs         (new std::vector<std::vector<float> >());
    std::auto_ptr<std::vector<std::vector<float> > >    vb_hitXs            (new std::vector<std::vector<float> >());
    std::auto_ptr<std::vector<std::vector<float> > >    vb_hitYs            (new std::vector<std::vector<float> >());
    std::auto_ptr<std::vector<std::vector<float> > >    vb_hitZs            (new std::vector<std::vector<float> >());
    std::auto_ptr<std::vector<float> >                  vb_trigDist         (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vb_trigOffset       (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vb_trigPos          (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vb_trigBend         (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vb_separation       (new std::vector<float>());
    std::auto_ptr<std::vector<int> >                    vb_detWindow        (new std::vector<int>());
    std::auto_ptr<std::vector<bool> >                   vb_isGenuine        (new std::vector<bool>());
    std::auto_ptr<std::vector<bool> >                   vb_isUnknown        (new std::vector<bool>());
    std::auto_ptr<std::vector<bool> >                   vb_isCombinatoric   (new std::vector<bool>());
    std::auto_ptr<std::vector<int> >                    vb_tpId             (new std::vector<int>());
    std::auto_ptr<std::vector<int> >                    vb_trkId            (new std::vector<int>());
    std::auto_ptr<std::vector<int> >                    vb_pdgId            (new std::vector<int>());
    std::auto_ptr<std::vector<float> >                  vb_simPt            (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vb_simEta           (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vb_simPhi           (new std::vector<float>());
    std::auto_ptr<unsigned>                             vb_size             (new unsigned(0));

    std::auto_ptr<std::vector<float> >                  vt_px               (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vt_py               (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vt_pz               (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vt_pt               (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vt_eta              (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vt_phi              (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vt_vx               (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vt_vy               (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vt_vz               (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vt_rinv             (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vt_chi2             (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vt_ptconsistency    (new std::vector<float>());
    std::auto_ptr<std::vector<unsigned> >               vt_nstubs           (new std::vector<unsigned>());
    std::auto_ptr<std::vector<bool> >                   vt_isGenuine        (new std::vector<bool>());
    std::auto_ptr<std::vector<bool> >                   vt_isUnknown        (new std::vector<bool>());
    std::auto_ptr<std::vector<bool> >                   vt_isCombinatoric   (new std::vector<bool>());
    std::auto_ptr<std::vector<int> >                    vt_tpId             (new std::vector<int>());
    std::auto_ptr<std::vector<int> >                    vt_pdgId            (new std::vector<int>());
    std::auto_ptr<std::vector<float> >                  vt_simPt            (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vt_simEta           (new std::vector<float>());
    std::auto_ptr<std::vector<float> >                  vt_simPhi           (new std::vector<float>());
    std::auto_ptr<unsigned>                             vt_size             (new unsigned(0));

    //__________________________________________________________________________
    /// Get all the handles!

    /// TTCluster, TTStub, TTTrack
    edm::Handle<edmNew::DetSetVector<TTCluster<Ref_PixelDigi_> > > pixelDigiTTClusters;
    iEvent.getByLabel(inputTagClus_, pixelDigiTTClusters);
    edm::Handle<edmNew::DetSetVector<TTStub<Ref_PixelDigi_> > > pixelDigiTTStubs;
    iEvent.getByLabel(inputTagStub_, pixelDigiTTStubs);
    edm::Handle<std::vector<TTTrack<Ref_PixelDigi_> > > pixelDigiTTTracks;
    iEvent.getByLabel(inputTagTrack_, pixelDigiTTTracks);

    /// MC truth association maps
    edm::Handle<TTClusterAssociationMap<Ref_PixelDigi_> > mcAssocTTClusters;
    iEvent.getByLabel(inputTagClusMCAssoc_, mcAssocTTClusters);
    edm::Handle<TTStubAssociationMap<Ref_PixelDigi_> >    mcAssocTTStubs;
    iEvent.getByLabel(inputTagStubMCAssoc_, mcAssocTTStubs);
    edm::Handle<TTTrackAssociationMap<Ref_PixelDigi_> >   mcAssocTTTracks;
    iEvent.getByLabel(inputTagTrackMCAssoc_, mcAssocTTTracks);

    /// SimPixelDigis
    edm::Handle<edm::DetSetVector<PixelDigi> > pixelDigis;
    iEvent.getByLabel(inputTagDigi_, pixelDigis);
    edm::Handle<edm::DetSetVector<PixelDigiSimLink> >  pixelDigiSimLinks;
    iEvent.getByLabel(inputTagDigi_, pixelDigiSimLinks);


    //__________________________________________________________________________
    /// TTClusters
    if (pixelDigiTTClusters.isValid() && mcAssocTTClusters.isValid()) {
        edm::LogInfo("NTupleClusters") << "Size: " << pixelDigiTTClusters->size();

        typedef typename edmNew::DetSetVector<TTCluster<Ref_PixelDigi_> >::const_iterator const_dsv_iter;
        typedef typename edmNew::DetSet      <TTCluster<Ref_PixelDigi_> >::const_iterator const_ds_iter;
        typedef edm::Ref<edmNew::DetSetVector<TTCluster<Ref_PixelDigi_> >, TTCluster<Ref_PixelDigi_> > reference;

        unsigned n = 0;
        for (const_dsv_iter itv = pixelDigiTTClusters->begin(); itv != pixelDigiTTClusters->end(); ++itv) {
            unsigned dsv_detId = itv->detId();
            unsigned dsv_size = itv->size();
            edm::LogInfo("NTupleClusters") << "dsv_detId: " << dsv_detId << " dsv_size: " << dsv_size;

            for (const_ds_iter it = itv->begin(); it != itv->end(); ++it) {
                const StackedTrackerDetId detId = it->getDetId();  // same as stub detId
                uint32_t stackId = detId.rawId();
                unsigned iLayer = detId.iLayer();
                unsigned iRing = detId.iRing();
                unsigned iSide = detId.iSide();
                if (iLayer == 999999 && (iRing == 999999 || iSide == 999999)) {
                    edm::LogError("NTupleClusters") << "Inconsistent detId with iLayer: " << iLayer << ", iRing: " << iRing << ", iSide: " << iSide;
                }
                unsigned iPhi = detId.iPhi();
                unsigned iZ = detId.iZ();
                bool isBarrel = detId.isBarrel();
                bool isEndcap = detId.isEndcap();
                if (isBarrel == isEndcap) {
                    edm::LogError("NTupleClusters") << "Inconsistent detId with isBarrel: " << isBarrel << ", isEndcap: " << isEndcap;
                }

                /// Module ID
                unsigned stack = it->getStackMember();  // inner sensor=0, outer sensor=1
                unsigned width = it->findWidth();  // cluster width in xy
                const DetId geoId = theStackedGeometry->idToDet(detId, stack)->geographicalId();
                unsigned iModuleLayer = getModuleLayer(geoId);
                unsigned iModuleLadder = getModuleLadder(geoId);
                unsigned iModuleModule = getModuleModule(geoId);
                unsigned moduleId = getModuleId(geoId);

                /// Find pixel pitch and topology related information
                const GeomDetUnit* geoUnit = theStackedGeometry->idToDetUnit(detId, stack);
                const PixelGeomDetUnit* pixUnit = dynamic_cast<const PixelGeomDetUnit*>(geoUnit);
                const PixelTopology* pixTopo = dynamic_cast<const PixelTopology*>(&(pixUnit->specificTopology()) );
                //const Bounds& bounds = pixUnit->specificSurface().bounds();

                uint32_t iModuleCols = pixTopo->ncolumns();
                uint32_t iModuleRows = pixTopo->nrows();
                uint32_t iModulePitchX = pixTopo->pitch().first;
                uint32_t iModulePitchY = pixTopo->pitch().second;
                //uint32_t iModuleROCsX = pixTopo->rocsX();
                //uint32_t iModuleROCsY = pixTopo->rocsY();
                //uint32_t iModuleRowsPerROC = pixTopo->rowsperroc();
                //uint32_t iModuleColsPerROC = pixTopo->colsperroc();
                //float detThickness = bounds.thickness();
                //float detLength = bounds.length();
                //float detWidth = bounds.width();
                //float detWidthAtHalfLength = bounds.widthAtHalfLength();

                /// digis a.k.a. hits
                const std::vector<Ref_PixelDigi_>& hits = it->getHits();
                const std::vector<int>& hitCols = it->getCols();
                const std::vector<int>& hitRows = it->getRows();
                std::vector<int> hitADCs;  // should be unsigned, but use int for uniformity
                std::vector<int> hitChans;
                std::vector<int> hitTrkIds;
                std::vector<unsigned> hitEvtIds;
                std::vector<float> hitFracs;
                std::vector<float> hitXs;
                std::vector<float> hitYs;
                std::vector<float> hitZs;
                for (unsigned ii = 0; ii < hits.size(); ++ii) {
                    hitADCs.push_back(hits.at(ii)->adc());
                    hitChans.push_back(hits.at(ii)->channel() );
                    hitTrkIds.push_back(-1);      // to be updated when sim links are found
                    hitEvtIds.push_back(999999);  // to be updated when sim links are found
                    hitFracs.push_back(0.);       // to be updated when sim links are found

                    // Don't save local position as it's nothing more but (row+0.5, col+0.5)
                    //const LocalPoint& hitLocalPosition = theStackedGeometry->findHitLocalPosition(&(*it), ii);
                    const GlobalPoint& hitPosition = theStackedGeometry->findHitGlobalPosition(&(*it), ii);
                    hitXs.push_back(hitPosition.x());
                    hitYs.push_back(hitPosition.y());
                    hitZs.push_back(hitPosition.z());
                }

                const MeasurementPoint& localcoord = it->findAverageLocalCoordinates();
                const LocalPoint&  localposition = theStackedGeometry->findAverageLocalPosition(&(*it));
                const GlobalPoint& position = theStackedGeometry->findAverageGlobalPosition(&(*it));
                const Surface::PositionType& surfposition = theStackedGeometry->idToDet(detId, stack)->position();
                //const Surface::RotationType& surfsrotation = theStackedGeometry->idToDet(detId, stack)->rotation();

                edm::LogInfo("NTupleClusters") << "stackId: " << stackId << " iLayer: " << iLayer << " iRing: " << iRing << " iSide: " << iSide << " iPhi: " << iPhi << " iZ: " << iZ << " isBarrel: " << isBarrel << " geoId: " << geoId.rawId() << " moduleId: " << moduleId << " stack: " << stack << " width: " << width << " position: " << position << " localposition: " << localposition << " localcoord: " << localcoord << " surfposition: " << surfposition;


                vc_x->push_back(position.x());                  // sviret/HL_LHC: CLUS_x
                vc_y->push_back(position.y());                  // sviret/HL_LHC: CLUS_y
                vc_z->push_back(position.z());                  // sviret/HL_LHC: CLUS_z
                vc_r->push_back(position.perp());
                vc_phi->push_back(position.phi());
                vc_localx->push_back(localposition.x());
                vc_localy->push_back(localposition.y());
                vc_coordx->push_back(localcoord.x());           // sviret/HL_LHC: CLUS_strip
                vc_coordy->push_back(localcoord.y());           // sviret/HL_LHC: CLUS_seg
                vc_surfx->push_back(surfposition.x());
                vc_surfy->push_back(surfposition.y());
                vc_surfz->push_back(surfposition.z());
                vc_stackId->push_back(stackId);  // this is stub raw Id
                vc_iLayer->push_back(iLayer);
                vc_iRing->push_back(iRing);
                vc_iSide->push_back(iSide);
                vc_iPhi->push_back(iPhi);
                vc_iZ->push_back(iZ);
                vc_iModLayer->push_back(iModuleLayer);          // sviret/HL_LHC: CLUS_layer
                vc_iModLadder->push_back(iModuleLadder);        // sviret/HL_LHC: CLUS_ladder
                vc_iModModule->push_back(iModuleModule);        // sviret/HL_LHC: CLUS_module
                vc_iModCols->push_back(iModuleCols);            // sviret/HL_LHC: CLUS_PS
                vc_iModRows->push_back(iModuleRows);            // sviret/HL_LHC: CLUS_nrows
                vc_iModPitchX->push_back(iModulePitchX);
                vc_iModPitchY->push_back(iModulePitchY);
                vc_barrel->push_back(isBarrel);
                vc_psmodule->push_back(theStackedGeometry->isPSModule(detId));
                vc_modId->push_back(moduleId);
                vc_geoId->push_back(geoId.rawId());
                vc_stack->push_back(stack);
                vc_width->push_back(width);                     // sviret/HL_LHC: CLUS_nstrips
                vc_nhits->push_back(hits.size());
                vc_hitCols->push_back(hitCols);
                vc_hitRows->push_back(hitRows);
                vc_hitADCs->push_back(hitADCs);
                vc_hitChans->push_back(hitChans);
                vc_hitXs->push_back(hitXs);
                vc_hitYs->push_back(hitYs);
                vc_hitZs->push_back(hitZs);

                /// Sim links for the hits
                if (pixelDigiSimLinks.isValid()) {
                    const edm::DetSet<PixelDigiSimLink>& simlink = (*pixelDigiSimLinks)[geoId];
                    for (unsigned ii = 0; ii < hits.size(); ++ii) {
                        for (edm::DetSet<PixelDigiSimLink>::const_iterator itsim = simlink.data.begin(); itsim != simlink.data.end(); ++itsim) {
                            if (hitChans.at(ii) == (int) itsim->channel()) {
                                hitTrkIds.at(ii) = itsim->SimTrackId();
                                hitEvtIds.at(ii) = itsim->eventId().rawId();
                                hitFracs.at(ii) = itsim->fraction();
                                break;
                            }
                        }
                    }
                }
                vc_hitTrkIds->push_back(hitTrkIds);
                vc_hitEvtIds->push_back(hitEvtIds);
                vc_hitFracs->push_back(hitFracs);

                /// Make the reference to retrieve from MC matching
                reference iref = edmNew::makeRefTo(pixelDigiTTClusters, it);
                vc_isGenuine->push_back(mcAssocTTClusters->isGenuine(iref));
                vc_isUnknown->push_back(mcAssocTTClusters->isUnknown(iref));
                vc_isCombinatoric->push_back(mcAssocTTClusters->isCombinatoric(iref));
                if (vc_isGenuine->back()) {
                    const edm::Ptr<TrackingParticle>& tpptr = mcAssocTTClusters->findTrackingParticlePtr(iref);
                    vc_tpId->push_back(tpptr.key());            // sviret/HL_LHC: CLUS_tp
                    vc_trkId->push_back(tpptr->g4Track_begin()->trackId());
                    vc_pdgId->push_back(tpptr->pdgId());        // sviret/HL_LHC: CLUS_pdgID
                    vc_simPt->push_back(tpptr->pt());           // sviret/HL_LHC: CLUS_ptGEN
                    vc_simEta->push_back(tpptr->eta());
                    vc_simPhi->push_back(tpptr->phi());
                } else {
                    vc_tpId->push_back(-1);
                    vc_trkId->push_back(-1);
                    vc_pdgId->push_back(-1);
                    vc_simPt->push_back(-99.);
                    vc_simEta->push_back(-99.);
                    vc_simPhi->push_back(-99.);
                }

                n++;
            }
        }
        *vc_size = vc_x->size();

    } else {
        edm::LogError("NTupleClusters") << "Cannot get the product: " << inputTagClus_;
    }


    /// TTStubs
    if (pixelDigiTTStubs.isValid() && mcAssocTTStubs.isValid()) {
        edm::LogInfo("NTupleStubs") << "Size: " << pixelDigiTTStubs->size();

        typedef typename edmNew::DetSetVector<TTStub<Ref_PixelDigi_> >::const_iterator const_dsv_iter;
        typedef typename edmNew::DetSet      <TTStub<Ref_PixelDigi_> >::const_iterator const_ds_iter;
        typedef edm::Ref<edmNew::DetSetVector<TTStub<Ref_PixelDigi_> >, TTStub<Ref_PixelDigi_> > reference;

        unsigned n = 0;
        for (const_dsv_iter itv = pixelDigiTTStubs->begin(); itv != pixelDigiTTStubs->end(); ++itv) {
            unsigned dsv_detId = itv->detId();
            unsigned dsv_size = itv->size();
            edm::LogInfo("NTupleStubs") << "dsv_detId: " << dsv_detId << " dsv_size: " << dsv_size;

            for (const_ds_iter it = itv->begin(); it != itv->end(); ++it) {
                const StackedTrackerDetId detId = it->getDetId();
                uint32_t stackId = detId.rawId();
                unsigned iLayer = detId.iLayer();
                unsigned iRing = detId.iRing();
                unsigned iSide = detId.iSide();
                if (iLayer == 999999 && (iRing == 999999 || iSide == 999999)) {
                    edm::LogError("NTupleStubs") << "Inconsistent detId with iLayer: " << iLayer << ", iRing: " << iRing << ", iSide: " << iSide;
                }
                unsigned iPhi = detId.iPhi();
                unsigned iZ = detId.iZ();
                bool isBarrel = detId.isBarrel();
                bool isEndcap = detId.isEndcap();
                if (isBarrel == isEndcap) {
                    edm::LogError("NTupleStubs") << "Inconsistent detId with isBarrel: " << isBarrel << ", isEndcap: " << isEndcap;
                }

                // Get the clusters (inner=0, outer=1)
                //const std::vector<edm::Ref<edmNew::DetSetVector<TTCluster<Ref_PixelDigi_> >, TTCluster<Ref_PixelDigi_> > >& clusterRefs = it->getClusterRefs();
                //if (clusterRefs.size() != 2) {
                //    edm::LogError("NTupleStubs") << "# clusters != 2: " << clusterRefs.size();
                //}
                //const DetId geoId0 = theStackedGeometry->idToDet(clusterRefs[0]->getDetId(), clusterRefs[0]->getStackMember())->geographicalId();
                //const DetId geoId1 = theStackedGeometry->idToDet(clusterRefs[1]->getDetId(), clusterRefs[1]->getStackMember())->geographicalId();
                const DetId geoId0 = theStackedGeometry->idToDet(detId, 0)->geographicalId();
                const DetId geoId1 = theStackedGeometry->idToDet(detId, 1)->geographicalId();


                // Get the cluster indices in this ntuple
                if (vc_geoId->empty()) {
                    edm::LogError("NTupleStubs") << "TTClusters have not been filled!" << std::endl;
                }
                unsigned clusId0 = findById(*vc_geoId, geoId0.rawId());
                unsigned clusId1 = findById(*vc_geoId, geoId1.rawId());

                // For moduleId, we use the inner cluster (geoId0)
                unsigned iModuleLayer = getModuleLayer(geoId0);
                unsigned iModuleLadder = getModuleLadder(geoId0);
                unsigned iModuleModule = getModuleModule(geoId0);
                unsigned moduleId = getModuleId(geoId0);

                const GeomDetUnit* geoUnit = theStackedGeometry->idToDetUnit(detId, 0);
                const PixelGeomDetUnit* pixUnit = dynamic_cast<const PixelGeomDetUnit*>(geoUnit);
                const PixelTopology* pixTopo = dynamic_cast<const PixelTopology*>(&(pixUnit->specificTopology()) );

                uint32_t iModuleCols = pixTopo->ncolumns();  // like # subladders
                uint32_t iModuleRows = pixTopo->nrows();     // like # submodules
                uint32_t iModulePitchX = pixTopo->pitch().first;
                uint32_t iModulePitchY = pixTopo->pitch().second;
                //uint32_t iModuleROCsX = pixTopo->rocsX();
                //uint32_t iModuleROCsY = pixTopo->rocsY();
                //uint32_t iModuleRowsPerROC = pixTopo->rowsperroc();
                //uint32_t iModuleColsPerROC = pixTopo->colsperroc();

                /// digis a.k.a. hits
                const edm::Ref<edmNew::DetSetVector<TTCluster<Ref_PixelDigi_> >, TTCluster<Ref_PixelDigi_> >& clusterRef = it->getClusterRef(0);
                //const DetId testgeoId0 = theStackedGeometry->idToDet(clusterRef->getDetId(), clusterRef->getStackMember())->geographicalId();
                //assert(testgeoId0 == geoId0);
                const std::vector<Ref_PixelDigi_>& hits = clusterRef->getHits();
                const std::vector<int>& hitCols = clusterRef->getCols();
                const std::vector<int>& hitRows = clusterRef->getRows();
                std::vector<int> hitADCs;  // should be unsigned, but use int for uniformity
                std::vector<int> hitChans;
                std::vector<int> hitTrkIds;
                std::vector<unsigned> hitEvtIds;
                std::vector<float> hitFracs;
                std::vector<float> hitXs;
                std::vector<float> hitYs;
                std::vector<float> hitZs;
                for (unsigned ii = 0; ii < hits.size(); ++ii) {
                    hitADCs.push_back(hits.at(ii)->adc());
                    hitChans.push_back(hits.at(ii)->channel() );
                    hitTrkIds.push_back(-1);      // to be updated when sim links are found
                    hitEvtIds.push_back(999999);  // to be updated when sim links are found
                    hitFracs.push_back(0.);       // to be updated when sim links are found

                    // Don't save local position as it's nothing more but (row+0.5, col+0.5)
                    //const LocalPoint& hitLocalPosition = theStackedGeometry->findHitLocalPosition(&(*clusterRef), ii);
                    const GlobalPoint& hitPosition = theStackedGeometry->findHitGlobalPosition(&(*clusterRef), ii);
                    hitXs.push_back(hitPosition.x());
                    hitYs.push_back(hitPosition.y());
                    hitZs.push_back(hitPosition.z());
                }

                // Global position and rough pt
                double magfieldStrength = theMagneticField->inTesla(GlobalPoint(0,0,0)).z();
                double roughPt = theStackedGeometry->findRoughPt(magfieldStrength, &(*it));
                const GlobalPoint& position = theStackedGeometry->findGlobalPosition(&(*it));
                const GlobalVector& direction = theStackedGeometry->findGlobalDirection(&(*it));
                const Surface::PositionType& surfposition0 = theStackedGeometry->idToDet(detId, 0)->position();
                const Surface::PositionType& surfposition1 = theStackedGeometry->idToDet(detId, 1)->position();
                double separation = isBarrel ? fabs(surfposition0.perp() - surfposition1.perp()) : fabs(surfposition0.z() - surfposition1.z());

                const MeasurementPoint& localcoord = clusterRef->findAverageLocalCoordinates();

                edm::LogInfo("NTupleStubs") << "stackId: " << stackId << " iLayer: " << iLayer << " iRing: " << iRing << " iSide: " << iSide << " iPhi: " << iPhi << " iZ: " << iZ << " isBarrel: " << isBarrel << " roughPt: " << roughPt << " position: " << position << " direction: " << direction << " geoId0: " << geoId0.rawId() << " geoId1: " << geoId1.rawId() << " clusId0: " << clusId0 << " clusId1: " << clusId1 << " moduleId: " << moduleId;

                vb_x->push_back(position.x());                  // sviret/HL_LHC: STUB_x
                vb_y->push_back(position.y());                  // sviret/HL_LHC: STUB_y
                vb_z->push_back(position.z());                  // sviret/HL_LHC: STUB_z
                vb_r->push_back(position.perp());
                vb_phi->push_back(position.phi());
                vb_coordx->push_back(localcoord.x());           // sviret/HL_LHC: STUB_strip
                vb_coordy->push_back(localcoord.y());           // sviret/HL_LHC: STUB_seg
                vb_dirx->push_back(direction.x());
                vb_diry->push_back(direction.y());
                vb_dirz->push_back(direction.z());
                vb_roughPt->push_back(roughPt);                 // sviret/HL_LHC: STUB_pt
                vb_stackId->push_back(stackId);
                vb_iLayer->push_back(iLayer);
                vb_iRing->push_back(iRing);
                vb_iSide->push_back(iSide);
                vb_iPhi->push_back(iPhi);
                vb_iZ->push_back(iZ);
                vb_iModLayer->push_back(iModuleLayer);          // sviret/HL_LHC: STUB_layer
                vb_iModLadder->push_back(iModuleLadder);        // sviret/HL_LHC: STUB_ladder
                vb_iModModule->push_back(iModuleModule);        // sviret/HL_LHC: STUB_module
                vb_iModCols->push_back(iModuleCols);
                vb_iModRows->push_back(iModuleRows);
                vb_iModPitchX->push_back(iModulePitchX);
                vb_iModPitchY->push_back(iModulePitchY);
                vb_barrel->push_back(isBarrel);
                vb_psmodule->push_back(theStackedGeometry->isPSModule(detId));
                vb_modId->push_back(moduleId);
                vb_geoId0->push_back(geoId0.rawId());
                vb_geoId1->push_back(geoId1.rawId());
                vb_clusId0->push_back(clusId0);                 // sviret/HL_LHC: STUB_clust1
                vb_clusId1->push_back(clusId1);                 // sviret/HL_LHC: STUB_clust2
                vb_nhits->push_back(hits.size());
                vb_hitCols->push_back(hitCols);                 // sviret/HL_LHC: STUB_pix
                vb_hitRows->push_back(hitRows);
                vb_hitADCs->push_back(hitADCs);
                vb_hitChans->push_back(hitChans);
                vb_hitXs->push_back(hitXs);
                vb_hitYs->push_back(hitYs);
                vb_hitZs->push_back(hitZs);

                vb_trigDist->push_back(it->getTriggerDisplacement());   // (in full-strip units)
                vb_trigOffset->push_back(it->getTriggerOffset());       // sviret/HL_LHC: STUB_cor    (in full-strip units)
                vb_trigPos->push_back(it->getTriggerPosition());        // (in full-strip units)
                vb_trigBend->push_back(it->getTriggerBend());           // sviret/HL_LHC: STUB_deltas (in full-strip units)
                vb_separation->push_back(separation);
                vb_detWindow->push_back(theStackedGeometry->getDetUnitWindow(detId));

                /// Sim links for the hits
                if (pixelDigiSimLinks.isValid()) {
                    const edm::DetSet<PixelDigiSimLink>& simlink = (*pixelDigiSimLinks)[geoId0];
                    for (unsigned ii = 0; ii < hits.size(); ++ii) {
                        for (edm::DetSet<PixelDigiSimLink>::const_iterator itsim = simlink.data.begin(); itsim != simlink.data.end(); ++itsim) {
                            if (hitChans.at(ii) == (int) itsim->channel()) {
                                hitTrkIds.at(ii) = itsim->SimTrackId();
                                hitEvtIds.at(ii) = itsim->eventId().rawId();
                                hitFracs.at(ii) = itsim->fraction();
                                break;
                            }
                        }
                    }
                }
                vb_hitTrkIds->push_back(hitTrkIds);
                vb_hitEvtIds->push_back(hitEvtIds);
                vb_hitFracs->push_back(hitFracs);

                /// Make the reference to retrieve from MC matching
                reference iref = edmNew::makeRefTo(pixelDigiTTStubs, it);
                vb_isGenuine->push_back(mcAssocTTStubs->isGenuine(iref));
                vb_isUnknown->push_back(mcAssocTTStubs->isUnknown(iref));
                vb_isCombinatoric->push_back(mcAssocTTStubs->isCombinatoric(iref));
                if (vb_isGenuine->back()) {
                    const edm::Ptr<TrackingParticle>& tpptr = mcAssocTTStubs->findTrackingParticlePtr(iref);
                    vb_tpId->push_back(tpptr.key());            // sviret/HL_LHC: STUB_tp
                    vb_trkId->push_back(tpptr->g4Track_begin()->trackId());
                    vb_pdgId->push_back(tpptr->pdgId());        // sviret/HL_LHC: STUB_pdgID
                    vb_simPt->push_back(tpptr->pt());
                    vb_simEta->push_back(tpptr->eta());         // sviret/HL_LHC: STUB_etaGEN
                    vb_simPhi->push_back(tpptr->phi());
                } else {
                    vb_tpId->push_back(-1);
                    vb_trkId->push_back(-1);
                    vb_pdgId->push_back(-1);
                    vb_simPt->push_back(-99.);
                    vb_simEta->push_back(-99.);
                    vb_simPhi->push_back(-99.);
                }

                n++;
            }
        }
        *vb_size = vb_x->size();

    } else {
        edm::LogError("NTupleStubs") << "Cannot get the product: " << inputTagStub_;
    }


    /// TTTracks
    if (pixelDigiTTTracks.isValid() && mcAssocTTTracks.isValid()) {
        edm::LogInfo("NTupleTracks") << "Size: " << pixelDigiTTTracks->size();

        typedef typename std::vector<TTTrack<Ref_PixelDigi_> >::const_iterator const_iterator;
        typedef edm::Ptr<TTTrack<Ref_PixelDigi_> > pointer;

        unsigned n = 0;
        for (const_iterator it = pixelDigiTTTracks->begin(); it != pixelDigiTTTracks->end(); ++it) {
            GlobalVector momentum = it->getMomentum();
            GlobalPoint  poca = it->getPOCA();  // point of closest approach

            edm::LogInfo("NTupleTracks") << "pt: " << momentum.perp() << " vz: " << poca.z() << std::endl;

            vt_px->push_back(momentum.x());
            vt_py->push_back(momentum.y());
            vt_pz->push_back(momentum.z());
            vt_pt->push_back(momentum.perp());
            vt_eta->push_back(momentum.eta());
            vt_phi->push_back(momentum.phi());
            vt_vx->push_back(poca.x());
            vt_vy->push_back(poca.y());
            vt_vz->push_back(poca.z());
            vt_rinv->push_back(it->getRInv());
            vt_chi2->push_back(it->getChi2());
            vt_ptconsistency->push_back(it->getStubPtConsistency());
            vt_nstubs->push_back(it->getStubRefs().size());

            /// Make the pointer to retrieve from MC matching
            pointer iptr(pixelDigiTTTracks, it - pixelDigiTTTracks->begin());
            vt_isGenuine->push_back(mcAssocTTTracks->isGenuine(iptr));
            vt_isUnknown->push_back(mcAssocTTTracks->isUnknown(iptr));
            vt_isCombinatoric->push_back(mcAssocTTTracks->isCombinatoric(iptr));
            if (vt_isGenuine->back()) {
                edm::Ptr<TrackingParticle> tpptr = mcAssocTTTracks->findTrackingParticlePtr(iptr);
                vt_tpId->push_back(tpptr.key());
                vt_pdgId->push_back(tpptr->pdgId());
                vt_simPt->push_back(tpptr->pt());
                vt_simEta->push_back(tpptr->eta());
                vt_simPhi->push_back(tpptr->phi());
            } else {
                vt_tpId->push_back(-1);
                vt_pdgId->push_back(-99);
                vt_simPt->push_back(-99.);
                vt_simEta->push_back(-99.);
                vt_simPhi->push_back(-99.);
            }

            n++;
        }
        *vt_size = vt_px->size();

    } else {
        edm::LogError("NTupleTracks") << "Cannot get the product: " << inputTagTrack_;
    }


    /// SimPixelDigis -- skipped for now
    /// Add position, modId, clusId, stubId and sim info?
    //if (pixelDigis.isValid() && pixelDigiSimLinks.isValid()) {
    //    edm::LogInfo("NTupleSimPixelDigis") << "Size: " << pixelDigis->size();
    //
    //    typedef edm::DetSetVector<PixelDigi>::const_iterator const_dsv_iterator;
    //    typedef edm::DetSet      <PixelDigi>::const_iterator const_ds_iterator;
    //    unsigned n = 0;
    //    for (const_dsv_iterator itv = pixelDigis->begin(); itv != pixelDigis->end(); ++itv) {
    //        uint32_t rawId = itv->detId();
    //        const DetId detId(rawId);
    //        DetId::Detector det = detId.det();  // Tracker=1,Muon=2,Ecal=3,Hcal=4,Calo=5,Forward=6
    //        int subdetId = detId.subdetId();  // PXB=1,PXF=2,...
    //        if (subdetId != PixelSubdetector::PixelBarrel && subdetId != PixelSubdetector::PixelEndcap)
    //            continue;
    //
    //        edm::LogInfo("NTupleSimPixelDigis") << "rawId: " << detId.rawId() << " det: " << det << " subdetId: " << subdetId << " size: " << itv->size();
    //
    //        for (const_ds_iterator it = itv->begin(); it != itv->end(); ++it) {
    //            int row = it->row();
    //            int col = it->column();
    //            int channel = it->channel();
    //            unsigned short adc = it->adc();
    //            edm::LogInfo("NTupleSimPixelDigis") << "row: " << row << " col: " << col << " channel: " << channel << " adc: " << adc;
    //        }
    //
    //        n++;
    //    }
    //
    //} else {
    //    edm::LogError("NTupleSimPixelDigis") << "Cannot get the product: " << inputTagPixelDigi_;
    //}


    //__________________________________________________________________________
    iEvent.put(vc_x             , prefixClus_ + "x"              + suffix_);
    iEvent.put(vc_y             , prefixClus_ + "y"              + suffix_);
    iEvent.put(vc_z             , prefixClus_ + "z"              + suffix_);
    iEvent.put(vc_r             , prefixClus_ + "r"              + suffix_);
    iEvent.put(vc_phi           , prefixClus_ + "phi"            + suffix_);
    iEvent.put(vc_localx        , prefixClus_ + "localx"         + suffix_);
    iEvent.put(vc_localy        , prefixClus_ + "localy"         + suffix_);
    iEvent.put(vc_coordx        , prefixClus_ + "coordx"         + suffix_);
    iEvent.put(vc_coordy        , prefixClus_ + "coordy"         + suffix_);
    iEvent.put(vc_surfx         , prefixClus_ + "surfx"          + suffix_);
    iEvent.put(vc_surfy         , prefixClus_ + "surfy"          + suffix_);
    iEvent.put(vc_surfz         , prefixClus_ + "surfz"          + suffix_);
    iEvent.put(vc_stackId       , prefixClus_ + "stackId"        + suffix_);
    iEvent.put(vc_iLayer        , prefixClus_ + "iLayer"         + suffix_);
    iEvent.put(vc_iRing         , prefixClus_ + "iRing"          + suffix_);
    iEvent.put(vc_iSide         , prefixClus_ + "iSide"          + suffix_);
    iEvent.put(vc_iPhi          , prefixClus_ + "iPhi"           + suffix_);
    iEvent.put(vc_iZ            , prefixClus_ + "iZ"             + suffix_);
    iEvent.put(vc_iModLayer     , prefixClus_ + "iModLayer"      + suffix_);
    iEvent.put(vc_iModLadder    , prefixClus_ + "iModLadder"     + suffix_);
    iEvent.put(vc_iModModule    , prefixClus_ + "iModModule"     + suffix_);
    iEvent.put(vc_iModCols      , prefixClus_ + "iModCols"       + suffix_);
    iEvent.put(vc_iModRows      , prefixClus_ + "iModRows"       + suffix_);
    iEvent.put(vc_iModPitchX    , prefixClus_ + "iModPitchX"     + suffix_);
    iEvent.put(vc_iModPitchY    , prefixClus_ + "iModPitchY"     + suffix_);
    iEvent.put(vc_barrel        , prefixClus_ + "barrel"         + suffix_);
    iEvent.put(vc_psmodule      , prefixClus_ + "psmodule"       + suffix_);
    iEvent.put(vc_modId         , prefixClus_ + "modId"          + suffix_);
    iEvent.put(vc_geoId         , prefixClus_ + "geoId"          + suffix_);
    iEvent.put(vc_stack         , prefixClus_ + "stack"          + suffix_);
    iEvent.put(vc_width         , prefixClus_ + "width"          + suffix_);
    iEvent.put(vc_nhits         , prefixClus_ + "nhits"          + suffix_);
    iEvent.put(vc_hitCols       , prefixClus_ + "hitCols"        + suffix_);
    iEvent.put(vc_hitRows       , prefixClus_ + "hitRows"        + suffix_);
    iEvent.put(vc_hitADCs       , prefixClus_ + "hitADCs"        + suffix_);
    iEvent.put(vc_hitChans      , prefixClus_ + "hitChans"       + suffix_);
    iEvent.put(vc_hitTrkIds     , prefixClus_ + "hitTrkIds"      + suffix_);
    iEvent.put(vc_hitEvtIds     , prefixClus_ + "hitEvtIds"      + suffix_);
    iEvent.put(vc_hitFracs      , prefixClus_ + "hitFracs"       + suffix_);
    iEvent.put(vc_hitXs         , prefixClus_ + "hitXs"          + suffix_);
    iEvent.put(vc_hitYs         , prefixClus_ + "hitYs"          + suffix_);
    iEvent.put(vc_hitZs         , prefixClus_ + "hitZs"          + suffix_);
    iEvent.put(vc_isGenuine     , prefixClus_ + "isGenuine"      + suffix_);
    iEvent.put(vc_isUnknown     , prefixClus_ + "isUnknown"      + suffix_);
    iEvent.put(vc_isCombinatoric, prefixClus_ + "isCombinatoric" + suffix_);
    iEvent.put(vc_tpId          , prefixClus_ + "tpId"           + suffix_);
    iEvent.put(vc_trkId         , prefixClus_ + "trkId"          + suffix_);
    iEvent.put(vc_pdgId         , prefixClus_ + "pdgId"          + suffix_);
    iEvent.put(vc_simPt         , prefixClus_ + "simPt"          + suffix_);
    iEvent.put(vc_simEta        , prefixClus_ + "simEta"         + suffix_);
    iEvent.put(vc_simPhi        , prefixClus_ + "simPhi"         + suffix_);
    iEvent.put(vc_size          , prefixClus_ + "size"           + suffix_);

    iEvent.put(vb_x             , prefixStub_ + "x"              + suffix_);
    iEvent.put(vb_y             , prefixStub_ + "y"              + suffix_);
    iEvent.put(vb_z             , prefixStub_ + "z"              + suffix_);
    iEvent.put(vb_r             , prefixStub_ + "r"              + suffix_);
    iEvent.put(vb_phi           , prefixStub_ + "phi"            + suffix_);
    iEvent.put(vb_coordx        , prefixStub_ + "coordx"         + suffix_);
    iEvent.put(vb_coordy        , prefixStub_ + "coordy"         + suffix_);
    iEvent.put(vb_dirx          , prefixStub_ + "dirx"           + suffix_);
    iEvent.put(vb_diry          , prefixStub_ + "diry"           + suffix_);
    iEvent.put(vb_dirz          , prefixStub_ + "dirz"           + suffix_);
    iEvent.put(vb_roughPt       , prefixStub_ + "roughPt"        + suffix_);
    iEvent.put(vb_stackId       , prefixStub_ + "stackId"        + suffix_);
    iEvent.put(vb_iLayer        , prefixStub_ + "iLayer"         + suffix_);
    iEvent.put(vb_iRing         , prefixStub_ + "iRing"          + suffix_);
    iEvent.put(vb_iSide         , prefixStub_ + "iSide"          + suffix_);
    iEvent.put(vb_iPhi          , prefixStub_ + "iPhi"           + suffix_);
    iEvent.put(vb_iZ            , prefixStub_ + "iZ"             + suffix_);
    iEvent.put(vb_iModLayer     , prefixStub_ + "iModLayer"      + suffix_);
    iEvent.put(vb_iModLadder    , prefixStub_ + "iModLadder"     + suffix_);
    iEvent.put(vb_iModModule    , prefixStub_ + "iModModule"     + suffix_);
    iEvent.put(vb_iModCols      , prefixStub_ + "iModCols"       + suffix_);
    iEvent.put(vb_iModRows      , prefixStub_ + "iModRows"       + suffix_);
    iEvent.put(vb_iModPitchX    , prefixStub_ + "iModPitchX"     + suffix_);
    iEvent.put(vb_iModPitchY    , prefixStub_ + "iModPitchY"     + suffix_);
    iEvent.put(vb_barrel        , prefixStub_ + "barrel"         + suffix_);
    iEvent.put(vb_psmodule      , prefixStub_ + "psmodule"       + suffix_);
    iEvent.put(vb_modId         , prefixStub_ + "modId"          + suffix_);
    iEvent.put(vb_geoId0        , prefixStub_ + "geoId0"         + suffix_);
    iEvent.put(vb_geoId1        , prefixStub_ + "geoId1"         + suffix_);
    iEvent.put(vb_clusId0       , prefixStub_ + "clusId0"        + suffix_);
    iEvent.put(vb_clusId1       , prefixStub_ + "clusId1"        + suffix_);
    iEvent.put(vb_nhits         , prefixStub_ + "nhits"          + suffix_);
    iEvent.put(vb_hitCols       , prefixStub_ + "hitCols"        + suffix_);
    iEvent.put(vb_hitRows       , prefixStub_ + "hitRows"        + suffix_);
    iEvent.put(vb_hitADCs       , prefixStub_ + "hitADCs"        + suffix_);
    iEvent.put(vb_hitChans      , prefixStub_ + "hitChans"       + suffix_);
    iEvent.put(vb_hitTrkIds     , prefixStub_ + "hitTrkIds"      + suffix_);
    iEvent.put(vb_hitEvtIds     , prefixStub_ + "hitEvtIds"      + suffix_);
    iEvent.put(vb_hitFracs      , prefixStub_ + "hitFracs"       + suffix_);
    iEvent.put(vb_hitXs         , prefixStub_ + "hitXs"          + suffix_);
    iEvent.put(vb_hitYs         , prefixStub_ + "hitYs"          + suffix_);
    iEvent.put(vb_hitZs         , prefixStub_ + "hitZs"          + suffix_);
    iEvent.put(vb_trigDist      , prefixStub_ + "trigDist"       + suffix_);
    iEvent.put(vb_trigOffset    , prefixStub_ + "trigOffset"     + suffix_);
    iEvent.put(vb_trigPos       , prefixStub_ + "trigPos"        + suffix_);
    iEvent.put(vb_trigBend      , prefixStub_ + "trigBend"       + suffix_);
    iEvent.put(vb_separation    , prefixStub_ + "separation"     + suffix_);
    iEvent.put(vb_detWindow     , prefixStub_ + "detWindow"      + suffix_);
    iEvent.put(vb_isGenuine     , prefixStub_ + "isGenuine"      + suffix_);
    iEvent.put(vb_isUnknown     , prefixStub_ + "isUnknown"      + suffix_);
    iEvent.put(vb_isCombinatoric, prefixStub_ + "isCombinatoric" + suffix_);
    iEvent.put(vb_tpId          , prefixStub_ + "tpId"           + suffix_);
    iEvent.put(vb_trkId         , prefixStub_ + "trkId"          + suffix_);
    iEvent.put(vb_pdgId         , prefixStub_ + "pdgId"          + suffix_);
    iEvent.put(vb_simPt         , prefixStub_ + "simPt"          + suffix_);
    iEvent.put(vb_simEta        , prefixStub_ + "simEta"         + suffix_);
    iEvent.put(vb_simPhi        , prefixStub_ + "simPhi"         + suffix_);
    iEvent.put(vb_size          , prefixStub_ + "size"           + suffix_);

    iEvent.put(vt_px            , prefixTrack_ + "px"             + suffix_);
    iEvent.put(vt_py            , prefixTrack_ + "py"             + suffix_);
    iEvent.put(vt_pz            , prefixTrack_ + "pz"             + suffix_);
    iEvent.put(vt_pt            , prefixTrack_ + "pt"             + suffix_);
    iEvent.put(vt_eta           , prefixTrack_ + "eta"            + suffix_);
    iEvent.put(vt_phi           , prefixTrack_ + "phi"            + suffix_);
    iEvent.put(vt_vx            , prefixTrack_ + "vx"             + suffix_);
    iEvent.put(vt_vy            , prefixTrack_ + "vy"             + suffix_);
    iEvent.put(vt_vz            , prefixTrack_ + "vz"             + suffix_);
    iEvent.put(vt_rinv          , prefixTrack_ + "rinv"           + suffix_);
    iEvent.put(vt_chi2          , prefixTrack_ + "chi2"           + suffix_);
    iEvent.put(vt_ptconsistency , prefixTrack_ + "ptconsistency"  + suffix_);
    iEvent.put(vt_nstubs        , prefixTrack_ + "nstubs"         + suffix_);
    iEvent.put(vt_isGenuine     , prefixTrack_ + "isGenuine"      + suffix_);
    iEvent.put(vt_isUnknown     , prefixTrack_ + "isUnknown"      + suffix_);
    iEvent.put(vt_isCombinatoric, prefixTrack_ + "isCombinatoric" + suffix_);
    iEvent.put(vt_tpId          , prefixTrack_ + "tpId"           + suffix_);
    iEvent.put(vt_pdgId         , prefixTrack_ + "pdgId"          + suffix_);
    iEvent.put(vt_simPt         , prefixTrack_ + "simPt"          + suffix_);
    iEvent.put(vt_simEta        , prefixTrack_ + "simEta"         + suffix_);
    iEvent.put(vt_simPhi        , prefixTrack_ + "simPhi"         + suffix_);
    iEvent.put(vt_size          , prefixTrack_ + "size"           + suffix_);
}
