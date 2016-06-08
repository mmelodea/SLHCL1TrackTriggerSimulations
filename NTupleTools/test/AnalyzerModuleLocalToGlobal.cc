#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
//#include "FWCore/ServiceRegistry/interface/Service.h"
//#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerDetUnit.h"
#include "DataFormats/SiPixelDetId/interface/StackedTrackerDetId.h"
//#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"

#include "SLHCL1TrackTriggerSimulations/NTupleTools/interface/ModuleIdFunctor.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include "TString.h"

class LinearRegression {
  public:
    LinearRegression() : Sx(0.), Sy(0.), Sxx(0.), Sxy(0.), Syy(0.), n(0) {}
    ~LinearRegression() {}

    void fill(double x, double y) {
        Sx += x;
        Sy += y;
        Sxx += x * x;
        Sxy += x * y;
        Syy += y * y;
        n += 1;
    }

    void compute(double& alpha, double& beta) {
        double nn = n;
        alpha = (Sy * Sxx - Sx * Sxy) / (nn * Sxx - Sx * Sx);
        beta = (nn * Sxy - Sx * Sy) / (nn * Sxx - Sx * Sx);
    }

  private:
    double Sx, Sy, Sxx, Sxy, Syy;
    int n;
};


class AnalyzerModuleLocalToGlobal : public edm::EDAnalyzer {
  public:
    /// Constructor/destructor
    explicit AnalyzerModuleLocalToGlobal(const edm::ParameterSet&);
    virtual ~AnalyzerModuleLocalToGlobal();

  private:
    virtual void beginRun(const edm::Run&, const edm::EventSetup&);
    virtual void endRun(const edm::Run&, const edm::EventSetup&);

    virtual void beginJob();
    virtual void endJob();
    virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  private:
    /// For event setup
    const TrackerGeometry * theGeometry;
    const StackedTrackerGeometry * theStackedGeometry;

    /// Maps
    std::map<uint32_t, uint32_t> moduleId0ToStackId;
    std::map<uint32_t, uint32_t> moduleId1ToStackId;
    std::map<uint32_t, uint32_t> moduleId0ToGeoId;
    std::map<uint32_t, uint32_t> moduleId1ToGeoId;

    /// Configurations
    std::string csvfile_;
    int verbose_;
};


AnalyzerModuleLocalToGlobal::AnalyzerModuleLocalToGlobal(const edm::ParameterSet& iConfig)
: csvfile_(iConfig.getParameter<std::string>("csv") ),
  verbose_(iConfig.getParameter<int>("verbosity") ) {}

AnalyzerModuleLocalToGlobal::~AnalyzerModuleLocalToGlobal() {}


// Here we make the layout of the detector
void AnalyzerModuleLocalToGlobal::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {
    /// Geometry setup
    edm::ESHandle<TrackerGeometry> geometryHandle;
    iSetup.get<TrackerDigiGeometryRecord>().get(geometryHandle);
    theGeometry = geometryHandle.product();

    edm::ESHandle<StackedTrackerGeometry> stackedGeometryHandle;
    iSetup.get<StackedTrackerGeometryRecord>().get(stackedGeometryHandle);
    theStackedGeometry = stackedGeometryHandle.product();

    /// Clear maps
    moduleId0ToStackId.clear();
    moduleId1ToStackId.clear();
    moduleId0ToGeoId.clear();
    moduleId1ToGeoId.clear();

    /// Prepare detId -> moduleId
    ModuleIdFunctor getModuleId;

    /// Loop over the detector elements
    StackedTrackerGeometry::StackContainerIterator  stkIterator;
    for (stkIterator = theStackedGeometry->stacks().begin();
         stkIterator != theStackedGeometry->stacks().end();
         ++stkIterator) {

        StackedTrackerDetUnit* stackDetUnit = *stkIterator;
        StackedTrackerDetId stackDetId = stackDetUnit->Id();
        assert(stackDetUnit == theStackedGeometry->idToStack(stackDetId));

        /// GeomDet and GeomDetUnit are needed to access each
        /// DetId and topology and geometric features
        /// Convert to specific DetId
        const GeomDet* det0 = theStackedGeometry->idToDet(stackDetId, 0);
        const GeomDet* det1 = theStackedGeometry->idToDet(stackDetId, 1);
        //const GeomDetUnit* detUnit0 = theStackedGeometry->idToDetUnit(stackDetId, 0);
        //const GeomDetUnit* detUnit1 = theStackedGeometry->idToDetUnit(stackDetId, 1);

        const DetId geoId0 = det0->geographicalId();
        const DetId geoId1 = det1->geographicalId();

        uint32_t moduleId0 = getModuleId(geoId0);
        uint32_t moduleId1 = getModuleId(geoId1);
        assert(moduleId0 == moduleId1);

        if (moduleId0ToGeoId.find(moduleId0) == moduleId0ToGeoId.end() )
            moduleId0ToGeoId.insert(std::make_pair(moduleId0, geoId0.rawId()) );
        else
            std::cout << "Error: This pair already exists in map! moduleId0: " << moduleId0 << " geoId0: " << geoId0.rawId() << " existing value in map: " << moduleId0ToGeoId.at(moduleId0) << std::endl;

        if (moduleId1ToGeoId.find(moduleId1) == moduleId1ToGeoId.end() )
            moduleId1ToGeoId.insert(std::make_pair(moduleId1, geoId1.rawId()) );
        else
            std::cout << "Error: This pair already exists in map! moduleId1: " << moduleId1 << " geoId1: " << geoId1.rawId() << " existing value in map: " << moduleId1ToGeoId.at(moduleId1) << std::endl;

        if (moduleId0ToStackId.find(moduleId0) == moduleId0ToStackId.end() )
            moduleId0ToStackId.insert(std::make_pair(moduleId0, stackDetId.rawId()) );
        else
            std::cout << "Error: This pair already exists in map! moduleId0: " << moduleId0 << " stackId: " << stackDetId.rawId() << " existing value in map: " << moduleId0ToStackId.at(moduleId0) << std::endl;

        if (moduleId1ToStackId.find(moduleId1) == moduleId1ToStackId.end() )
            moduleId1ToStackId.insert(std::make_pair(moduleId1, stackDetId.rawId()) );
        else
            std::cout << "Error: This pair already exists in map! moduleId1: " << moduleId1 << " stackId: " << stackDetId.rawId() << " existing value in map: " << moduleId1ToStackId.at(moduleId1) << std::endl;

    }  // end loop over stacks

}

void AnalyzerModuleLocalToGlobal::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {}


void AnalyzerModuleLocalToGlobal::beginJob() {}


// Here we write into the text file
void AnalyzerModuleLocalToGlobal::endJob() {
    std::cout << "Map size: " << moduleId0ToStackId.size() << ", " << moduleId1ToStackId.size()
              << ", " << moduleId0ToGeoId.size() << ", " << moduleId1ToGeoId.size() << std::endl;

    // Open text file
    ofstream csvfile(csvfile_);

    // Write text file
    //int i=0;
    csvfile << "moduleId/I, chipId/I, x_phi0/D, x_phi/D, x_z0/D, x_z/D, x_r0/D, x_r/D" << std::endl;

    for (const auto& kv : moduleId0ToStackId) {
        //i++;  if (i>20)  break;

        const uint32_t moduleId(kv.first);
        const StackedTrackerDetId stackDetId(kv.second);

        for (unsigned i=0; i<1; ++i) {  // only the bottom sensor
            const PixelGeomDetUnit* pixUnit = dynamic_cast<const PixelGeomDetUnit*>(theStackedGeometry->idToDetUnit(stackDetId, i));
            const PixelTopology* pixTopo = dynamic_cast<const PixelTopology*>(&(pixUnit->specificTopology()) );

            int nrows = pixTopo->nrows();
            int ncols = pixTopo->ncolumns();

            // Get phi, z, r conversions
            unsigned nchips = 8;
            unsigned nstrips = 256;
            assert(nchips * nstrips == 2048);
            unsigned nsegments = ncols;
            for (unsigned ichip=0; ichip<nchips; ichip++) {
                LinearRegression phiRegression;
                LinearRegression zRegression;
                LinearRegression rRegression;

                if (stackDetId.isBarrel()) {
                    for (unsigned istrip=0; istrip<nstrips; istrip++) {
                        /// Add 0.5 to get the center of the pixel
                        const MeasurementPoint mp(0.5 + (ichip*nstrips + istrip)*0.5, 0);
                        const GlobalPoint& gp = pixUnit->surface().toGlobal(pixUnit->topology().localPosition(mp));
                        phiRegression.fill(istrip, gp.phi());
                        rRegression.fill(istrip, gp.perp());
                    }

                    for (unsigned isegment=0; isegment<nsegments; isegment++) {
                        /// Add 0.5 to get the center of the pixel
                        const MeasurementPoint mp(0, 0.5 + isegment);
                        const GlobalPoint& gp = pixUnit->surface().toGlobal(pixUnit->topology().localPosition(mp));
                        zRegression.fill(isegment, gp.z());
                    }

                } else {
                    for (unsigned istrip=0; istrip<nstrips; istrip++) {
                        /// Add 0.5 to get the center of the pixel
                        const MeasurementPoint mp(0.5 + (ichip*nstrips + istrip)*0.5, 0);
                        const GlobalPoint& gp = pixUnit->surface().toGlobal(pixUnit->topology().localPosition(mp));
                        phiRegression.fill(istrip, gp.phi());
                        zRegression.fill(istrip, gp.z());
                    }

                    for (unsigned isegment=0; isegment<nsegments; isegment++) {
                        /// Add 0.5 to get the center of the pixel
                        const MeasurementPoint mp(0, 0.5 + isegment);
                        const GlobalPoint& gp = pixUnit->surface().toGlobal(pixUnit->topology().localPosition(mp));
                        rRegression.fill(isegment, gp.perp());
                    }
                }

                double x_phi0 = 0., x_phi = 0., x_z0 = 0., x_z = 0., x_r0 = 0., x_r = 0.;
                phiRegression.compute(x_phi0, x_phi);
                zRegression.compute(x_z0, x_z);
                rRegression.compute(x_r0, x_r);

                // Positions are in unit of centimeter
                csvfile.unsetf(std::ios_base::floatfield);
                csvfile << moduleId << ", " << ichip << ", "
                        << std::scientific << std::setprecision(12)
                        << x_phi0 << ", " << x_phi << ", "
                        << x_z0   << ", " << x_z   << ", "
                        << x_r0   << ", " << x_r
                        << std::endl;
            }
        }
    }

    // Close text file
    csvfile.close();

    std::cout << ">>> " << csvfile_ << " is written." << std::endl;
}


// ANALYZE
void AnalyzerModuleLocalToGlobal::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    // Do nothing
}

// DEFINE THIS AS A PLUG-IN
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(AnalyzerModuleLocalToGlobal);
