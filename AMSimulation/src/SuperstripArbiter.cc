#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/SuperstripArbiter.h"

#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/Helper.h"
using namespace slhcl1tt;

#include <algorithm>
#include <cassert>
#include <iostream>
#include <stdexcept>

static const unsigned MAX_NSTRIPS   = 1024;
static const unsigned MAX_NSEGMENTS = 32;
static const unsigned MODULE_NBITS  = 8;

namespace {
bool is_power_of_two(int x) {
    return x > 0 && !(x & (x-1));
}

unsigned most_sig_bit(unsigned v) {
    unsigned r = 0;
    while (v >>= 1)
        r++;
    return r;
}

unsigned round_to_uint(float x) {
    return std::floor(x + 0.5);
}
}


// _____________________________________________________________________________
SuperstripArbiter::SuperstripArbiter()
: sstype_(SuperstripType::UNKNOWN),
  nsuperstripsPerLayer_(999999),
  useGlobalCoord_(false),
  fixedwidth_nstrips_(0),
  fixedwidth_nz_(0),
  fixedwidth_bit_lshift1_(0),
  fixedwidth_bit_lshift2_(0),
  fixedwidth_bit_rshift1_(0),
  fixedwidth_bit_rshift2_(0),
  projective_nx_(0),
  projective_nz_(0),
  projective_max_nx_(0),
  fountain_sf_(0.),
  fountain_nz_(0),
  fountain_max_nx_(0) {

    // phiWidths for 6 barrel layer [0..5], 5 +z endcap disks [6..10], 5 -z endcap disks [11..15]
    // CUIDADO: dummy values are used for the endcap layers
    phiWidths_ = {
        0.00381*2, 0.00439, 0.00459, 0.00485, 0.00523, 0.00575,
        9.99999, 9.99999, 9.99999, 9.99999, 9.99999,
        9.99999, 9.99999, 9.99999, 9.99999, 9.99999
    };
    assert(phiWidths_.size() == 16);

    phiWidths_opt_lowpt_ = {
        0.01091, 0.00518, 0.00409, 0.00309, 0.00275, 0.00571,
        9.99999, 9.99999, 9.99999, 9.99999, 9.99999,
        9.99999, 9.99999, 9.99999, 9.99999, 9.99999
    };
    assert(phiWidths_opt_lowpt_.size() == 16);

    phiWidths_opt_highpt_ = {
        0.01079, 0.00612, 0.00439, 0.00330, 0.00260, 0.00450,
        9.99999, 9.99999, 9.99999, 9.99999, 9.99999,
        9.99999, 9.99999, 9.99999, 9.99999, 9.99999
    };
    assert(phiWidths_opt_highpt_.size() == 16);

    // Average radii [cm] in the 6 barrel layers
    rMeans_ = {
        22.5913, 35.4772, 50.5402, 68.3101, 88.5002, 107.71
    };
    assert(rMeans_.size() == 6);

    // Strip pitchs [um] in the 6 barrel layers
    float pitches[6] = {100., 100., 100., 90., 90., 90.};
    // Module thicknesses [cm] in the 6 barrel layers
    float thicknesses[6] = {0.26, 0.16, 0.16, 0.18, 0.18, 0.18};

    // Constants used for global phi corrections
    // r * dphi / dr ~= pitch * ds / thickness
    // dphi ~= (pitch [um] / thickness [cm] / r [cm]) * (1e2/1e6) [cm/um] * ds * dr [cm]
    // --> drCorr constant = pitch * thickness /r * 1e-4
    drCorrs_.clear();
    drCorrs_.resize(6);
    for (unsigned i=0; i<6; ++i) {
        drCorrs_[i] = pitches[i] / thicknesses[i] / rMeans_[i] * 1e-4;
    }
}

// _____________________________________________________________________________
void SuperstripArbiter::setDefinition(TString definition, unsigned tt, const TriggerTowerMap* ttmap) {

    // Parse the definition string
    if (definition.Contains("_")) {
        unsigned pos = definition.Index("_");
        unsigned len = definition.Length();
        TString token1 = definition(0,2);
        TString token2 = definition(2,pos-2);
        TString token3 = definition(pos+1,2);
	// TString token4 = definition(pos+3,len-(pos+3));
        TString token4 = definition(pos+3,1);
	
	TString token5, token6;

	unsigned posL5 = definition.Contains("L5x");
	if (posL5>0) token5  = definition(posL5+3,1);
	unsigned posL10 = definition.Contains("L10x");
	if (posL10>0) token6  = definition(posL10+4,1);
	
        if (token1 == "ss") {
            sstype_ = SuperstripType::FIXEDWIDTH;
        } else if (token1 == "nx") {
            sstype_ = SuperstripType::PROJECTIVE;
        } else if (token1 == "sf") {
            sstype_ = SuperstripType::FOUNTAIN;
        } else if (token1 == "op") {
            sstype_ = SuperstripType::FOUNTAINOPT;
        } else {
            throw std::invalid_argument("Incorrect superstrip definition.");
        }

        token2.ReplaceAll("p", ".");
        if (token2.IsFloat()) {
            // Do nothing
        } else {
            throw std::invalid_argument("Incorrect superstrip definition.");
        }
        float token2f = token2.Atof();

        if (token3 == "nz") {
            // Do nothing
        } else {
            throw std::invalid_argument("Incorrect superstrip definition.");
        }

        if (token4.IsFloat()) {
            // Do nothing
        } else {
            throw std::invalid_argument("Incorrect superstrip definition.");
        }
        float token4f = token4.Atof();

	if (token5.IsFloat()) {
            // Do nothing
        } else {
	  token5="1";
        }
        float token5f = token5.Atof();

	if (token6.IsFloat()) {
            // Do nothing
        } else {
	  token6="1";
        }
        float token6f = token6.Atof();

        switch (sstype_) {
        case SuperstripType::FIXEDWIDTH:
            fixedwidth_nstrips_     = token2f;
            fixedwidth_nz_          = token4f;
            useGlobalCoord_         = false;

            if (!is_power_of_two(fixedwidth_nstrips_) || !is_power_of_two(fixedwidth_nz_) ||
                fixedwidth_nstrips_ == 0 || fixedwidth_nz_ == 0 ||
                fixedwidth_nstrips_ > MAX_NSTRIPS || fixedwidth_nz_ > MAX_NSEGMENTS)
                throw std::invalid_argument("Incorrect fixedwidth superstrip definition.");
            break;

        case SuperstripType::PROJECTIVE:
            projective_nx_  = token2f;
            projective_nz_  = token4f;
            useGlobalCoord_ = true;

            if (projective_nx_ == 0 || projective_nz_ == 0)
                throw std::invalid_argument("Incorrect projective superstrip definition.");
            break;

        case SuperstripType::FOUNTAIN:
            fountain_sf_    = token2f;
            fountain_nz_    = token4f;
            useGlobalCoord_ = true;
	    fountain_xfactor_ .push_back( token5f );
	    fountain_xfactor_ .push_back( 1 );
	    fountain_xfactor_ .push_back( 1 );
	    fountain_xfactor_ .push_back( 1 );
	    fountain_xfactor_ .push_back( 1 );
	    fountain_xfactor_ .push_back( token6f );

            if (fountain_sf_ <= 0. || fountain_nz_ == 0)
                throw std::invalid_argument("Incorrect fountain superstrip definition.");
            break;

	case SuperstripType::FOUNTAINOPT:
            fountainopt_pt_ = token2f;
            fountain_nz_    = token4f;
            useGlobalCoord_ = true;

            if (fountainopt_pt_ <= 0. || fountain_nz_ == 0)
                throw std::invalid_argument("Incorrect fountain superstrip definition.");
            break;

        case SuperstripType::UNKNOWN:
        default:
            throw std::invalid_argument("Incorrect superstrip definition.");
            break;
        }

    } else {
        throw std::invalid_argument("Incorrect superstrip definition.");
    }

    // Learn about the trigger tower geometry
    if (!useGlobalCoord()) {
        towerModules_.clear();
        towerModules_.resize(16);

        const std::vector<unsigned>& moduleIds = ttmap->getTriggerTowerModules(tt);
        for (unsigned i=0; i<moduleIds.size(); ++i) {
            unsigned moduleId = moduleIds.at(i);
            unsigned lay16    = compressLayer(decodeLayer(moduleId));
            assert(lay16 < 16);

            towerModules_.at(lay16).push_back(moduleId);
        }

        // Sanity check
        for (std::vector<std::vector<unsigned> >::const_iterator it = towerModules_.begin();
             it != towerModules_.end(); ++it) {
            assert(towerModules_.size() < (1u << MODULE_NBITS));  // make sure number of bits is enough
        }

    } else {
        phiMins_.clear();
        phiMaxs_.clear();
        zMins_.clear();
        zMaxs_.clear();

        phiMins_.resize(16);
        phiMaxs_.resize(16);
        zMins_.resize(16);
        zMaxs_.resize(16);

        const std::vector<LayerBounds>& boundaries = ttmap->getTriggerTowerBoundaries(tt);
        for (unsigned i=0; i<boundaries.size(); ++i) {
            const LayerBounds& lb = boundaries.at(i);
            unsigned lay16        = compressLayer(lb.layer);
            assert(lay16 < 16);

            phiMins_.at(lay16) = lb.phiMin;
            phiMaxs_.at(lay16) = lb.phiMax;
            zMins_  .at(lay16) = lb.zMin;
            zMaxs_  .at(lay16) = lb.zMax;
        }
    }

    // Set number of possible superstrips per layer
    switch(sstype_) {
    case SuperstripType::FIXEDWIDTH:
        fixedwidth_bit_rshift1_ = most_sig_bit(fixedwidth_nstrips_);
        fixedwidth_bit_rshift2_ = most_sig_bit(MAX_NSEGMENTS / fixedwidth_nz_);

        fixedwidth_bit_lshift1_ = most_sig_bit(MAX_NSTRIPS / fixedwidth_nstrips_);
        fixedwidth_bit_lshift2_ = most_sig_bit(fixedwidth_nz_);

        nsuperstripsPerLayer_ = (1u << (MODULE_NBITS + fixedwidth_bit_lshift2_ + fixedwidth_bit_lshift1_));
        break;

    case SuperstripType::PROJECTIVE:
        projective_phiBins_.clear();
        projective_zBins_.clear();

        projective_phiBins_.resize(16);
        projective_zBins_.resize(16);

        for (unsigned i=0; i<phiMins_.size(); ++i) {
            projective_phiBins_.at(i) = M_PI / 4. / projective_nx_;  // constant
            //projective_phiBins_.at(i) = (phiMaxs_.at(i) - phiMins_.at(i)) / projective_nx_;  // variable
            projective_zBins_  .at(i) = (zMaxs_.at(i) - zMins_.at(i)) / projective_nz_;

            unsigned nx = round_to_uint((phiMaxs_.at(i) - phiMins_.at(i)) / projective_phiBins_.at(i));
            if (projective_max_nx_ < nx)
                projective_max_nx_ = nx;
        }
        nsuperstripsPerLayer_ = projective_max_nx_ * projective_nz_;
        break;

    case SuperstripType::FOUNTAIN:
        fountain_phiBins_.clear();
        fountain_zBins_.clear();

        fountain_phiBins_.resize(16);
        fountain_zBins_.resize(16);

        for (unsigned i=0; i<phiMins_.size(); ++i) {
	  fountain_phiBins_.at(i) = phiWidths_.at(i) * fountain_xfactor_[i] * fountain_sf_; 
	  fountain_zBins_  .at(i) = (zMaxs_.at(i) - zMins_.at(i)) / fountain_nz_;

	  unsigned nx = round_to_uint((phiMaxs_.at(i) - phiMins_.at(i)) / fountain_phiBins_.at(i));
	  if (fountain_max_nx_ < nx)
	    fountain_max_nx_ = nx;
        }
        //nsuperstripsPerLayer_ = fountain_max_nx_ * fountain_nz_;
        nsuperstripsPerLayer_ = 1<<11;  // 12-bit superstrip
	print();
        break;

    case SuperstripType::FOUNTAINOPT:
        fountain_phiBins_.clear();
        fountain_zBins_.clear();

        fountain_phiBins_.resize(16);
        fountain_zBins_.resize(16);

        for (unsigned i=0; i<phiMins_.size(); ++i) {
	  if (fountainopt_pt_ == 1)
	    fountain_phiBins_.at(i) = phiWidths_opt_lowpt_.at(i); 
	  else
	    fountain_phiBins_.at(i) = phiWidths_opt_highpt_.at(i); 
	  
	  fountain_zBins_  .at(i) = (zMaxs_.at(i) - zMins_.at(i)) / fountain_nz_;

	  unsigned nx = round_to_uint((phiMaxs_.at(i) - phiMins_.at(i)) / fountain_phiBins_.at(i));
	  if (fountain_max_nx_ < nx)
	    fountain_max_nx_ = nx;
        }
        //nsuperstripsPerLayer_ = fountain_max_nx_ * fountain_nz_;
        nsuperstripsPerLayer_ = 1<<11;  // 12-bit superstrip
	print();
        break;

    default:
        break;
    }

}

// _____________________________________________________________________________
unsigned SuperstripArbiter::superstripLocal(unsigned moduleId, float strip, float segment) const {
    switch (sstype_) {
    case SuperstripType::FIXEDWIDTH:
        return superstripFixedwidth(moduleId, strip, segment);
        break;

    default:
        throw std::logic_error("Incompatible superstrip type.");
        break;
    }
}

// _____________________________________________________________________________
unsigned SuperstripArbiter::superstripGlobal(unsigned moduleId, float r, float phi, float z, float ds, unsigned nDSbins, std::string DSM) const {
    switch (sstype_) {
    case SuperstripType::PROJECTIVE:
      return superstripProjective(moduleId, r, phi, z, ds);
        break;

    case SuperstripType::FOUNTAIN:
      return superstripFountain(moduleId, r, phi, z, ds, nDSbins, DSM);
        break;

    case SuperstripType::FOUNTAINOPT:
        return superstripFountain(moduleId, r, phi, z, ds);
        break;

    default:
        throw std::logic_error("Incompatible superstrip type.");
        break;
    }
}

// _____________________________________________________________________________
unsigned SuperstripArbiter::compressModuleId(unsigned moduleId) const {
    unsigned lay16    = compressLayer(decodeLayer(moduleId));
    //assert(lay16 < 16);

    std::vector<unsigned>::const_iterator found = std::lower_bound(towerModules_.at(lay16).begin(), towerModules_.at(lay16).end(), moduleId);
    if (found == towerModules_.at(lay16).end() || *found != moduleId)
        throw std::logic_error("Unexpected module ID.");

    return found - towerModules_.at(lay16).begin();
}


// _____________________________________________________________________________
unsigned SuperstripArbiter::superstripFixedwidth(unsigned moduleId, float strip, float segment) const {
    // The superstrip ID is an unsigned integer with 3 fields: moduleId|segment|strip
    // The length of the 'strip' field is controlled by fixedwidth_bitshift1
    // The length of the 'segment' field is controlled by (fixedwidth_bitshift2 - fixedwidth_bitshift1)

    unsigned moduleCode = compressModuleId(moduleId);
    strip    -= 0.25;  // round down the half int precision
    segment  -= 0.25;  // round down the half int precision

    // Use strip in the least significant field
    unsigned ss = round_to_uint(strip) >> fixedwidth_bit_rshift1_;
    //std::cout << strip+0.25 << " " << round_to_uint(strip) << " " << ss << std::endl;

    // Use segment in the 2nd most sig field
    if (isPSModule(moduleId))
        ss |= ((round_to_uint(segment) >> fixedwidth_bit_rshift2_) << fixedwidth_bit_lshift1_);
    else
        ss |= (((round_to_uint(segment)*MAX_NSEGMENTS/2) >> fixedwidth_bit_rshift2_) << fixedwidth_bit_lshift1_);
    //std::cout << segment+0.25 << " " << round_to_uint(segment) << " " << ss << std::endl;

    // Use module ID in the most sig field
    ss |= (moduleCode << (fixedwidth_bit_lshift1_ + fixedwidth_bit_lshift2_));
    //std::cout << moduleId << " " << moduleCode << " " << ss << std::endl;

    return ss;
}

// _____________________________________________________________________________
unsigned SuperstripArbiter::superstripProjective(unsigned moduleId, float r, float phi, float z, float ds) const {
    unsigned lay16    = compressLayer(decodeLayer(moduleId));

    int n_phi = projective_max_nx_;
    int n_z   = projective_nz_;

    int i_phi = std::floor((phi - phiMins_.at(lay16)) / projective_phiBins_.at(lay16));
    int i_z   = std::floor((z - zMins_.at(lay16)) / projective_zBins_.at(lay16));

    i_phi     = (i_phi < 0) ? 0 : (i_phi >= n_phi) ? (n_phi - 1) : i_phi;  // proper range
    i_z       = (i_z   < 0) ? 0 : (i_z   >= n_z  ) ? (n_z   - 1) : i_z;    // proper range

    unsigned ss = i_z * n_phi + i_phi;
    return ss;
}

// _____________________________________________________________________________
unsigned SuperstripArbiter::superstripFountain(unsigned moduleId, float r, float phi, float z, float ds, unsigned nDSbins, std::string DSM) const {
    unsigned lay16    = compressLayer(decodeLayer(moduleId));

    // For barrel, correct phi based on ds and dr
    //if (lay16 < 6)
    //    phi = phi + drCorrs_.at(lay16) * ds * (r - rMeans_.at(lay16));

    int n_phi = fountain_max_nx_;
    int n_z   = fountain_nz_;

    int i_phi = std::floor((phi - phiMins_.at(lay16)) / fountain_phiBins_.at(lay16));
    int i_z   = std::floor((z - zMins_.at(lay16)) / fountain_zBins_.at(lay16));
    //std::cout<<"i_phi: "<<i_phi<<"\t i_z: "<<i_z<<std::endl;

    i_phi     = (i_phi < 0) ? 0 : (i_phi >= n_phi) ? (n_phi - 1) : i_phi;  // proper range
    i_z       = (i_z   < 0) ? 0 : (i_z   >= n_z  ) ? (n_z   - 1) : i_z;    // proper range
    //std::cout<<"i_phi: "<<i_phi<<"\t i_z: "<<i_z<<std::endl;

    int i_ds = 0; //if not request to apply deltaS approach to a layer

    /*3bins, central_bin=7bins wide, original Sergo's approach I implemented (the outer bins are not divided)*/
    if(DSM == "ASYM" && nDSbins > 1){
      float bin_limit = ((nDSbins-1)/2.)*0.5;//0.5 = half strip                                                                                    
      if (ds < -bin_limit) i_ds = 1;
      if (ds > bin_limit) i_ds = 2;
    }


    //This approach is applied by layer
    if(DSM == "SYM" && nDSbins > 1){
      /*3bins_9binWide*/
      if(nDSbins == 3){
	if(ds <= -2.5) i_ds = 1;
	if(ds >= 2.5) i_ds = 2;
      }

    
      /*5bins_7binWide*/
      if(nDSbins == 5){
	if(ds <= -2.0 && ds >= -5.0) i_ds = 1;
	if(ds <= -5.5) i_ds = 2;
	if(ds >= 2.0 && ds <= 5.0) i_ds = 3;
	if(ds >= 5.5) i_ds = 4;      
      }


      /*7bins_5binWide*/
      if(nDSbins == 7){
	if(ds <= -1.5 && ds >= -3.5) i_ds = 1;
	if(ds <= -4.0 && ds >= -6.0) i_ds = 2;
	if(ds <= -6.5) i_ds = 3;
	if(ds >= 1.5 && ds <= 3.5) i_ds = 4;
	if(ds >= 4.0 && ds <= 6.0) i_ds = 5;
	if(ds >= 6.5) i_ds = 6;
      }
    
    
      /*9bins_3binWide config*/
      if(nDSbins == 9){
	if(ds <= -1.0 && ds >= -2) i_ds = 1;
	if(ds <= -2.5 && ds >= -3.5) i_ds = 2;
	if(ds <= -4.0 && ds >= -5.0) i_ds = 3;
	if(ds <= -5.5) i_ds = 4;
	if(ds >= 1.0 && ds <= 2) i_ds = 5;
	if(ds >= 2.5 && ds <= 3.5) i_ds = 6;
	if(ds >= 4.0 && ds <= 5.0) i_ds = 7;
	if(ds >= 5.5) i_ds = 8;
      }
    }
    
    if((nDSbins != 1 && nDSbins != 3 && nDSbins != 5 && nDSbins != 7 && nDSbins != 9) || (DSM != "ASYM" && DSM != "SYM")){
      throw std::invalid_argument( "[ERROR] Not allowed deltaS ou deltaSM! Allowed values for deltaS are [1,3,5,7,9] and for deltaSM are [ASYM,SYM]." );
    }


    unsigned ss = i_ds * n_phi + i_phi;
    //unsigned ss = ((i_z & 0x7) << 9) | (i_phi & 0x1ff);  // use magic number of 512 (= 1<<9)
    // unsigned ss =(((i_ds >> 4 ) & 0x1) << 9) | (i_phi & 0x1ff);  // use magic number of 512 (= 1<<9)
    //unsigned ss =((i_ds & 0x3) << 9) | (i_phi & 0x1ff);  // use magic number of 512 (= 1<<9)

    //std::cout<<"ds = "<<ds<<"\tss = "<<ss<<std::endl;
    return ss;
}

// _____________________________________________________________________________
void SuperstripArbiter::print() {
    switch (sstype_) {
    case SuperstripType::FIXEDWIDTH:
        std::cout << "Using fixedwidth superstrip with nstrips: " << fixedwidth_nstrips_ << ", nz: " << fixedwidth_nz_
                  << ", left shifts: " << fixedwidth_bit_lshift1_ << "," << fixedwidth_bit_lshift2_
                  << ", right shifts: " << fixedwidth_bit_rshift1_ << "," << fixedwidth_bit_rshift2_
                  << ", nsuperstrips per layer: " << nsuperstripsPerLayer_ << std::endl;
        break;

    case SuperstripType::PROJECTIVE:
        std::cout << "Using projective superstrip with nx: " << projective_nx_ << ", nz: " << projective_nz_
                  << ", phi bins: " << projective_phiBins_.at(0) << "," << projective_phiBins_.at(1) << "," << projective_phiBins_.at(2)
                  << "," << projective_phiBins_.at(3) << "," << projective_phiBins_.at(4) << "," << projective_phiBins_.at(5)
                  << ", z bins: " << projective_zBins_.at(0) << "," << projective_zBins_.at(1) << "," << projective_zBins_.at(2)
                  << "," << projective_zBins_.at(3) << "," << projective_zBins_.at(4) << "," << projective_zBins_.at(5)
                  << ", nsuperstrips per layer: " << nsuperstripsPerLayer_ << std::endl;
        break;

    case SuperstripType::FOUNTAIN:
        std::cout << "Using fountain superstrip with sf: " << fountain_sf_ << ", nz: " << fountain_nz_
                  << ", phi bins: " << fountain_phiBins_.at(0) << "," << fountain_phiBins_.at(1) << "," << fountain_phiBins_.at(2)
                  << "," << fountain_phiBins_.at(3) << "," << fountain_phiBins_.at(4) << "," << fountain_phiBins_.at(5)
                  << ", z bins: " << fountain_zBins_.at(0) << "," << fountain_zBins_.at(1) << "," << fountain_zBins_.at(2)
                  << "," << fountain_zBins_.at(3) << "," << fountain_zBins_.at(4) << "," << fountain_zBins_.at(5)
                  << ", nsuperstrips per layer: " << nsuperstripsPerLayer_ << std::endl;
        break;

    case SuperstripType::FOUNTAINOPT:
        std::cout << "Using optimized pt range: " << fountainopt_pt_ << ", nz: " << fountain_nz_
                  << ", phi bins: " << fountain_phiBins_.at(0) << "," << fountain_phiBins_.at(1) << "," << fountain_phiBins_.at(2)
                  << "," << fountain_phiBins_.at(3) << "," << fountain_phiBins_.at(4) << "," << fountain_phiBins_.at(5)
                  << ", z bins: " << fountain_zBins_.at(0) << "," << fountain_zBins_.at(1) << "," << fountain_zBins_.at(2)
                  << "," << fountain_zBins_.at(3) << "," << fountain_zBins_.at(4) << "," << fountain_zBins_.at(5)
                  << ", nsuperstrips per layer: " << nsuperstripsPerLayer_ << std::endl;
        break;

    default:
        break;
    }
}
