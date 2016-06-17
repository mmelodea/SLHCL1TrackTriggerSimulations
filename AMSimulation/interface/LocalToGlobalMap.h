#ifndef AMSimulation_LocalToGlobalMap_h_
#define AMSimulation_LocalToGlobalMap_h_

#include <map>
#include <string>
#include <vector>
#include <stdint.h>
#include "TString.h"


namespace slhcl1tt {

// Define local-to-global coefficients
struct LocalToGlobal {
  public:
    float x_phi0;
    float x_phi;
    float x_z0;
    float x_z;
    float x_r0;
    float x_r;

    void set(const std::vector<float>& values) {
        x_phi0 = values.at(0);
        x_phi  = values.at(1);
        x_z0   = values.at(2);
        x_z    = values.at(3);
        x_r0   = values.at(4);
        x_r    = values.at(5);
    }
};

struct LocalToGlobalInt {
  public:
    int64_t i_phi0;
    int64_t i_phi;
    int64_t i_z0;
    int64_t i_z;
    int64_t i_r0;
    int64_t i_r;
};

// Define local-to-global conversion
class LocalToGlobalMap {

  public:
    // Constructor
    LocalToGlobalMap();

    // Destructor
    ~LocalToGlobalMap() {}

    // Functions
    // Read local-to-global conversion csv file
    void read(TString datadir);

    // Read local-to-global conversion csv file
    void readLocalToGlobalMap(TString csvfile);

    void convert(const unsigned moduleId, const float strip, const float segment,
                 float& conv_r, float& conv_phi, float& conv_z, LocalToGlobal& conv_l2g);

    void convertInt(const unsigned moduleId, const float strip, const float segment, const unsigned tt, const LocalToGlobal& conv_l2g,
                    int64_t& conv_r, int64_t& conv_phi, int64_t& conv_z, LocalToGlobalInt& conv_l2g_int);

    // Debug
    void print();

  private:
    std::map<std::pair<unsigned, unsigned>, LocalToGlobal> l2gmap_;
};


}  // namespace slhcl1tt

#endif
