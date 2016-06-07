#ifndef AMSimulation_Local2GlobalMap_h_
#define AMSimulation_Local2GlobalMap_h_

#include <map>
#include <string>
#include <vector>
#include "TString.h"


namespace slhcl1tt {

// Define local-to-global coefficients
struct Local2Global {
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

// Define local-to-global conversion
class Local2GlobalMap {

  public:
    // Constructor
    Local2GlobalMap();

    // Destructor
    ~Local2GlobalMap() {}

    // Functions
    // Read local-to-global conversion csv file
    void read(TString datadir);

    // Read local-to-global conversion csv file
    void readLocal2GlobalMap(TString csvfile);

    void convert(const unsigned moduleId, const float strip, const float segment,
                 float& conv_r, float& conv_phi, float& conv_z, Local2Global& conv_l2g);

    // Debug
    void print();

  private:
    std::map<std::pair<unsigned, unsigned>, Local2Global> l2gmap_;

};


}  // namespace slhcl1tt

#endif
