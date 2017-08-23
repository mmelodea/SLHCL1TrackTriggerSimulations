#ifndef AMSimulation_CombinationBuilderFactory_h_
#define AMSimulation_CombinationBuilderFactory_h_

#include <vector>
#include <memory>
//#include <iostream>
#include "LinearizedTrackFit/LinearizedTrackFit/interface/SimpleCombinationBuilder.h"
#include "LinearizedTrackFit/LinearizedTrackFit/interface/AdvancedCombinationBuilder.h"

namespace slhcl1tt {

class CombinationBuilderFactory
{
 public:
  // Constructor
  CombinationBuilderFactory(const bool advancedCombinationBuilder) : verbose_(true)
  {
    // bool advancedCombinationBuilder = false;
    bool fillEmptyStubs = true;
    if (advancedCombinationBuilder){
      std::cout<<"Using advanced combination builder..."<<std::endl;
      combinationBuilder_ = std::make_shared<AdvancedCombinationBuilder>(fillEmptyStubs);
    }
    else{
      std::cout<<"Using simple combination builder..."<<std::endl;
      combinationBuilder_ = std::make_shared<SimpleCombinationBuilder>(fillEmptyStubs);
    }
  }

  // Destructor
  ~CombinationBuilderFactory() {}

  // Enum
  enum Flag { BAD=999999999 };

  // Return all possible combinations
  template<typename T> std::vector<std::vector<T> > combine(const std::vector<std::vector<T> >& stubRefs)
  {
    std::vector<T> combination;
    std::vector<std::vector<T> > combinations;

    // Build a Road to pass to the CB
    Road road(stubRefs);
    // std::cout << "filled road" << std::endl;
    combinationBuilder_->initialize(road);
    // std::cout << "initialized CB" << std::endl;
    // Could use an R-value reference in the CB for this?
    // combinationBuilder_.initialize(Road(m_hits));

    // std::cout << "total combinations for this road = " << combinationBuilder_->totalCombinations() << std::endl;

    for (int comb = 0; comb < combinationBuilder_->totalCombinations(); ++comb) {
      const StubsCombination& stubsCombination(combinationBuilder_->nextCombination());
      // std::cout << "loaded combination number " << comb << std::endl;
      combination.clear();
      for (const auto& s : stubsCombination) {
	// std::cout << "filling stubRef = " << s.stubRef() << std::endl;
	combination.push_back(s.stubRef());
	// if (s.size())
	//   combination.push_back(s.at(i));
	// else // empty layer
	//   combination.push_back(CombinationBuilderFactory::BAD);
      }
      combinations.push_back(combination);
    }
    return combinations;
  }

  // Debug
  void print();

 private:
  // Member data
  int verbose_;
  std::shared_ptr<CombinationBuilderBase> combinationBuilder_;
};

}

#endif
