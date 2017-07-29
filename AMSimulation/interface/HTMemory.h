//
//  HTMemory.h
//  HTMSimulation
//
//  Created by Luciano Ristori on 10/15/15.
//  Copyright (c) 2015 Luciano Ristori. All rights reserved.
//
//
//  2017/07/29: Included in this package and modified with permission from
//              Luciano Ristori.
//

#ifndef __HTMSimulation__HTMemory__
#define __HTMSimulation__HTMemory__

#include <string>
#include <vector>

#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"


class HTMemory {

public:

    typedef double parType;
    typedef double xType;

    static const unsigned int six_layers = 6;

    ///////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////

    //
    // class Hit
    //

    class Hit {
    public:
        int indHit;
        int iLayer;
        xType x;

        void print() const;
    };

    struct SortHits {
        bool operator()(const Hit &a, const Hit &b) const {
            return a.indHit < b.indHit;
        }
    };

    ///////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////

    //
    // class Candidate
    //

    class Candidate {
    public:
        unsigned nTotalHits; // total number of hits
        unsigned nHitLayers; // number of layers with at least one hit
        std::vector<Hit> hitList[six_layers];  // holds a list of Hits for each detector layer

        void clear(){
            nTotalHits = 0;
            nHitLayers = 0;
            for (unsigned iL = 0; iL != six_layers; ++iL) {
                hitList[iL].clear();
            }
        }

        void print(const std::vector<int>& tpIdList) const;
    };

    void printCandidateList(const std::vector<Candidate> &list, const std::vector<int>& tpIdList) const;

    ///////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////

    //
    // class Cell
    //

    class Cell {
    public:
        parType a;
        parType b;
        //std::vector<xType> xCenter; // cell center (where ideal track goes through)
        std::vector<xType> xMin; // low cut per layer
        std::vector<xType> xMax; // high cut per layer
        std::vector<int> nHits; // number of hits in each layer
        unsigned nHitLayers; // number of layers with at least one hit
        bool suppressed; // cell is suppressed by clustering

        Candidate candidate; // track candidate

        Cell(int _nLayers, parType _a, parType _b); // constructor
    };

    ///////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////

    //
    // class HTMemory
    //

    // ATTRIBUTES

public:

    std::string name;
    int nLayers;
    int na;
    parType aMin;
    parType aMax;
    int nb;
    parType bMin;
    parType bMax;

    // histograms

    TH1D* h_a_train; // training track distribution
    TH1D* h_b_train;

    TH1D* h_n6; // number of clusters - 6 layers
    TH1D* h_n5; // number of clusters - 5 layers
    TH1D* h_n56; // number of clusters - 5 plus 6 layers

    TH2D* gridTrain; // training track parameter grid
    TH2D* grid; // track parameter grid for HT

    TCanvas* gridCanvas; // to plot grid (see drawGrid())

    ///////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////
    //    histogram array
    //    na x nb x nLayers

    std::vector<std::vector<std::vector<TH1D*> > > histograms;

    ///////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////
    //    cell array
    //    na x nb

    std::vector<std::vector<Cell> > cells;

    ///////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////

    // METHODS

private:

    void init(const std::string _name, int _nLayers, int _na, parType _aMin, parType _aMax, int _nb, parType _bMin, parType _bMax);

public:

    // HTMemory default constructor
    HTMemory();

    // HTMemory constructor, all parameters given as arguments
    HTMemory(const std::string _name, int _nLayers, int _na, parType _aMin, parType _aMax, int _nb, parType _bMin, parType _bMax);

    // HTMemory constructor, configuration file given as argument
    HTMemory(const std::string _name, const std::string fileName);

    ~HTMemory();

    void train(parType a, parType b, const std::vector<Hit>& hitList);

    void setCuts(double coverage);

    void input(const std::vector<Hit>& hitList);

    void cluster();

    //void output(std::vector<Candidate>& candList);

    void writeHists(const std::string fileName, bool writeHistogramArray) const;

    void writeCuts(const std::string fileName) const;

    void writeConfig(const std::string filename) const;

    void drawGrid(const std::string fileName) const;

    Candidate getBestCandidate() const;

    // Extra function for use in amsim
    Candidate getBestCandidate_amsim(
        // input
        const std::string htmconf, const std::vector<std::vector<unsigned> >& road_stubRefs, const std::vector<float>& stubs_z,
        // output
        std::vector<std::vector<unsigned> >& new_road_stubRefs
    );
};


#endif /* defined(__HTMSimulation__HTMemory__) */
