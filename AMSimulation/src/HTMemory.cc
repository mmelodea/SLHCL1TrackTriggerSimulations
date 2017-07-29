//
//  HTMemory.cpp
//  HTMSimulation
//
//  Created by Luciano Ristori on 10/15/15.
//  Copyright (c) 2015 Luciano Ristori. All rights reserved.
//
//
//  2017/07/29: Included in this package and modified with permission from
//              Luciano Ristori.
//

#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/HTMemory.h"

#include <cassert>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <set>
#include <random>

#include "TFile.h"


///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

void HTMemory::Hit::print() const {
    std::cout << "hit: " << indHit << " layer: " << iLayer << " x: " << x << std::endl;
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

void HTMemory::Candidate::print(const std::vector<int>& tpIdList) const {
    for (unsigned iL = 0; iL != six_layers; ++iL) {
        std::cout << "L" << iL << " ";
        unsigned nHits = hitList[iL].size();
        for (unsigned iH = 0; iH != nHits; ++iH) {
            std::cout << tpIdList.at(hitList[iL][iH].indHit) << " ";
        }
    }
    std::cout << std::endl;
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

void HTMemory::printCandidateList(const std::vector<Candidate> &list, const std::vector<int>& tpIdList) const {
    for (unsigned iC = 0; iC != list.size(); ++iC) {
        std::cout << "C" << iC << " ";
        list[iC].print(tpIdList);
    }
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

HTMemory::Cell::Cell(int _nLayers, parType _a, parType _b) {
    a = _a;
    b = _b;

    // initialize vectors of limits
    for (int i = 0; i != _nLayers; ++i) {
        xMin.push_back(0.);
        xMax.push_back(0.);
        nHits.push_back(0);
    }

    nHitLayers = 0;
    suppressed = false;
    candidate.clear();
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

// HTMemory default constructor

HTMemory::HTMemory()
    : h_a_train(0), h_b_train(0), h_n6(0), h_n5(0), h_n56(0), gridTrain(0), grid(0), gridCanvas(0)
{

}


///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

// HTMemory constructor, all parameters given as arguments

HTMemory::HTMemory(const std::string _name, int _nLayers, int _na, parType _aMin, parType _aMax, int _nb, parType _bMin, parType _bMax)
    : h_a_train(0), h_b_train(0), h_n6(0), h_n5(0), h_n56(0), gridTrain(0), grid(0), gridCanvas(0)
{

    // initialize htm

    init(_name, _nLayers, _na, _aMin, _aMax, _nb, _bMin, _bMax);

    // initialize training track distribution histograms

    //h_a_train = new TH1D((name + "h_a_train").c_str(), (name + "h_a_train").c_str(), 100, 0.,-1.);
    //h_b_train = new TH1D((name + "h_b_train").c_str(), (name + "h_b_train").c_str(), 100, 0.,-1.);
    h_a_train = new TH1D((name + "h_a_train").c_str(), (name + "h_a_train").c_str(), 100, aMin,aMax);
    h_b_train = new TH1D((name + "h_b_train").c_str(), (name + "h_b_train").c_str(), 100, bMin,bMax);

    // initialize number of clusters histograms

    h_n6 = new TH1D((name + "h_n6").c_str(), (name + "h_n6").c_str(), 100, 0.,-1.);
    h_n5 = new TH1D((name + "h_n5").c_str(), (name + "h_n5").c_str(), 100, 0.,-1.);
    h_n56 = new TH1D((name + "h_n56").c_str(), (name + "h_n56").c_str(), 100, 0.,-1.);

    // initialize histogram array histograms[na][nb][nLayers]

    // 3-vector with na x nb x nLayers histograms

    const int nHistoBins = 100;  // number of bins in histograms

    for (int ia = 0; ia != na; ++ia) {

        // 2-vector with nb x nLayers histograms
        std::vector<std::vector<TH1D*> > hvec2;
        for (int ib = 0; ib != nb; ++ib) {

            // 1-vector with nLayers histograms
            std::vector<TH1D*> hvec1;
            for (int il = 0; il != nLayers; ++il) {
                std::stringstream ss;
                ss << name << ia << "_" << ib << "_" << il;
                std::string sss = ss.str();
                TString title = TString(sss.c_str());
                TH1D* h = new TH1D(title, title, nHistoBins, 0., -1.);
                hvec1.push_back(h);
            }

            hvec2.push_back(hvec1);
        }

        histograms.push_back(hvec2);
    } // end initialize histogram array

}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

// HTMemory constructor, configuration file given as argument

HTMemory::HTMemory(const std::string _name, const std::string configFileName)
    : h_a_train(0), h_b_train(0), h_n6(0), h_n5(0), h_n56(0), gridTrain(0), grid(0), gridCanvas(0)
{

    if (configFileName == "")  return;  // do nothing

    int _nLayers;
    int _na;
    parType _aMin;
    parType _aMax;
    int _nb;
    parType _bMin;
    parType _bMax;

    // open file for input of htm configuration data

    std::ifstream configFile;
    configFile.open(configFileName);
    if (!configFile) {
        std::cout << "error opening " << configFileName << std::endl;
        throw std::runtime_error("error opening configFile");
        return;
    }

    // read the first line with parameter array configuration

    configFile
    >> _nLayers
    >> _na
    >> _aMin
    >> _aMax
    >> _nb
    >> _bMin
    >> _bMax;

    // initialize htm

    init(_name, _nLayers, _na, _aMin, _aMax, _nb, _bMin, _bMax);

    // loop on all cells and load xMin and xMax

    int _ia, _ib, _il;
    parType _a, _b;
    xType _xMin, _xMax;

    for (int ia = 0; ia != na; ++ia) {
        for (int ib = 0; ib != nb; ++ib) {
            Cell* cellPointer = &cells[ia][ib];
            for (int il = 0; il != nLayers; ++il) {
                configFile
                >> _ia
                >> _ib
                >> _a
                >> _b
                >> _il
                >> _xMin
                >> _xMax;

                assert(ia == _ia);
                assert(ib == _ib);
                assert(il == _il);

                cellPointer->a = _a;
                cellPointer->b = _b;
                cellPointer->xMin[il] = _xMin;
                cellPointer->xMax[il] = _xMax;
            }
        }
    }  // end loop on all cells

    configFile.close();
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

HTMemory::~HTMemory() {
    if (h_a_train)  delete h_a_train;
    if (h_b_train)  delete h_b_train;
    if (h_n6)       delete h_n6;
    if (h_n5)       delete h_n5;
    if (h_n56)      delete h_n56;
    if (gridTrain)  delete gridTrain;
    if (grid)       delete grid;
    if (gridCanvas) delete gridCanvas;

    for (auto it1 = histograms.cbegin(); it1 != histograms.cend(); ++it1)
        for (auto it2 = it1->cbegin(); it2 != it1->cend(); ++it2)
            for (auto it3 = it2->cbegin(); it3 != it2->cend(); ++it3)
                if (*it3)  delete *it3;
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

void HTMemory::init(const std::string _name, int _nLayers, int _na, parType _aMin, parType _aMax, int _nb, parType _bMin, parType _bMax) {

    // sanity check

    assert(_nLayers <= 10);
    assert(_na <= 1000);
    assert(_nb <= 1000);

    // copy parameters to local

    name = _name;
    nLayers = _nLayers;
    na = _na;
    aMin = _aMin;
    aMax = _aMax;
    nb = _nb;
    bMin = _bMin;
    bMax = _bMax;

    // initialize training track distribution histograms

    h_a_train = new TH1D((name + "h_a_train").c_str(), (name + "h_a_train").c_str(), 100, 0.,-1.);
    h_b_train = new TH1D((name + "h_b_train").c_str(), (name + "h_b_train").c_str(), 100, 0.,-1.);

    // initialize number of clusters histograms

    h_n6 = new TH1D((name + "h_n6").c_str(), (name + "h_n6").c_str(), 100, 0.,-1.);
    h_n5 = new TH1D((name + "h_n5").c_str(), (name + "h_n5").c_str(), 100, 0.,-1.);
    h_n56 = new TH1D((name + "h_n56").c_str(), (name + "h_n56").c_str(), 100, 0.,-1.);

    // initialize parameter grid

    grid = new TH2D((name + "grid").c_str(), (name + "grid").c_str(), na, aMin, aMax, nb, bMin, bMax);
    gridTrain = new TH2D((name + "gridTrain").c_str(), (name + "gridTrain").c_str(), na, aMin, aMax, nb, bMin, bMax);

    // initialize cell array cells[na][nb]

    for (int ia = 0; ia != na; ++ia) {
        std::vector<Cell> tempCells;
        for (int ib = 0; ib != nb; ++ib) {
            parType a = grid->GetXaxis()->GetBinCenter(ia+1);
            parType b = grid->GetYaxis()->GetBinCenter(ib+1);
            Cell tempCell(nLayers, a, b);
            tempCells.push_back(tempCell);
        }
        cells.push_back(tempCells);
    }

    // canvas position and size

    int x1 = 50;
    int y1 = 50;
    int width = 800;
    int height = 1000*(double)nb/(double)na;

    gridCanvas = new TCanvas((name + "gridCanvas").c_str(), (name + "gridCanvas").c_str(), x1, y1, width, height);

} // end HTMemory::init

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

void HTMemory::writeHists(const std::string _fileName, bool _writeHistogramArray) const {

    TString fileName = _fileName.c_str();
    TFile* histFile = TFile::Open(fileName,"RECREATE");  // histogram file

    // create directory for histograms related to training

    TDirectory *trainDir = histFile->mkdir("training");
    trainDir->cd();
    gridTrain->SetStats(0);
    gridTrain->Write(); // write 2Dim training track parameter
    h_a_train->Write(); // write 1Dim training track parameter
    h_b_train->Write(); // write 1Dim training track parameter

    // create directory for histograms related to simulation

    TDirectory *simDir = histFile->mkdir("simulation");
    simDir->cd();

    // cluster multiplicity

    h_n6->Write();
    h_n5->Write();
    h_n56->Write();


    if (!_writeHistogramArray) {
        histFile->Write();
        return;
    }


    // big histogram array

    TDirectory *gridDir = trainDir->mkdir("histogram grid");
    gridDir->cd();

    // loop on histogram array

    for (int ia = 0; ia != na; ++ia) {
        for (int ib = 0; ib != nb; ++ib) {

            // create a subdirectory named h_ia_ib

            std::stringstream ss;
            ss << "h_" << ia << "_" << ib;
            std::string sss = ss.str();
            TString title = sss.c_str();
            TDirectory *innerDir = gridDir->mkdir(title);
            innerDir->cd();

            // write all histograms for that (ia,ib) bin in the grid

            for (int il = 0; il != nLayers; ++il) {
                histograms[ia][ib][il]->Write();
            }
        }
    }  // end loop on histogram array

    histFile->Write();

} // end HTMemory::writeHists

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

void HTMemory::train(parType a, parType b, const std::vector<Hit>& hitList) {

    gridTrain->Fill(a,b);

    h_a_train->Fill(a);
    h_b_train->Fill(b);

    int ia = grid->GetXaxis()->FindFixBin(a) - 1;
    int ib = grid->GetYaxis()->FindFixBin(b) - 1;

    ia = std::max(ia,0);
    ia = std::min(ia,na-1);
    ib = std::max(ib,0);
    ib = std::min(ib,nb-1);

    assert(ia >= 0 && ia <= (na -1));
    assert(ib >= 0 && ib <= (nb -1));

    for (long int ih = 0; ih != (long int)hitList.size(); ++ ih) {
        histograms[ia][ib][hitList[ih].iLayer]->Fill(hitList[ih].x);
    }

} // end HTMemory::train

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

void HTMemory::setCuts(double coverage) {

    // loop on all cells and find xMin and xMax to satisfy coverage

    for (int ia = 0; ia != na; ++ia) {
        for (int ib = 0; ib != nb; ++ib) {
            for (int il = 0; il != nLayers; ++il) {
                TH1D* h = histograms[ia][ib][il];
                int nHistoBins = h->GetNbinsX();
                long int nEntries = h->GetEntries();
                if (nEntries < 10) continue;// skip layer if not enough statistics
                long int nLowTail = 0.5*nEntries*(1. - coverage);
                long int nHighTail = nEntries - nLowTail;
                int lowBin = 0;
                int highBin = -1;
                long int integral = 0;
                for (int iBin = 0; iBin != nHistoBins + 1; ++ iBin) {
                    integral += h->GetBinContent(iBin);
                    if (integral <= nLowTail) lowBin = iBin;
                    highBin = iBin;
                    if (integral >= nHighTail) break;
                }
                //assert(lowBin != -1 && highBin != -1);
                //assert(lowBin > 0);
                assert(highBin <= nHistoBins);

                //if (lowBin < 0) std::cout << "*********** " << ia << " " << ib << " " << il << " " << lowBin << std::endl;

                cells[ia][ib].xMin[il] = h->GetXaxis()->GetBinUpEdge(lowBin);
                cells[ia][ib].xMax[il] = h->GetXaxis()->GetBinUpEdge(highBin);
            }
        }
    }

} // end HTMemory::setCuts

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

void HTMemory::writeCuts(const std::string cutFileName) const {

    // open file for output of cut data

    std::ofstream cutFile;
    cutFile.open(cutFileName);
    if (!cutFile) {
        std::cout << "error opening " << cutFileName << std::endl;
        throw std::runtime_error("error opening cutFile");
        return;
    }

    // loop on all cells and print xMin and xMax

    for (int ia = 0; ia != na; ++ia) {
        for (int ib = 0; ib != nb; ++ib) {
            const Cell* cellPointer = &cells[ia][ib];
            cutFile << "--- ia " << ia << " ib " << ib << " a " << cellPointer->a << " b " << cellPointer->b << std::endl;
            for (int il = 0; il != nLayers; ++il) {
                cutFile
                << "layer " << il
                << " xmin " << cellPointer->xMin[il]
                << " xmax " << cellPointer->xMax[il]
                << " delta " << cellPointer->xMax[il] - cellPointer->xMin[il]
                << std::endl;
            }
        }
    }  // end loop on all cells

    cutFile.close();

} // end HTMemory::writeCuts

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

void HTMemory::writeConfig(const std::string configFileName) const {

    // open file for output of htm configuration data

    std::ofstream configFile;
    configFile.open(configFileName);
    if (!configFile) {
        std::cout << "error opening " << configFileName << std::endl;
        throw std::runtime_error("error opening configFile");
        return;
    }

    configFile
    << " " << nLayers
    << " " << na
    << " " << aMin
    << " " << aMax
    << " " << nb
    << " " << bMin
    << " " << bMax
    << std::endl;

    // loop on all cells and print xMin and xMax

    for (int ia = 0; ia != na; ++ia) {
        for (int ib = 0; ib != nb; ++ib) {
            const Cell* cellPointer = &cells[ia][ib];
            for (int il = 0; il != nLayers; ++il) {
                configFile
                << " " << ia
                << " " << ib
                << " " << cellPointer->a
                << " " << cellPointer->b
                << " " << il
                << " " << cellPointer->xMin[il]
                << " " << cellPointer->xMax[il]
                << std::endl;
            }
        }
    }  // end loop on all cells

    configFile.close();

} // end HTMemory::writeConfig

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

void HTMemory::input(const std::vector<Hit>& hitList) {

    // clear grid

    for (int ia = 0; ia != na; ++ia) {
        for (int ib = 0; ib != nb; ++ib) {
            Cell* cellPointer = &cells[ia][ib];
            cellPointer->nHitLayers = 0;
            cellPointer->suppressed = false;
            cellPointer->candidate.clear();
            for (int il = 0; il != nLayers; ++il) {
                cellPointer->nHits[il] = 0;
            }
        }
    }  // end clear grid

    // process hitList

    for (long int iHit = 0; iHit != (long int)hitList.size(); ++iHit) { // loop on all hits
        xType x = hitList[iHit].x;
        int iLayer = hitList[iHit].iLayer;

        for (int ia = 0; ia != na; ++ia) { // loop on all cells
            for (int ib = 0; ib != nb; ++ib) {
                Cell* cellPointer = &cells[ia][ib];
                if (cellPointer->xMin[iLayer] <= x && x <= cellPointer->xMax[iLayer]) { // test x within cuts
                    cellPointer->candidate.hitList[iLayer].push_back(hitList[iHit]); // add hit to candidate
                    cellPointer->candidate.nTotalHits++;// increment total number of hits
                    if (cellPointer->nHits[iLayer]==0) cellPointer->nHitLayers++; // count hit layers
                    if (cellPointer->nHits[iLayer]==0) cellPointer->candidate.nHitLayers++; // count hit layers
                    cellPointer->nHits[iLayer]++; // increment hit counter per layer
                }
            }
        }  // end loop on all cells
    }  // end loop on all hits

} // end HTMemory::input

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

void HTMemory::drawGrid(const std::string pictName) const {

    std::string prefix = "plots/"; // folder where the plots will be stored

    // fill grid with number of layers hit
    // for unsuppressed cells

    grid->Reset();

    for (int ia = 0; ia != na; ++ia) {
        for (int ib = 0; ib != nb; ++ib) {
            if (!cells[ia][ib].suppressed) {
                unsigned content = cells[ia][ib].nHitLayers;
                grid->SetBinContent(ia+1, ib+1, (double) content);
            }
        }
    }

    gridCanvas->cd();
    grid->SetStats(0);
    //grid->SetMarkerSize(0.6);
    grid->SetMarkerSize(1.0);
    grid->SetTitle(pictName.c_str());
    grid->DrawCopy("TEXT COL");
    std::string fullFileName = prefix + pictName + ".pdf";
    gridCanvas->SaveAs(fullFileName.c_str());

} // end HTMemory::drawGrid

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

void HTMemory::cluster() {

    for (int ia = 0; ia != na; ++ia) { // loop on all cells
        for (int ib = 0; ib != nb; ++ib) {

            //std::cout << "array cell " << ia << " " << ib << std::endl;

            Cell* c = &cells[ia][ib]; // pointer to center cell

            if (c->nHitLayers < 5) { // suppress this cell and do nothing
                c->suppressed = true;
                continue;
            }

            // skip border cells

            if (ia == 0) continue;
            if (ia == (na - 1)) continue;
            if (ib == 0) continue;
            if (ib == (nb - 1)) continue;


            // treat 8 adjacent cells

            Cell* adj; // pointer to adjacent cell

            adj = &cells[ia+1][ib+1]; // upper right
            if (adj->nHitLayers <= c->nHitLayers) adj->suppressed = true;

            adj = &cells[ia+1][ib]; // center right
            if (adj->nHitLayers <= c->nHitLayers) adj->suppressed = true;

            adj = &cells[ia+1][ib-1]; // lower right
            if (adj->nHitLayers <= c->nHitLayers) adj->suppressed = true;

            adj = &cells[ia][ib-1]; // lower center
            if (adj->nHitLayers <= c->nHitLayers) adj->suppressed = true;

            adj = &cells[ia-1][ib-1]; // lower left
            if (adj->nHitLayers < c->nHitLayers) adj->suppressed = true;

            adj = &cells[ia-1][ib]; // center left
            if (adj->nHitLayers < c->nHitLayers) adj->suppressed = true;

            adj = &cells[ia-1][ib+1]; // upper left
            if (adj->nHitLayers < c->nHitLayers) adj->suppressed = true;

            adj = &cells[ia][ib+1]; // upper center
            if (adj->nHitLayers < c->nHitLayers) adj->suppressed = true;
        }
    }

    // count number of cells remaining after clustering

    int n6 = 0; // number of cells with 6 layers
    int n5 = 0; // number of cells with 5 layers

    for (int ia = 0; ia != na; ++ia) {
        for (int ib = 0; ib != nb; ++ib) {

            //std::cout << "counting "<< ia << " " << ib << std::endl;

            if (!cells[ia][ib].suppressed) {
                unsigned content = cells[ia][ib].nHitLayers;
                if (content == 6) ++n6;
                if (content == 5) ++n5;
            }
        }
    }

    //std::cout << "Filling " << n6 << " " << n5 << std::endl;

    h_n6->Fill(n6); // fill histograms with number of cells
    h_n5->Fill(n5);
    h_n56->Fill(n5+n6);

} // end HTMemory::cluster

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

HTMemory::Candidate HTMemory::getBestCandidate() const {

    // find best cell in grid and return Candidate
    // best cell is defined as max number of layers hit (5 or 6)

    std::vector<const Cell*> cells5, cells6;

    // create two lists of all cells with 5 and 6 layers

    for (int ia = 0; ia != na; ++ia) {
        for (int ib = 0; ib != nb; ++ib) {
            const Cell* c = &cells[ia][ib]; // pointer to cell
            if (c->nHitLayers == 5) cells5.push_back(c);
            if (c->nHitLayers == 6) cells6.push_back(c);
        }
    }

    // create 6 sets of stubs, one per layer (to avoid repetitions).
    // fill the sets with stubs from all the selected cells (with 5 or 6 layers hit)

    std::set<Hit, SortHits> stubs[six_layers];

    // do cells with 5 layers

    for (unsigned iCell = 0; iCell != cells5.size(); ++iCell) {
        const Candidate& thisCandidate = cells5[iCell]->candidate;
        for (unsigned iLayer = 0; iLayer != six_layers; ++iLayer) {
            for (unsigned iHit = 0; iHit != thisCandidate.hitList[iLayer].size(); ++iHit) {
                stubs[iLayer].insert(thisCandidate.hitList[iLayer][iHit]);
            }
        }
    }

    // do cells with 6 layers

    for (unsigned iCell = 0; iCell != cells6.size(); ++iCell) {
        const Candidate& thisCandidate = cells6[iCell]->candidate;
        for (unsigned iLayer = 0; iLayer != six_layers; ++iLayer) {
            for (unsigned iHit = 0; iHit != thisCandidate.hitList[iLayer].size(); ++iHit) {
                stubs[iLayer].insert(thisCandidate.hitList[iLayer][iHit]);
            }
        }
    }

    // return one candidate with all hits compatible with the best cells (5 or 6 layers)

    Candidate returnCand;
    returnCand.clear();

    // copy sets into stub lists for the returned candidate

    for (unsigned iL = 0; iL != six_layers; ++iL) {
        if (stubs[iL].size()) ++returnCand.nHitLayers;
        for (auto iH = stubs[iL].cbegin(); iH != stubs[iL].cend(); iH++) {
            ++returnCand.nTotalHits;
            returnCand.hitList[iL].push_back(*iH);
        }
    }

    return returnCand;

} // end HTMemory::getBestCandidate

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

// Extra function for use in amsim

HTMemory::Candidate HTMemory::getBestCandidate_amsim(
    const std::string htmconf, const std::vector<std::vector<unsigned> >& road_stubRefs, const std::vector<float>& stubs_z,
    std::vector<std::vector<unsigned> >& new_road_stubRefs
) {

    static std::default_random_engine generator;

    Candidate tempCand;
    tempCand.clear();

    // prepare list of hits for HTM input with r-z coordinates

    std::vector<Hit> hitList;

    for (unsigned iL = 0; iL != road_stubRefs.size(); ++iL) {
        unsigned nS = road_stubRefs.at(iL).size();
        for (unsigned iS = 0; iS != nS; ++iS) {
            unsigned istub = road_stubRefs.at(iL).at(iS);
            float stub_z = stubs_z.at(istub);

            Hit tempHit;
            tempHit.iLayer = iL; // layer number [0,5]
            tempHit.indHit = istub;
            tempHit.x = stub_z; // pick z cordinate
            hitList.push_back(tempHit);

            if (tempCand.hitList[iL].size()==0) ++tempCand.nHitLayers;
            tempCand.hitList[iL].push_back(tempHit);
            tempCand.nTotalHits++;
        }
    }

    // select filtering mode

    TString htmconf_tstring   = htmconf.c_str();
    bool kHTFilter      = !cells.empty();
    bool kHTCandidate   = htmconf_tstring.Contains("C");
    bool kHTCombination = htmconf_tstring.Contains("b");

    if (kHTFilter && !kHTCandidate && kHTCombination) {
        std::cout << "[WARNING] kHTCandidate is set to false but kHTCombination is set to true. Ignoring kHTCombination." << std::endl;
        kHTCombination = false;
    }

    if (kHTFilter) {

        // input hit list to htm

        input(hitList);

        // get best candidate

        Candidate bestCandidate = getBestCandidate();

        if (bestCandidate.nTotalHits > 0) {
            if (!kHTCandidate && !kHTCombination) {
                // output tempCand
                //tempCand = tempCand;
            } else if (kHTCandidate && !kHTCombination) {
                // output bestCandidate
                tempCand = bestCandidate;
            } else if (kHTCandidate &&  kHTCombination) {
                // output modified bestCandidate
                unsigned tempNHits = 0;
                for (unsigned iL=0; iL != six_layers; ++iL) {
                    std::vector<Hit> tempHitList;
                    unsigned nHits = bestCandidate.hitList[iL].size();
                    if(nHits) ++tempNHits;
                    if(nHits > 1) {
                        std::uniform_int_distribution<int> flat(0,nHits-1);
                        unsigned iHrand = flat(generator);
                        tempHitList.push_back(bestCandidate.hitList[iL][iHrand]);
                        bestCandidate.hitList[iL] = tempHitList;
                    }
                }
                bestCandidate.nTotalHits = tempNHits;
                tempCand = bestCandidate;
            }

        } else { // if bestCandidate.nTotalHits == 0
            // output bestCandidate
            tempCand = bestCandidate;
        }

    } else {  // if !kHTFilter
        // output tempCand
        //tempCand = tempCand;
    }

    // create new stubRefs

    new_road_stubRefs.clear();

    for (unsigned iL=0; iL != six_layers; ++iL) {
        std::vector<unsigned> temp;
        unsigned nS = tempCand.hitList[iL].size();
        for (unsigned iS = 0; iS != nS; ++iS) {
            unsigned istub = tempCand.hitList[iL][iS].indHit;
            temp.push_back(istub);
        }
        new_road_stubRefs.push_back(temp);
    }

    assert(road_stubRefs.size() == new_road_stubRefs.size());
    //if (!kHTFilter)  assert(road_stubRefs == new_road_stubRefs);
    //if (kHTFilter && !kHTCandidate && !kHTCombination && tempCand.nTotalHits > 0)  assert(road_stubRefs == new_road_stubRefs);

    return tempCand;
}
