#!/usr/bin/env python

import os
import sys
from math import sqrt, pi
from itertools import izip

from ROOT import TTree, TFile, TH1F, TLegend, gROOT, gSystem, gPad

# ______________________________________________________________________________
# Options

class Options:
    def __init__(self):
        pass

options = Options()

options.nevents = 9999

options.event_list = []
#options.event_list = [7, 43, 59, 152, 202, 221, 235, 375, 439,]

options.outdir = "figures/"
options.donotdelete = []

tlegend = TLegend(0.70,0.74,0.96,0.94)
tlegend.SetFillStyle(0)
tlegend.SetLineColor(0)
tlegend.SetShadowColor(0)
tlegend.SetBorderSize(0)


# ______________________________________________________________________________
# Speed up TTree loading

def set_branch_status(tree):
    tree.SetBranchStatus("*", 0)
    tree.SetBranchStatus("AMTTRoads_patternRef"     , 1)
    #tree.SetBranchStatus("AMTTRoads_patternInvPt"   , 1)
    #tree.SetBranchStatus("AMTTRoads_superstripIds"  , 1)
    #tree.SetBranchStatus("AMTTRoads_stubRefs"       , 1)
    #tree.SetBranchStatus("AMTTTracks_pt"            , 1)
    #tree.SetBranchStatus("AMTTTracks_eta"           , 1)
    #tree.SetBranchStatus("AMTTTracks_phi0"          , 1)
    #tree.SetBranchStatus("AMTTTracks_z0"            , 1)
    #tree.SetBranchStatus("AMTTTracks_invPt"         , 1)
    #tree.SetBranchStatus("AMTTTracks_chi2"          , 1)
    #tree.SetBranchStatus("AMTTTracks_ndof"          , 1)
    tree.SetBranchStatus("AMTTTracks_patternRef"    , 1)
    #tree.SetBranchStatus("AMTTTracks_roadRef"       , 1)
    #tree.SetBranchStatus("AMTTTracks_combRef"       , 1)
    #tree.SetBranchStatus("AMTTTracks_stubRefs"      , 1)
    #tree.SetBranchStatus("AMTTTracks_tpId"          , 1)
    tree.SetBranchStatus("AMTTTracks_synTpId"       , 1)
    #tree.SetBranchStatus("trkParts_primary"         , 1)
    #tree.SetBranchStatus("trkParts_pt"              , 1)
    #tree.SetBranchStatus("trkParts_eta"             , 1)
    #tree.SetBranchStatus("trkParts_phi"             , 1)
    #tree.SetBranchStatus("trkParts_vz"              , 1)
    #tree.SetBranchStatus("trkParts_charge"          , 1)

def reset_branch_status(tree):
    tree.SetBranchStatus("*", 1)


# ______________________________________________________________________________
# Counters

def count_nroads(evt):
    return evt.AMTTRoads_patternRef.size()

def count_ntracks(evt):
    return len(filter(lambda x: x != -1 and x != -2, evt.AMTTTracks_synTpId))

# ______________________________________________________________________________
# Drawing

def style(h):
    h.SetLineWidth(2)

def draw(h_nroads1, h_nroads2, h_dntracks):
    # nroads
    if True:
        h1 = h_nroads1
        h2 = h_nroads2
        style(h1)
        style(h2)

        h1.SetMaximum(h2.GetMaximum() * 1.4)
        h1.Draw()

        h2.SetLineColor(632)  # kRed
        h2.Draw("same")

        tlegend.AddEntry(h1, "before (#mu=%.1f)" % h1.GetMean(), "l")
        tlegend.AddEntry(h2, "after (#mu=%.1f)" % h2.GetMean(), "l")
        tlegend.Draw()
        gPad.Print(options.outdir+"h_nroads.png")

    # ntracks
    if True:
        h1 = h_dntracks
        style(h1)

        h1.SetMaximum(h1.GetMaximum() * 1.4)
        h1.SetMarkerSize(2)
        h1.Draw("hist text")
        gPad.Print(options.outdir+"h_dntracks.png")

    options.donotdelete.append([h_nroads1, h_nroads2, h_dntracks])
    return

# ______________________________________________________________________________
# Main functions

def init():

    base = os.getenv("CMSSW_BASE")
    gROOT.LoadMacro(base+"/src/SLHCL1TrackTriggerSimulations/AMSimulationIO/src/AMSimulationIOLinkDef.h")
    gROOT.LoadMacro(base+"/src/SLHCL1TrackTriggerSimulations/AMSimulation/res1/tdrstyle.C")
    gROOT.ProcessLine("setTDRStyle()")

    if not os.path.exists(options.outdir):
        os.makedirs(options.outdir)

def main():

    # Get input file name
    if len(sys.argv) != 3:
        print("Usage: %s %s %s" % (sys.argv[0], "tracks1.root", "tracks2.root"))
        raise Exception("Cannot parse arguments")

    fname1 = sys.argv[1]
    fname2 = sys.argv[2]

    tfile1 = TFile.Open(fname1)
    ttree1 = tfile1.Get("ntupler/tree")
    if not ttree1:
        raise Exception("Cannot retrieve the first tree")

    tfile2 = TFile.Open(fname2)
    ttree2 = tfile2.Get("ntupler/tree")
    if not ttree2:
        raise Exception("Cannot retrieve the second tree")

    # Speed up TTree loading
    set_branch_status(ttree1)
    set_branch_status(ttree2)

    # Book histograms
    h_nroads1 = TH1F("nroads1", "; # roads/tower/BX; Entries", 60, 0, 300)
    h_nroads2 = TH1F("nroads2", "; # roads/tower/BX; Entries", 60, 0, 300)
    h_dntracks = TH1F("dntracks", "; change in # good tracks/tower/BX; Entries", 8, -5.5, 2.5)

    # Loop over events
    for ievt, (evt1, evt2) in enumerate(izip(ttree1, ttree2)):

        if ievt == options.nevents:
            break

        if options.event_list:
            if ievt not in options.event_list:
                continue

        # AM roads
        nroads1 = count_nroads(evt1)
        nroads2 = count_nroads(evt2)

        #print("evt {0:3} -- Found {1:3} vs {2:3} roads".format(ievt, nroads1, nroads2))
        h_nroads1.Fill(nroads1)
        h_nroads2.Fill(nroads2)

        # AM tracks
        ntracks1 = count_ntracks(evt1)
        ntracks2 = count_ntracks(evt2)

        #print("evt {0:3} -- Found {1:3} vs {2:3} tracks".format(ievt, ntracks1, ntracks2))
        h_dntracks.Fill(ntracks2 - ntracks1)

        #print("-" * 16)

    draw(h_nroads1, h_nroads2, h_dntracks)

    set_branch_status(ttree1)
    set_branch_status(ttree2)
    return 0


# ______________________________________________________________________________
if __name__ == "__main__":

    init()

    main()
