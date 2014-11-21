#!/usr/bin/env python

from rootdrawing import *
from roothelper import *
from array import array
import numpy as np
import os

# FIXME: find quantiles

# ______________________________________________________________________________
# Configurations

sections = {}
sections["fixed"     ] = False
sections["projective"] = True

drawerInit = DrawerInit()
gStyle.SetPadRightMargin(0.1)

EOS = "/eos/uscms/store/user/l1upgrades/SLHC/GEN/620_SLHC12p1_results_PU140/"
DATE = "20141108"

fixed_settings_tt27 = [
    ("ss32"  ,   42318379, 0.9026),
    ("ss64"  ,   17008261, 0.9729),
    ("ss128" ,    6253924, 0.9930),
    ("ss256" ,    2039355, 0.9981),
    ("ss512" ,     628771, 0.9996),
    ("ss1024",     213156, 0.9999),
]

fixed_settings_tt27_pt3 = [
    ("ss32"  ,   29343005,  0.9055),
    ("ss64"  ,   11854473,  0.9734),
    ("ss128" ,    4377734,  0.9929),
    ("ss256" ,    1453179,  0.9980),
    ("ss512" ,     463751,  0.9994),
    ("ss1024",     158710,  0.9999),
]

projective_settings_tt27 = [
    ("600x0" ,   15808528,  0.9763),
    ("400x0" ,    4285974,  0.9952),
    ("200x0" ,     541342,  0.9997),
    ("200x1" ,    3205364,  0.9955),
    ("100x2" ,    1615816,  0.9981),
    ("20x10" ,    1129484,  0.9989),
]

projective_settings_tt27_pt3 = [
    ("600x0" ,    5806521,  0.9901),
    ("400x0" ,    1717472,  0.9980),
    ("200x0" ,     263322,  0.9998),
    ("200x1" ,    1883636,  0.9970),
    ("100x2" ,    1031066,  0.9985),
    ("20x10" ,     763513,  0.9990),
]

imgdir = "figures_roads/"

# Number of events
nentries = 1000
#nentries = 100

chain = TChain("ntupler/tree", "")

if not imgdir.endswith("/"):  imgdir += "/"
if gSystem.AccessPathName(imgdir):
    gSystem.mkdir(imgdir)

# ______________________________________________________________________________
# Load

def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d) if f.endswith(".root")]


# ______________________________________________________________________________
# Fixed

if sections["fixed"]:

    gStyle.SetTitleSize(0.05, "Y")
    latex.SetTextSize(0.05)

    def bookRate(settings, results):
        in_quantiles = array('d', [0.0015, 0.025, 0.16, 0.5, 0.84, 0.975, 0.9985])
        quantiles = array('d', [0., 0., 0., 0., 0., 0., 0.])
        h1 = TH1F("h1", "", 10000, 0, 10000)
        h2 = TH1F("h2", "", 10000, 0, 10000)

        superstrips = []
        graphs = {}

        for ss, npatterns, coverage in settings:
            chain.Reset()
            infiles = listdir_fullpath(EOS + "/" + (results % ss) + "_" + DATE)
            for f in infiles:
                chain.Add(f)

            # Sample 200 points
            npoints = 200
            every = npatterns / npoints
            xvalues, yvalues1, yvalues2 = [], [], []
            for i in xrange(0,npatterns,every):
                xvalues.append(i); yvalues1.append([]); yvalues2.append([])
            xvalues.append(npatterns-1); yvalues1.append([]); yvalues2.append([])
            npoints = len(xvalues)

            # Loop over events
            tree = chain
            tree.SetBranchStatus("*", 0)
            tree.SetBranchStatus("AMTTRoads_nsuperstrips"     , 1)
            tree.SetBranchStatus("AMTTRoads_bankIndex"        , 1)
            tree.SetBranchStatus("AMTTRoads_stubSuperstripIds", 1)

            for i_ievt, ievt in enumerate(tree):
                if (i_ievt == nentries):  break

                # Init counters
                for i in xrange(npoints):
                    yvalues1[i].append(0)
                    yvalues2[i].append(0)

                for nsuperstrips, bankIndex, superstripIds in izip(ievt.AMTTRoads_nsuperstrips, ievt.AMTTRoads_bankIndex, ievt.AMTTRoads_stubSuperstripIds):
                    # Calculate # of combinations
                    ssidmap = {}
                    for ssid in superstripIds:
                        ssidmap[ssid] = ssidmap.get(ssid, 0) + 1
                    assert(nsuperstrips == len(ssidmap))

                    ncombinations_per_road = 1
                    for k, v in ssidmap.iteritems():
                        ncombinations_per_road *= v

                    # Sum # of roads and # of combinations
                    for i in xrange(npoints):
                        if bankIndex < xvalues[i]:
                            yvalues1[i][-1] += 1
                            yvalues2[i][-1] += ncombinations_per_road

                assert(yvalues1[-1][-1] == len(ievt.AMTTRoads_bankIndex))

            # Find mean, 2 sigma quantiles
            yvalues1_mean, yvalues1_p2sigma = [], []
            yvalues2_mean, yvalues2_p2sigma = [], []
            for i in xrange(npoints):
                h1.Reset(); h2.Reset()
                for j in yvalues1[i]:
                    h1.Fill(j)
                for j in yvalues2[i]:
                    h2.Fill(j)

                h1.GetQuantiles(len(quantiles), quantiles, in_quantiles)
                yvalues1_mean.append(h1.GetMean())
                yvalues1_p2sigma.append(quantiles[5])

                h2.GetQuantiles(len(quantiles), quantiles, in_quantiles)
                yvalues2_mean.append(h2.GetMean())
                yvalues2_p2sigma.append(quantiles[5])

            for gname, yvalues in [("gr1_mean_%s"    % ss, yvalues1_mean),
                                   ("gr1_p2sigma_%s" % ss, yvalues1_p2sigma),
                                   ("gr2_mean_%s"    % ss, yvalues2_mean),
                                   ("gr2_p2sigma_%s" % ss, yvalues2_p2sigma)]:
                gr = TGraph(npoints, array('d', xvalues), array('d', yvalues))
                gr.SetName(gname)
                graphs[gname] = gr

            superstrips.append((ss, npatterns, yvalues1_mean[-1], yvalues1_p2sigma[-1], yvalues2_mean[-1], yvalues2_p2sigma[-1]))
            print superstrips[-1]

            tree.SetBranchStatus("*", 0)
        return (superstrips, graphs)

    def drawRate(superstrips, graphs, xmin=0, xmax=1e8, tower="tt27"):
        ymax1, ymax2, ymax3, ymax4 = superstrips[3][2:6]  # ss256

        hframe = TH1F("hframe", "; # of patterns; <# of roads>", 100, xmin, xmax)
        hframe.SetStats(0); hframe.SetMinimum(0); hframe.SetMaximum(ymax1 * 1.2)
        hframe.SetNdivisions(510, "Y")

        # Style
        for k, v in graphs.iteritems():
            v.SetLineWidth(2); v.SetLineStyle(1); v.SetMarkerSize(0)
            for i, ss in enumerate(superstrips):
                if ss[0] in k:
                    v.SetLineColor(paletteSet1[i])

        # Legend
        moveLegend(0.66,0.15,0.96,0.45); legend.Clear()
        for i, ss in enumerate(superstrips):
            gr = graphs["gr1_mean_%s" % ss[0]]
            legend.AddEntry(gr, ss[0], "l")

        # Draw
        hframe.Draw()
        for i, ss in enumerate(superstrips):
            gr = graphs["gr1_mean_%s" % ss[0]]
            gr.Draw("C")
        legend.Draw()
        CMS_label()
        save(imgdir, "fixed_roads_mean_%s" % tower, dot_root=True)

        # Zoom in
        hframe.Draw()
        hframe.GetXaxis().SetRangeUser(0, xmax/50)
        for i, ss in enumerate(superstrips):
            gr = graphs["gr1_mean_%s" % ss[0]]
            gr.Draw("C")
        legend.Draw()
        CMS_label()
        save(imgdir, "fixed_roads_mean_zoom_%s" % tower)

        # Draw
        hframe.GetYaxis().SetTitle("+2 #sigma quantile of # of roads"); hframe.SetMaximum(ymax2 * 1.2); hframe.GetXaxis().UnZoom()
        hframe.Draw()
        for i, ss in enumerate(superstrips):
            gr = graphs["gr1_p2sigma_%s" % ss[0]]
            gr.Draw("C")
        legend.Draw()
        CMS_label()
        save(imgdir, "fixed_roads_p2sigma_%s" % tower, dot_root=True)

        # Zoom in
        hframe.Draw()
        hframe.GetXaxis().SetRangeUser(0, xmax/50)
        for i, ss in enumerate(superstrips):
            gr = graphs["gr1_p2sigma_%s" % ss[0]]
            gr.Draw("C")
        legend.Draw()
        CMS_label()
        save(imgdir, "fixed_roads_p2sigma_zoom_%s" % tower)

        # Draw
        hframe.GetYaxis().SetTitle("<# of fit combinations>"); hframe.SetMaximum(ymax3 * 1.2); hframe.GetXaxis().UnZoom()
        hframe.Draw()
        for i, ss in enumerate(superstrips):
            gr = graphs["gr2_mean_%s" % ss[0]]
            gr.Draw("C")
        legend.Draw()
        CMS_label()
        save(imgdir, "fixed_fits_mean_%s" % tower, dot_root=True)

        # Zoom in
        hframe.Draw()
        hframe.GetXaxis().SetRangeUser(0, xmax/50)
        for i, ss in enumerate(superstrips):
            gr = graphs["gr2_mean_%s" % ss[0]]
            gr.Draw("C")
        legend.Draw()
        CMS_label()
        save(imgdir, "fixed_fits_mean_zoom_%s" % tower)

        # Draw
        hframe.GetYaxis().SetTitle("+2 #sigma quantile of # of fit combinations"); hframe.SetMaximum(ymax4 * 1.2); hframe.GetXaxis().UnZoom()
        hframe.Draw()
        for i, ss in enumerate(superstrips):
            gr = graphs["gr2_p2sigma_%s" % ss[0]]
            gr.Draw("C")
        legend.Draw()
        CMS_label()
        save(imgdir, "fixed_fits_p2sigma_%s" % tower, dot_root=True)

        # Zoom in
        hframe.Draw()
        hframe.GetXaxis().SetRangeUser(0, xmax/50)
        for i, ss in enumerate(superstrips):
            gr = graphs["gr2_p2sigma_%s" % ss[0]]
            gr.Draw("C")
        legend.Draw()
        CMS_label()
        save(imgdir, "fixed_fits_p2sigma_zoom_%s" % tower)

        donotdelete = [hframe]
        return donotdelete


    # Barrel 2 GeV
    settings = fixed_settings_tt27
    results = "Neutrino_sp16_%s_tt27"
    (superstrips, graphs) = bookRate(settings, results)
    d = drawRate(superstrips, graphs, xmax=5e7, tower="tt27")

    # Barrel 3 GeV
    #settings = fixed_settings_tt27_pt3
    #results = "Neutrino_sp16_%s_tt27_pt3"
    #(superstrips, graphs) = bookRate(settings, results)
    #d = drawRate(superstrips, graphs, xmax=5e7, tower="tt27_pt3")


# ______________________________________________________________________________
# Projective

if sections["projective"]:

    gStyle.SetTitleSize(0.05, "Y")
    latex.SetTextSize(0.05)

    def bookRate(settings, results):
        in_quantiles = array('d', [0.0015, 0.025, 0.16, 0.5, 0.84, 0.975, 0.9985])
        quantiles = array('d', [0., 0., 0., 0., 0., 0., 0.])
        h1 = TH1F("h1", "", 10000, 0, 10000)
        h2 = TH1F("h2", "", 10000, 0, 10000)

        superstrips = []
        graphs = {}

        for ss, npatterns, coverage in settings:
            chain.Reset()
            infiles = listdir_fullpath(EOS + "/" + (results % ss) + "_" + DATE)
            for f in infiles:
                chain.Add(f)

            # Sample 200 points
            npoints = 200
            every = npatterns / npoints
            xvalues, yvalues1, yvalues2 = [], [], []
            for i in xrange(0,npatterns,every):
                xvalues.append(i); yvalues1.append([]); yvalues2.append([])
            xvalues.append(npatterns-1); yvalues1.append([]); yvalues2.append([])
            npoints = len(xvalues)

            # Loop over events
            tree = chain
            tree.SetBranchStatus("*", 0)
            tree.SetBranchStatus("AMTTRoads_nsuperstrips"     , 1)
            tree.SetBranchStatus("AMTTRoads_bankIndex"        , 1)
            tree.SetBranchStatus("AMTTRoads_stubSuperstripIds", 1)

            for i_ievt, ievt in enumerate(tree):
                if (i_ievt == nentries):  break

                # Init counters
                for i in xrange(npoints):
                    yvalues1[i].append(0)
                    yvalues2[i].append(0)

                for nsuperstrips, bankIndex, superstripIds in izip(ievt.AMTTRoads_nsuperstrips, ievt.AMTTRoads_bankIndex, ievt.AMTTRoads_stubSuperstripIds):
                    # Calculate # of combinations
                    ssidmap = {}
                    for ssid in superstripIds:
                        ssidmap[ssid] = ssidmap.get(ssid, 0) + 1
                    assert(nsuperstrips == len(ssidmap))

                    ncombinations_per_road = 1
                    for k, v in ssidmap.iteritems():
                        ncombinations_per_road *= v

                    # Sum # of roads and # of combinations
                    for i in xrange(npoints):
                        if bankIndex < xvalues[i]:
                            yvalues1[i][-1] += 1
                            yvalues2[i][-1] += ncombinations_per_road

                assert(yvalues1[-1][-1] == len(ievt.AMTTRoads_bankIndex))

            # Find mean, 2 sigma quantiles
            yvalues1_mean, yvalues1_p2sigma = [], []
            yvalues2_mean, yvalues2_p2sigma = [], []
            for i in xrange(npoints):
                h1.Reset(); h2.Reset()
                for j in yvalues1[i]:
                    h1.Fill(j)
                for j in yvalues2[i]:
                    h2.Fill(j)

                h1.GetQuantiles(len(quantiles), quantiles, in_quantiles)
                yvalues1_mean.append(h1.GetMean())
                yvalues1_p2sigma.append(quantiles[5])

                h2.GetQuantiles(len(quantiles), quantiles, in_quantiles)
                yvalues2_mean.append(h2.GetMean())
                yvalues2_p2sigma.append(quantiles[5])

            for gname, yvalues in [("gr1_mean_%s"    % ss, yvalues1_mean),
                                   ("gr1_p2sigma_%s" % ss, yvalues1_p2sigma),
                                   ("gr2_mean_%s"    % ss, yvalues2_mean),
                                   ("gr2_p2sigma_%s" % ss, yvalues2_p2sigma)]:
                gr = TGraph(npoints, array('d', xvalues), array('d', yvalues))
                gr.SetName(gname)
                graphs[gname] = gr

            superstrips.append((ss, npatterns, yvalues1_mean[-1], yvalues1_p2sigma[-1], yvalues2_mean[-1], yvalues2_p2sigma[-1]))
            print superstrips[-1]

            tree.SetBranchStatus("*", 0)
        return (superstrips, graphs)

    def drawRate(superstrips, graphs, xmin=0, xmax=1e8, tower="tt27"):
        ymax1, ymax2, ymax3, ymax4 = superstrips[1][2:6]  # 400x1

        hframe = TH1F("hframe", "; # of patterns; <# of roads>", 100, xmin, xmax)
        hframe.SetStats(0); hframe.SetMinimum(0); hframe.SetMaximum(ymax1 * 1.2)
        hframe.SetNdivisions(510, "Y")

        # Style
        for k, v in graphs.iteritems():
            v.SetLineWidth(2); v.SetLineStyle(1); v.SetMarkerSize(0)
            for i, ss in enumerate(superstrips):
                if ss[0] in k:
                    v.SetLineColor(paletteSet1[i])

        # Legend
        moveLegend(0.66,0.15,0.96,0.45); legend.Clear()
        for i, ss in enumerate(superstrips):
            gr = graphs["gr1_mean_%s" % ss[0]]
            legend.AddEntry(gr, ss[0], "l")

        # Draw
        hframe.Draw()
        for i, ss in enumerate(superstrips):
            gr = graphs["gr1_mean_%s" % ss[0]]
            gr.Draw("C")
        legend.Draw()
        CMS_label()
        save(imgdir, "projective_roads_mean_%s" % tower, dot_root=True)

        # Zoom in
        hframe.Draw()
        hframe.GetXaxis().SetRangeUser(0, xmax/50)
        for i, ss in enumerate(superstrips):
            gr = graphs["gr1_mean_%s" % ss[0]]
            gr.Draw("C")
        legend.Draw()
        CMS_label()
        save(imgdir, "projective_roads_mean_zoom_%s" % tower)

        # Draw
        hframe.GetYaxis().SetTitle("+2 #sigma quantile of # of roads"); hframe.SetMaximum(ymax2 * 1.2); hframe.GetXaxis().UnZoom()
        hframe.Draw()
        for i, ss in enumerate(superstrips):
            gr = graphs["gr1_p2sigma_%s" % ss[0]]
            gr.Draw("C")
        legend.Draw()
        CMS_label()
        save(imgdir, "projective_roads_p2sigma_%s" % tower, dot_root=True)

        # Zoom in
        hframe.Draw()
        hframe.GetXaxis().SetRangeUser(0, xmax/50)
        for i, ss in enumerate(superstrips):
            gr = graphs["gr1_p2sigma_%s" % ss[0]]
            gr.Draw("C")
        legend.Draw()
        CMS_label()
        save(imgdir, "projective_roads_p2sigma_zoom_%s" % tower)

        # Draw
        hframe.GetYaxis().SetTitle("<# of fit combinations>"); hframe.SetMaximum(ymax3 * 1.2); hframe.GetXaxis().UnZoom()
        hframe.Draw()
        for i, ss in enumerate(superstrips):
            gr = graphs["gr2_mean_%s" % ss[0]]
            gr.Draw("C")
        legend.Draw()
        CMS_label()
        save(imgdir, "projective_fits_mean_%s" % tower, dot_root=True)

        # Zoom in
        hframe.Draw()
        hframe.GetXaxis().SetRangeUser(0, xmax/50)
        for i, ss in enumerate(superstrips):
            gr = graphs["gr2_mean_%s" % ss[0]]
            gr.Draw("C")
        legend.Draw()
        CMS_label()
        save(imgdir, "projective_fits_mean_zoom_%s" % tower)

        # Draw
        hframe.GetYaxis().SetTitle("+2 #sigma quantile of # of fit combinations"); hframe.SetMaximum(ymax4 * 1.2); hframe.GetXaxis().UnZoom()
        hframe.Draw()
        for i, ss in enumerate(superstrips):
            gr = graphs["gr2_p2sigma_%s" % ss[0]]
            gr.Draw("C")
        legend.Draw()
        CMS_label()
        save(imgdir, "projective_fits_p2sigma_%s" % tower, dot_root=True)

        # Zoom in
        hframe.Draw()
        hframe.GetXaxis().SetRangeUser(0, xmax/50)
        for i, ss in enumerate(superstrips):
            gr = graphs["gr2_p2sigma_%s" % ss[0]]
            gr.Draw("C")
        legend.Draw()
        CMS_label()
        save(imgdir, "projective_fits_p2sigma_zoom_%s" % tower)

        donotdelete = [hframe]
        return donotdelete


    # Barrel 2 GeV
    #settings = projective_settings_tt27
    #results = "Neutrino_lu%s_tt27"
    #(superstrips, graphs) = bookRate(settings, results)
    #d = drawRate(superstrips, graphs, xmax=5e7, tower="tt27")

    # Barrel 3 GeV
    settings = projective_settings_tt27_pt3
    results = "Neutrino_lu%s_tt27_pt3"
    (superstrips, graphs) = bookRate(settings, results)
    d = drawRate(superstrips, graphs, xmax=5e7, tower="tt27_pt3")

