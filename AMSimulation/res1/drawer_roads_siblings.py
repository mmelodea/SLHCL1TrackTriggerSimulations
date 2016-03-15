#!/usr/bin/env python

from rootdrawing import *
from parser import *

# Configurations
#col  = TColor.GetColor("#1f78b4")  # mu0
#fcol = TColor.GetColor("#a6cee3")  # mu0

col  = TColor.GetColor("#e31a1c")  # nu140
fcol = TColor.GetColor("#fb9a99")  # nu140

#col  = TColor.GetColor("#6a3d9a")  # tttt140
#fcol = TColor.GetColor("#cab2d6")  # tttt140


# ______________________________________________________________________________
def getHitBits(stubRefs):
    for i in xrange(6):
        if stubRefs.at(i).size() == 0:
            return i + 1
    return 0

# ______________________________________________________________________________
def drawer_book(options):
    histos = {}

    hname = "nroads_per_event"
    nbins, xmin, xmax = 200, 0., 200.*options.xscale
    #nbins, xmin, xmax = 40, 0., 40.*options.xscale
    histos[hname] = TH1F(hname, "; # roads/tower/BX"                , nbins, xmin, xmax)

    hname = "ncombinations_per_road"
    nbins, xmin, xmax = 400, 0., 400.*options.xscale
    #nbins, xmin, xmax = 40, 0., 40.*options.xscale
    histos[hname] = TH1F(hname, "; # combinations/road/tower/BX"    , nbins, xmin, xmax)

    hname = "ncombinations_per_event"
    nbins, xmin, xmax = 800, 0., 1600.*options.xscale
    #nbins, xmin, xmax = 150, 0., 150.*options.xscale
    histos[hname] = TH1F(hname, "; # combinations/tower/BX"         , nbins, xmin, xmax)

    hname = "nsuperstrips_per_road"
    nbins, xmin, xmax = 20, 0., 20.
    #nbins, xmin, xmax = 20, 0., 20.
    histos[hname] = TH1F(hname, "; # superstrips/road/tower/BX"     , nbins, xmin, xmax)

    hname = "nstubs_per_superstrip"
    nbins, xmin, xmax = 50, 0., 50.
    #nbins, xmin, xmax = 20, 0., 20.
    histos[hname] = TH1F(hname, "; # stubs/superstrip/road/tower/BX", nbins, xmin, xmax)

    hname = "nstubs_per_road"
    nbins, xmin, xmax = 50, 0., 50.
    #nbins, xmin, xmax = 20, 0., 20.
    histos[hname] = TH1F(hname, "; # stubs/road/tower/BX"           , nbins, xmin, xmax)

    hname = "nstubs_per_event"
    nbins, xmin, xmax = 500, 0., 500.
    #nbins, xmin, xmax = 50, 0., 50.
    histos[hname] = TH1F(hname, "; # stubs/tower/BX"                , nbins, xmin, xmax)

    for i in xrange(6):
        hname = "nstubs_per_layer_%i" % i
        nbins, xmin, xmax = 50, 0., 50.
        #nbins, xmin, xmax = 20, 0., 20.
        histos[hname] = TH1F(hname, "; # stubs/layer/road/tower/BX" , nbins, xmin, xmax)

    hname = "road_hitbits"
    nbins, xmin, xmax = 10, 0, 10.
    histos[hname] = TH1F(hname, "; road hitbits" , nbins, xmin, xmax)

    hname = "combination_hitbits"
    nbins, xmin, xmax = 10, 0, 10.
    histos[hname] = TH1F(hname, "; combination hitbits" , nbins, xmin, xmax)

    hname = "nsiblings_per_event"
    nbins, xmin, xmax = 50, 0., 50.
    histos[hname] = TH1F(hname, "; # siblings/tower/BX"             , nbins, xmin, xmax)

    hname = "nfamilies_per_event"
    nbins, xmin, xmax = 100, 0., 100.
    histos[hname] = TH1F(hname, "; # families/tower/BX"             , nbins, xmin, xmax)
    hname = "nfamilies_sib1_per_event"
    nbins, xmin, xmax = 100, 0., 100.
    histos[hname] = TH1F(hname, "; # families/tower/BX"             , nbins, xmin, xmax)
    hname = "nfamilies_sib2_per_event"
    nbins, xmin, xmax = 100, 0., 100.
    histos[hname] = TH1F(hname, "; # families/tower/BX"             , nbins, xmin, xmax)

    hname = "road_hitbits_thesibling"
    nbins, xmin, xmax = 10, 0, 10.
    histos[hname] = TH1F(hname, "; road hitbits" , nbins, xmin, xmax)
    hname = "road_hitbits_nsiblings_7"
    nbins, xmin, xmax = 10, 0, 10.
    histos[hname] = TH1F(hname, "; road hitbits" , nbins, xmin, xmax)
    hname = "road_hitbits_nsiblings_6"
    nbins, xmin, xmax = 10, 0, 10.
    histos[hname] = TH1F(hname, "; road hitbits" , nbins, xmin, xmax)
    hname = "road_hitbits_nsiblings_5"
    nbins, xmin, xmax = 10, 0, 10.
    histos[hname] = TH1F(hname, "; road hitbits" , nbins, xmin, xmax)
    hname = "road_hitbits_nsiblings_4"
    nbins, xmin, xmax = 10, 0, 10.
    histos[hname] = TH1F(hname, "; road hitbits" , nbins, xmin, xmax)
    hname = "road_hitbits_nsiblings_3"
    nbins, xmin, xmax = 10, 0, 10.
    histos[hname] = TH1F(hname, "; road hitbits" , nbins, xmin, xmax)
    hname = "road_hitbits_nsiblings_2"
    nbins, xmin, xmax = 10, 0, 10.
    histos[hname] = TH1F(hname, "; road hitbits" , nbins, xmin, xmax)
    hname = "road_hitbits_nsiblings_1"
    nbins, xmin, xmax = 10, 0, 10.
    histos[hname] = TH1F(hname, "; road hitbits" , nbins, xmin, xmax)

    # Style
    for hname, h in histos.iteritems():
        h.SetLineWidth(2); h.SetMarkerSize(0)
        h.SetLineColor(col); h.SetFillColor(fcol)
        if h.ClassName() == "TH1F":
            binwidth = (h.GetXaxis().GetXmax() - h.GetXaxis().GetXmin())/h.GetNbinsX()
            h.SetYTitle("Entries / %.1f" % binwidth)

            if "road_hitbits" in hname or "combination_hitbits" in hname:
                h.logy = False
            else:
                h.logy = True
    donotdelete.append(histos)
    return histos


def drawer_project(tree, histos, options):
    tree.SetBranchStatus("*", 0)
    tree.SetBranchStatus("AMTTRoads_patternRef"   , 1)
    tree.SetBranchStatus("AMTTRoads_tower"        , 1)
    tree.SetBranchStatus("AMTTRoads_nstubs"       , 1)
    tree.SetBranchStatus("AMTTRoads_superstripIds", 1)
    tree.SetBranchStatus("AMTTRoads_stubRefs"     , 1)

    # Loop over events
    for ievt, evt in enumerate(tree):
        if (ievt == options.nentries):  break

        nroads_per_event = 0
        ncombinations_per_event = 0
        stubmap = {}  # per event

        if options.has_roads and evt.AMTTRoads_patternRef.size() == 0:
            continue

        # Loop over roads
        for patternRef, tower, nstubs, superstripIds, stubRefs in izip(evt.AMTTRoads_patternRef, evt.AMTTRoads_tower, evt.AMTTRoads_nstubs, evt.AMTTRoads_superstripIds, evt.AMTTRoads_stubRefs):

            if tower == options.tower and patternRef < options.npatterns:

                ssidmap = {}  # per road

                # superstripIds[i] is the i-th superstrip ID in the pattern (or road)
                # stubRefs[i][j] is the j-th stub REF in the i-th superstrip in the road
                l = 0  # layer i
                for ssid, ssid_stubRefs in izip(superstripIds, stubRefs):
                    for stub in ssid_stubRefs:
                        ssidmap[(l,ssid)] = ssidmap.get((l,ssid), 0) + 1

                        stubmap[stub] = stubmap.get(stub, 0) + 1

                    nstubs_per_layer = ssid_stubRefs.size()
                    histos["nstubs_per_layer_%i" % l].Fill(nstubs_per_layer)  # a superstrip can enter more than once

                    l += 1

                nsuperstrips_per_road = len(ssidmap)
                #assert(nsuperstrips_per_road == 6)
                histos["nsuperstrips_per_road"].Fill(nsuperstrips_per_road)

                road_hitbits = getHitBits(stubRefs)
                histos["road_hitbits"].Fill(road_hitbits)

                nstubs_per_road = 0
                ncombinations_per_road = 1

                # Loop over k=(l,ssid), v=count in ssidmap
                for k, v in ssidmap.iteritems():
                    nstubs_per_superstrip = v
                    nstubs_per_road += v
                    if v != 0:  # if no stub in the superstrip
                        ncombinations_per_road *= v

                    histos["nstubs_per_superstrip"].Fill(nstubs_per_superstrip)

                assert(nstubs_per_road == nstubs)
                histos["nstubs_per_road"].Fill(nstubs_per_road)

                histos["ncombinations_per_road"].Fill(ncombinations_per_road)

                histos["combination_hitbits"].Fill(road_hitbits, ncombinations_per_road)

                nroads_per_event += 1
                ncombinations_per_event += ncombinations_per_road

        assert(nroads_per_event <= evt.AMTTRoads_stubRefs.size())
        histos["nroads_per_event"].Fill(nroads_per_event)

        nstubs_per_event = len(stubmap)
        histos["nstubs_per_event"].Fill(nstubs_per_event)

        histos["ncombinations_per_event"].Fill(ncombinations_per_event)


        # Stupid clustering algorithm
        def is_sibling(isuperstripIds, jsuperstripIds):
            count = 0
            for i, j in izip(isuperstripIds, jsuperstripIds):
                if i != j:
                    count += 1
            if count == 0:
                print "ERROR: isuperstripIds == jsuperstripIds"
            return count <= 1

        # Pythonify
        roads_superstripIds = []
        roads_superstripIds_include = []
        roads_hitBits = []

        for patternRef, tower, nstubs, superstripIds, stubRefs in izip(evt.AMTTRoads_patternRef, evt.AMTTRoads_tower, evt.AMTTRoads_nstubs, evt.AMTTRoads_superstripIds, evt.AMTTRoads_stubRefs):
            if tower == options.tower and patternRef < options.npatterns:
                roads_superstripIds.append([])
                roads_superstripIds_include.append(True)
                roads_hitBits.append(getHitBits(stubRefs))

                for ssid, ssid_stubRefs in izip(superstripIds, stubRefs):
                    roads_superstripIds[-1].append(ssid)

        # Verbose
        #for iroad in xrange(len(roads_superstripIds)):
        #    print iroad, roads_superstripIds[iroad], roads_hitBits[iroad]

        # Build family
        roads_family = []
        k = 0
        while sum(roads_superstripIds_include) > 0:
            most_popular_iroad = -1
            most_popular_iroad_siblings = []

            for iroad in xrange(len(roads_superstripIds)):
                if not roads_superstripIds_include[iroad]: continue

                siblings = []

                for jroad in xrange(len(roads_superstripIds)):
                    if iroad == jroad:  continue
                    if not roads_superstripIds_include[jroad]: continue

                    if is_sibling(roads_superstripIds[iroad], roads_superstripIds[jroad]):
                        siblings.append(jroad)

                if most_popular_iroad == -1 or len(siblings) > len(most_popular_iroad_siblings):
                    most_popular_iroad = iroad
                    most_popular_iroad_siblings = siblings

            roads_family.append((most_popular_iroad, most_popular_iroad_siblings))
            roads_superstripIds_include[most_popular_iroad] = False
            for jroad in most_popular_iroad_siblings:
                roads_superstripIds_include[jroad] = False

            k += 1

        # Verbose
        #for kroad in xrange(len(roads_family)):
        #    print kroad, roads_family[kroad]
        print "Number of families: %i" % (len(roads_family))
        nsiblings = [len(k[1]) for k in roads_family[:3]]
        nsiblings += [-1, -1, -1]
        print "Number of siblings (most to least): %i, %i, %i, ..." % (nsiblings[0], nsiblings[1], nsiblings[2])

        # Check
        check_nroads = 0
        for kroad in xrange(len(roads_family)):
            check_nroads += 1
            check_nroads += len(roads_family[kroad][1])
        assert(check_nroads == len(roads_superstripIds))

        # Plot
        nfamilies = len(roads_family)
        histos["nfamilies_per_event"].Fill(nfamilies)

        nfamilies_sib1 = 0
        nfamilies_sib2 = 0
        for kroad in xrange(nfamilies):
            nsiblings = len(roads_family[kroad][1])
            if nsiblings > 0:
                nfamilies_sib1 += 1
            if nsiblings > 1:
                nfamilies_sib2 += 1
        histos["nfamilies_sib1_per_event"].Fill(nfamilies_sib1)
        histos["nfamilies_sib2_per_event"].Fill(nfamilies_sib2)

        for kroad in xrange(nfamilies):
            nsiblings = len(roads_family[kroad][1])
            histos["nsiblings_per_event"].Fill(nsiblings)

            sroad = roads_family[kroad][0]
            road_hitbits = roads_hitBits[sroad]
            histos["road_hitbits_thesibling"].Fill(road_hitbits)

            for sroad in roads_family[kroad][1]:
                road_hitbits = roads_hitBits[sroad]

                if nsiblings == 7:
                    histos["road_hitbits_nsiblings_7"].Fill(road_hitbits)
                elif nsiblings == 6:
                    histos["road_hitbits_nsiblings_6"].Fill(road_hitbits)
                elif nsiblings == 5:
                    histos["road_hitbits_nsiblings_5"].Fill(road_hitbits)
                elif nsiblings == 4:
                    histos["road_hitbits_nsiblings_4"].Fill(road_hitbits)
                elif nsiblings == 3:
                    histos["road_hitbits_nsiblings_3"].Fill(road_hitbits)
                elif nsiblings == 2:
                    histos["road_hitbits_nsiblings_2"].Fill(road_hitbits)
                elif nsiblings == 1:
                    histos["road_hitbits_nsiblings_1"].Fill(road_hitbits)

    tree.SetBranchStatus("*", 1)
    return


def drawer_draw(histos, options):
    def parse_ss(ss):
        if "lu" in ss:
            ss = ss.replace("lu", "")
        elif "sf" in ss:
            ss = ss.replace("sf", "sf=").replace("_nz", ",nz=")
            ss = ss.replace("0p", "0.")
        ss = ss.replace("_", " ")
        return ss

    def display_quantiles(h, in_quantiles=[0.95,0.99], scalebox=(1.,1.)):
        # Display one-sided confidence intervals, a.k.a quantiles
        n = len(in_quantiles)
        in_quantiles = array('d', in_quantiles)
        quantiles = array('d', [0. for i in xrange(n)])
        h.GetQuantiles(n, quantiles, in_quantiles)

        gPad.Modified(); gPad.Update()
        ps = h.FindObject("stats").Clone("mystats")

        newX1NDC = ps.GetX2NDC() - (ps.GetX2NDC() - ps.GetX1NDC()) * scalebox[0]
        newY1NDC = ps.GetY2NDC() - ((ps.GetY2NDC() - ps.GetY1NDC()) / 5 * (5 + n)) * scalebox[1]
        ps.SetX1NDC(newX1NDC)
        ps.SetY1NDC(newY1NDC)

        for iq, q in enumerate(in_quantiles):
            ps.AddText("%i%% CI = %6.4g" % (int(q*100), quantiles[iq]))
        h.stats = [h.GetMean()] + quantiles.tolist()

        h.SetStats(0)
        #gPad.Modified(); gPad.Update()
        ps.Draw()
        return ps

    def plot_quantiles(h):
        hquantile = h.Clone(h.GetName() + "_quantile")
        hquantile.Reset()
        hquantile.SetYTitle("Percentile of entries #geq x-value")
        for i in xrange(h.GetNbinsX(), 0, -1):
            if i == h.GetNbinsX():
                hquantile.SetBinContent(i, h.GetBinContent(i))
            else:
                hquantile.SetBinContent(i, hquantile.GetBinContent(i+1)+h.GetBinContent(i))
        hquantile.Scale(1.0/h.GetBinContent(1))
        return hquantile

    for hname, h in histos.iteritems():
        if "per_layer" in hname:
            continue

        if h.logy:
            h.SetMaximum(h.GetMaximum() * 14); h.SetMinimum(0.5)
        else:
            h.SetMaximum(h.GetMaximum() * 1.4); h.SetMinimum(0.)
        h.SetStats(1); h.Draw("hist")
        gPad.SetLogy(h.logy)
        ps = display_quantiles(h)

        tlatex.DrawLatex(0.6, 0.185, "%s [%.0fK bank]" % (parse_ss(options.ss), options.npatterns*1e-3))
        CMS_label()
        save(options.outdir, "%s_%s" % (hname, options.ss), dot_root=True)

        if hname in ["nroads_per_event", "ncombinations_per_event"]:
            hquantile = plot_quantiles(h)
            hquantile.SetMaximum(hquantile.GetMaximum() * 14); hquantile.SetMinimum(0.0001)
            hquantile.SetStats(0); hquantile.Draw("hist")
            gPad.SetLogy(True)
            ps.Draw()

            tlatex.DrawLatex(0.6, 0.185, "%s [%.0fK bank]" % (parse_ss(options.ss), options.npatterns*1e-3))
            CMS_label()
            save(options.outdir, "%s_%s" % (hquantile.GetName(), options.ss), dot_root=True)

    # Specialized: nstubs_per_layer_%i plots
    if True:
        c1 = TCanvas("c1", "c1")
        lmargin, rmargin, bmargin, tmargin = (gStyle.GetPadLeftMargin(), gStyle.GetPadRightMargin(), gStyle.GetPadBottomMargin(), gStyle.GetPadTopMargin())
        width, height = (1.0 - lmargin - rmargin) * 0.5, (1.0 - bmargin - tmargin) * 0.333
        pad1 = TPad("pad1", "pad1", lmargin+width*0, 1.0-tmargin-height*1, lmargin+width*1, 1.0-tmargin-height*0)
        pad1.SetNumber(1); pad1.SetMargin(0.0, 0.0, 0.0, 0.0)
        pad1.Draw()
        pad2 = TPad("pad2", "pad2", lmargin+width*0, 1.0-tmargin-height*2, lmargin+width*1, 1.0-tmargin-height*1)
        pad2.SetNumber(2); pad2.SetMargin(0.0, 0.0, 0.0, 0.0)
        pad2.Draw()
        pad3 = TPad("pad3", "pad3", lmargin+width*0, 1.0-tmargin-height*3, lmargin+width*1, 1.0-tmargin-height*2)
        pad3.SetNumber(3); pad3.SetMargin(0.0, 0.0, 0.0, 0.0)
        pad3.Draw()
        pad4 = TPad("pad4", "pad4", lmargin+width*1, 1.0-tmargin-height*1, lmargin+width*2, 1.0-tmargin-height*0)
        pad4.SetNumber(4); pad4.SetMargin(0.0, 0.0, 0.0, 0.0)
        pad4.Draw()
        pad5 = TPad("pad5", "pad5", lmargin+width*1, 1.0-tmargin-height*2, lmargin+width*2, 1.0-tmargin-height*1)
        pad5.SetNumber(5); pad5.SetMargin(0.0, 0.0, 0.0, 0.0)
        pad5.Draw()
        pad6 = TPad("pad6", "pad6", lmargin+width*1, 1.0-tmargin-height*3, lmargin+width*2, 1.0-tmargin-height*2)
        pad6.SetNumber(6); pad6.SetMargin(0.0, 0.0, 0.0, 0.0)
        pad6.Draw()
        c1.cd()
        c1.Modified(); c1.Update()

        temps = (tlatex.GetTextSize(),)
        tlatex.SetTextSize(0.08)

        ymax = -1

        for i in xrange(6):
            hname = "nstubs_per_layer_%i" % i
            h = histos[hname]
            c1.cd(i+1)

            # Style
            if i/3 == 0:
                #h.GetYaxis().SetTitleOffset(0.90)
                h.GetYaxis().SetTitleSize(h.GetYaxis().GetTitleSize() * 2)
                h.GetYaxis().SetLabelSize(h.GetYaxis().GetLabelSize() * 2)
            else:
                h.GetYaxis().SetTitle("")
                h.GetYaxis().SetLabelSize(0)

            if i%3 == 2:
                #h.GetXaxis().SetTitleOffset(1.10)
                h.GetXaxis().SetTitleSize(h.GetXaxis().GetTitleSize() * 2)
                h.GetXaxis().SetLabelSize(h.GetXaxis().GetLabelSize() * 2)
            else:
                h.GetXaxis().SetTitle("")
                h.GetXaxis().SetLabelSize(0)

            if ymax == -1:
                ymax = h.GetMaximum()

            if h.logy:
                h.SetMaximum(ymax * 14); h.SetMinimum(0.5)
            else:
                h.SetMaximum(ymax * 1.4); h.SetMinimum(0.)
            h.SetStats(1); h.Draw("hist")
            gPad.SetLogy(h.logy)
            display_quantiles(h, scalebox=(2.,2.))

            tlatex.DrawLatex(0.5, 0.260, "Layer %i" % i)
            tlatex.DrawLatex(0.5, 0.185, "%s [%.0fK bank]" % (parse_ss(options.ss), options.npatterns*1e-3))

        c1.cd(0)
        CMS_label()
        hname = "nstubs_per_layer"
        save(options.outdir, "%s_%s" % (hname, options.ss), dot_root=True)

        tlatex.SetTextSize(temps[0])
    return


def drawer_sitrep(histos, options):
    print "--- SITREP ---------------------------------------------------------"
    print "--- Using tt{0}, pu{1}, ss={2}, npatterns={3}".format(options.tower, options.pu, options.ss, options.npatterns)
    print "--- Variable, mean, 95%% CI, 99%% CI:"
    h = histos["nroads_per_event"]
    print "nroads per event\t{0:6.4g}\t{1:6.4g}\t{2:6.4g}".format(*h.stats)
    h = histos["nstubs_per_road"]
    print "nstubs per road \t{0:6.4g}\t{1:6.4g}\t{2:6.4g}".format(*h.stats)
    h = histos["ncombinations_per_road"]
    print "ncombs per road \t{0:6.4g}\t{1:6.4g}\t{2:6.4g}".format(*h.stats)
    h = histos["ncombinations_per_event"]
    print "ncombs per event\t{0:6.4g}\t{1:6.4g}\t{2:6.4g}".format(*h.stats)


# ______________________________________________________________________________
# Main function
def main(options):

    # Init
    drawerInit = DrawerInit()
    tchain = TChain("ntupler/tree", "")
    tchain.AddFileInfoList(options.tfilecoll.GetList())

    # Process
    histos = drawer_book(options)
    drawer_project(tchain, histos, options)
    drawer_draw(histos, options)
    drawer_sitrep(histos, options)


# ______________________________________________________________________________
if __name__ == '__main__':

    # Setup argument parser
    parser = argparse.ArgumentParser()

    # Add default arguments
    add_drawer_arguments(parser)

    # Add more arguments
    parser.add_argument("ss", help="short name of superstrip definition (e.g. ss256)")
    parser.add_argument("npatterns", type=int, help="number of patterns to reach the desired coverage")
    parser.add_argument("--coverage", type=float, default=0.95, help="desired coverage (default: %(default)s)")
    parser.add_argument("--minPt", type=float, default=2, help="min pT (default: %(default)s)")
    parser.add_argument("--xscale", type=float, default=1, help="scale factor for the x-axis range (default: %(default)s)")
    parser.add_argument("--has-roads", action="store_true", help="only events with at least 1 road (default: %(default)s)")

    # Parse default arguments
    options = parser.parse_args()
    parse_drawer_options(options)
    options.ptmin = options.minPt

    # Call the main function
    main(options)
