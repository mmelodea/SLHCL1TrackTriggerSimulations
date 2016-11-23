#!/usr/bin/env python

from rootdrawing import *
from parser import *
from TrackParametersToTT import *

col = TColor.GetColor("#BA0000")
fcol = TColor.GetColor("#E58080")
col_res = TColor.GetColor("#008000")


# ______________________________________________________________________________
def use_TEfficiency():
    del ROOT.TEfficiency
    gROOT.ProcessLine(".L ./fixes/TEfficiency.cxx+")
    #gSystem.Load("./fixes/TEfficiency_cxx.so")

    global TEfficiency
    TEfficiency = ROOT.TEfficiency


# ______________________________________________________________________________
def drawer_book(options):
    histos = {}

    prefix = "efficiency_"
    ytitle = "#varepsilon_{synthetic}"

    hname = prefix + "pt"  # pt, 0-500 GeV
    histos[hname] = TEfficiency(hname, "; p_{T} [GeV]; %s" % ytitle, 100, 0., 500.)

    hname = prefix + "ppt"  # pt, 0-50 GeV
    histos[hname] = TEfficiency(hname, "; p_{T} [GeV]; %s" % ytitle, 100, 0., 50.)

    hname = prefix + "pppt"  # pt, 0-10 GeV
    histos[hname] = TEfficiency(hname, "; p_{T} [GeV]; %s" % ytitle, 100, 0., 10.)

    hname = prefix + "eta"
    histos[hname] = TEfficiency(hname, "; #eta; %s" % ytitle,        70, -0.4, 1.)

    hname = prefix + "etaStar"
    histos[hname] = TEfficiency(hname, "; #eta*; %s" % ytitle,       70, -0.4, 1.)

    hname = prefix + "phi"
    histos[hname] = TEfficiency(hname, "; #phi; %s" % ytitle,        80, 0.4, 2.)

    hname = prefix + "phiStar"
    histos[hname] = TEfficiency(hname, "; #phi*; %s" % ytitle,       80, 0.4, 2.)

    hname = prefix + "vz"
    histos[hname] = TEfficiency(hname, "; v_{z} [cm]; %s" % ytitle,  100, -30, 30)

    hname = prefix + "charge"
    histos[hname] = TEfficiency(hname, "; charge; %s" % ytitle,      5, -2.5, 2.5)

    prefix = "resolution_"
    xtitle = "p_{T} resolution / p_{T}"

    hname = prefix + "pt"  # pt, for 8-100 GeV
    histos[hname] = TH1F(hname, "; %s {8-100 GeV}" % xtitle, 100, -0.20, 0.20)

    hname = prefix + "ppt"  # pt, for 3-8 GeV
    histos[hname] = TH1F(hname, "; %s {3-8 GeV}" % xtitle, 100, -0.20, 0.20)

    hname = prefix + "pt_vs_pt"
    histos[hname] = TH2F(hname, "; c*p_{T} [GeV]; p_{T} resolution/p_{T}", 40, -100, 100, 100, -0.20, 0.20)

    hname = prefix + "pt_vs_eta"
    histos[hname] = TH2F(hname, "; #eta; p_{T} resolution/p_{T}", 40, 0, 0.8, 100, -0.20, 0.20)

    # Change binning
    if options.low_stat:
        prefix = "efficiency_"
        histos[prefix + "pt"  ].SetBins(50, 0., 500.)
        histos[prefix + "ppt" ].SetBins(50, 0., 50.)
        histos[prefix + "pppt"].SetBins(50, 0., 10.)

        prefix = "resolution_"
        histos[prefix + "pt_vs_pt"]  = TH2F(hname, "; c*p_{T} [GeV]; p_{T} resolution/p_{T}", 20, -100, 100, 100, -0.20, 0.20)
        histos[prefix + "pt_vs_eta"] = TH2F(hname, "; #eta; p_{T} resolution/p_{T}", 20, 0, 0.8, 100, -0.20, 0.20)

    if options.low_low_stat:
        prefix = "efficiency_"
        histos[prefix + "pt"  ].SetBins(20, 0., 500.)
        histos[prefix + "ppt" ].SetBins(20, 0., 50.)
        histos[prefix + "pppt"].SetBins(20, 0., 10.)

        prefix = "resolution_"
        histos[prefix + "pt_vs_pt"]  = TH2F(hname, "; c*p_{T} [GeV]; p_{T} resolution/p_{T}", 10, -100, 100, 100, -0.20, 0.20)
        histos[prefix + "pt_vs_eta"] = TH2F(hname, "; #eta; p_{T} resolution/p_{T}", 10, 0, 0.8, 100, -0.20, 0.20)

    # Error statistics
    for hname, h in histos.iteritems():
        if h.ClassName() == "TEfficiency":
            h.SetConfidenceLevel(0.682689492137)  # 1-sigma
            h.SetStatisticOption(0)  # kFCP
            #h.SetStatisticOption(6)  # kBUniform

    # Style
    for hname, h in histos.iteritems():
        h.SetLineWidth(2); h.SetMarkerSize(0)
        h.SetMarkerStyle(24); h.SetMarkerSize(1)
        h.SetMarkerColor(col); h.SetLineColor(col); h.SetFillColor(fcol)
    donotdelete.append(histos)
    return histos

def repr_cppvector(v):
    s = ",".join([str(x) for x in v])
    s = "[" + s + "]"
    return s

def drawer_project(tree, histos, options):
    tree.SetBranchStatus("*", 0)
    tree.SetBranchStatus("trkParts_pt"     , 1)
    tree.SetBranchStatus("trkParts_eta"    , 1)
    tree.SetBranchStatus("trkParts_phi"    , 1)
    #tree.SetBranchStatus("trkParts_vx"     , 1)
    #tree.SetBranchStatus("trkParts_vy"     , 1)
    tree.SetBranchStatus("trkParts_vz"     , 1)
    tree.SetBranchStatus("trkParts_charge" , 1)
    tree.SetBranchStatus("trkParts_primary", 1)
    tree.SetBranchStatus("trkParts_intime" , 1)
    tree.SetBranchStatus("trkParts_signal" , 1)
    tree.SetBranchStatus("trkParts_pdgId"  , 1)
    #tree.SetBranchStatus("AMTTRoads_patternRef"  , 1)
    #tree.SetBranchStatus("AMTTRoads_stubRefs"    , 1)
    #tree.SetBranchStatus("AMTTTracks_invPt"      , 1)
    #tree.SetBranchStatus("AMTTTracks_phi0"       , 1)
    #tree.SetBranchStatus("AMTTTracks_cottheta"   , 1)
    #tree.SetBranchStatus("AMTTTracks_z0"         , 1)
    tree.SetBranchStatus("AMTTTracks_pt"         , 1)
    tree.SetBranchStatus("AMTTTracks_eta"        , 1)
    #tree.SetBranchStatus("AMTTTracks_chi2"       , 1)
    #tree.SetBranchStatus("AMTTTracks_ndof"       , 1)
    tree.SetBranchStatus("AMTTTracks_synMatchCat", 1)
    tree.SetBranchStatus("AMTTTracks_synTpId"    , 1)
    tree.SetBranchStatus("AMTTTracks_patternRef" , 1)
    #tree.SetBranchStatus("AMTTTracks_stubRefs"   , 1)

    # Loop over events
    numEntries, denomEntries = 0, 0

    for ievt, evt in enumerate(tree):
        if (ievt == options.nentries):  break

        if (ievt % 100 == 0):  print "Processing event: %i" % ievt

        # Loop over tracking particles
        nparts_all = evt.trkParts_primary.size()
        trkparts = {}

        for ipart in xrange(nparts_all):
            if options.pu == 0:  # single-track events
                if ipart > 0:
                    break

            charge  = evt.trkParts_charge [ipart]
            primary = evt.trkParts_primary[ipart]
            intime  = evt.trkParts_intime [ipart]
            signal  = evt.trkParts_signal [ipart]

            if not (charge!=0 and primary and intime):
                continue

            if options.signal and not signal:
                continue

            pt      = evt.trkParts_pt     [ipart]
            eta     = evt.trkParts_eta    [ipart]
            phi     = evt.trkParts_phi    [ipart]
            #vx      = evt.trkParts_vx     [ipart]
            #vy      = evt.trkParts_vy     [ipart]
            vz      = evt.trkParts_vz     [ipart]
            pdgId   = evt.trkParts_pdgId  [ipart]

            #aux_TT = TrackParametersToTT(phi, float(charge)/pt , eta, vz)
            #if aux_TT != options.tower:
            #    continue

            trkparts[ipart] = (pt, eta, phi, vz, charge, pdgId)
            if options.verbose:  print ievt, "part ", ipart, trkparts[ipart]

        # Loop over reconstructed tracks
        ntracks_all = evt.AMTTTracks_patternRef.size()
        trkparts_trigger = {}
        trkparts_trigger_vars = {}

        for itrack in xrange(ntracks_all):
            trigger = False

            patternRef  = evt.AMTTTracks_patternRef [itrack]
            synMatchCat = evt.AMTTTracks_synMatchCat[itrack]
            if patternRef < options.maxPatterns:
                if synMatchCat == 1:  # good track
                    trigger = True

            track_pt    = evt.AMTTTracks_pt         [itrack]
            track_eta   = evt.AMTTTracks_eta        [itrack]
            #track_chi2  = evt.AMTTTracks_chi2       [itrack]
            #track_ndof  = evt.AMTTTracks_ndof       [itrack]
            synTpId     = evt.AMTTTracks_synTpId    [itrack]

            if trigger:
                trkparts_trigger[synTpId] = True
                trkparts_trigger_vars[synTpId] = (track_pt, track_eta)

            if options.verbose > 1:  print ievt, "track", itrack, track_pt, track_eta, synMatchCat, synTpId
            if options.verbose and trigger:  print ievt, "good track", itrack, track_pt, track_eta, synMatchCat, synTpId

        #for k, v in trkparts_trigger.iteritems():
        #    assert k in trkparts

        # ______________________________________________________________________
        # Fill the histograms

        if options.verbose > 2: print "trigger tower", options.ptmin, options.ptmax, options.phimin, options.phimax, options.etamin, options.etamax, options.vzmin, options.vzmax

        for k, v in trkparts.iteritems():
            pt, eta, phi, vz, charge, pdgId = v
            phiStar = get_phiStar_from_phi(phi, float(charge)/pt)
            etaStar = get_etaStar_from_eta(eta, vz, float(charge)/pt)

            # Trigger has fired for this tracking particle
            trigger = k in trkparts_trigger

            if options.verbose:  print ievt, "part ", k, v, phiStar, etaStar, trigger

            # Check N-1 conditions
            prefix = "efficiency_"
            if options.phimin < phiStar < options.phimax and options.etamin < etaStar < options.etamax and options.vzmin < vz < options.vzmax:
                histos[prefix + "pt"      ].Fill(trigger, pt)
                histos[prefix + "ppt"     ].Fill(trigger, pt)
                histos[prefix + "pppt"    ].Fill(trigger, pt)
            if options.ptmin < pt < options.ptmax and options.phimin < phiStar < options.phimax and options.vzmin < vz < options.vzmax:
                histos[prefix + "eta"     ].Fill(trigger, eta)
                histos[prefix + "etaStar" ].Fill(trigger, etaStar)
            if options.ptmin < pt < options.ptmax and options.etamin < etaStar < options.etamax and options.vzmin < vz < options.vzmax:
                histos[prefix + "phi"     ].Fill(trigger, phi)
                histos[prefix + "phiStar" ].Fill(trigger, phiStar)
            if options.ptmin < pt < options.ptmax and options.phimin < phiStar < options.phimax and options.etamin < etaStar < options.etamax:
                histos[prefix + "vz"      ].Fill(trigger, vz)
            if options.ptmin < pt < options.ptmax and options.phimin < phiStar < options.phimax and options.etamin < etaStar < options.etamax and options.vzmin < vz < options.vzmax:
                histos[prefix + "charge"  ].Fill(trigger, charge)
                denomEntries += 1
                if trigger:  numEntries += 1

            # If trigger, also find pT resolution
            if trigger:
                track_pt, track_eta = trkparts_trigger_vars[k]

                prefix = "resolution_"
                if options.phimin < phiStar < options.phimax and options.etamin < etaStar < options.etamax and options.vzmin < vz < options.vzmax:
                    ptres = (track_pt - pt) / pt
                    histos[prefix + "pt_vs_pt"].Fill(charge*pt, ptres)
                    if 8 < pt < 100:
                        histos[prefix + "pt"      ].Fill(ptres)
                    if 3 < pt < 8:
                        histos[prefix + "ppt"     ].Fill(ptres)

                if options.ptmin < pt < options.ptmax and options.phimin < phiStar < options.phimax and options.vzmin < vz < options.vzmax:
                    ptres = (track_pt - pt) / pt
                    histos[prefix + "pt_vs_eta"].Fill(eta, ptres)

        if options.verbose:  print ievt, numEntries, denomEntries

    # Bookkeeping
    prefix = "efficiency_"
    histos[prefix + "pt"      ].numEntries = numEntries
    histos[prefix + "pt"      ].denomEntries = denomEntries

    tree.SetBranchStatus("*", 1)
    return

def drawer_draw(histos, options):
    def displayGaus(h, ignorebins=0, scalebox=(1.,1.)):
        # Display gaussian fit parameters
        fitxmin, fitxmax = h.GetBinCenter(ignorebins+1), h.GetBinCenter(h.GetNbinsX()-1-ignorebins+1)
        r = h.Fit("gaus", "INSQ", "", fitxmin, fitxmax)
        p1, e1, p2, e2 = r.Parameter(1), r.ParError(1), r.Parameter(2), r.ParError(2)

        gPad.Modified(); gPad.Update()
        ps = h.FindObject("stats")
        ps.SetName("mystats")

        n = 2
        newX1NDC = ps.GetX2NDC() - (ps.GetX2NDC() - ps.GetX1NDC()) * scalebox[0]
        newY1NDC = ps.GetY2NDC() - ((ps.GetY2NDC() - ps.GetY1NDC()) / 5 * (5 + n)) * scalebox[1]
        ps.SetX1NDC(newX1NDC)
        ps.SetY1NDC(newY1NDC)

        ps.AddText("Gaus #mu = %6.4g" % (p1))
        ps.AddText("Gaus #sigma = %6.4g" % (p2))
        #ps.AddText("Gaus #mu = %6.4g #pm %6.4g" % (p1, e1))
        #ps.AddText("Gaus #sigma = %6.4g #pm %6.4g" % (p2, e2))
        h.stats = [p1, e1, p2, e2]

        h.SetStats(0)
        #gPad.Modified(); gPad.Update()
        ps.Draw()

    tline.SetLineStyle(1)
    tline2 = TLine()
    tline2.SetLineColor(4)
    tline2.SetLineWidth(2)
    tline2.SetLineStyle(2)

    for hname, h in histos.iteritems():
        if h.ClassName() == "TEfficiency":
            ymin, ymax = 0., 1.2
            h1 = h.GetCopyTotalHisto(); h1.SetName(h1.GetName()+"_frame"); h1.Reset()
            h1.SetMinimum(ymin); h1.SetMaximum(ymax)
            h1.SetNdivisions(524, "Y")  # 24 grid lines from 0 to 1.2
            h1.SetStats(0); h1.Draw()

            xmin, xmax = h1.GetXaxis().GetXmin(), h1.GetXaxis().GetXmax()
            tline.DrawLine(xmin, 1.0, xmax, 1.0)
            if hname.endswith("pt"):
                tline2.DrawLine(options.ptmin, ymin, options.ptmin, ymax)
            elif hname.endswith("phi") or hname.endswith("phiStar"):
                tline2.DrawLine(options.phimin, ymin, options.phimin, ymax)
                tline2.DrawLine(options.phimax, ymin, options.phimax, ymax)
            elif hname.endswith("eta") or hname.endswith("etaStar"):
                tline2.DrawLine(options.etamin, ymin, options.etamin, ymax)
                tline2.DrawLine(options.etamax, ymin, options.etamax, ymax)
            elif hname.endswith("vz"):
                tline2.DrawLine(options.vzmin, ymin, options.vzmin, ymax)
                tline2.DrawLine(options.vzmax, ymin, options.vzmax, ymax)

            h.gr = h.CreateGraph()
            h.gr.Draw("p")
            gPad.RedrawAxis("g"); h1.SetNdivisions(512, "Y")  # only change num of labels, not the grid lines

        elif h.ClassName() == "TH1F":
            ymax = h.GetMaximum() * 1.4
            h.SetMinimum(0); h.SetMaximum(ymax)
            h.SetLineColor(col_res); h.SetFillStyle(0)
            h.Draw("hist")
            displayGaus(h, ignorebins=10)

        elif h.ClassName() == "TH2F":
            ymax = 0.05
            nbinsx = h.GetNbinsX()
            gr1 = TGraphAsymmErrors(nbinsx)
            gr1.SetName(hname + "_gr1"); setattr(h, "gr1", gr1)
            gr1.SetMarkerSize(1); gr1.SetMarkerColor(col_res)
            gr2 = TGraphAsymmErrors(nbinsx)
            gr2.SetName(hname + "_gr2"); setattr(h, "gr2", gr2)
            gr2.SetMarkerSize(1); gr2.SetMarkerColor(col_res)

            for ibin in xrange(nbinsx):
                hpy = h.ProjectionY(hname + "_py%i" % ibin, ibin+1, ibin+1); setattr(h, "py%i" % ibin, hpy)
                x = h.GetXaxis().GetBinCenter(ibin+1)

                if hpy.GetEntries() < 30:
                    print "WARNING: not enough entries in %s: %i" % (hpy.GetName(), hpy.GetEntries())
                    gr1.SetPoint(ibin, x, 999999)
                    gr2.SetPoint(ibin, x, 999999)
                    continue

                if hpy.GetEntries() < 100:
                    print "INFO: rebinning in %s: %i" % (hpy.GetName(), hpy.GetEntries())
                    hpy.Rebin(2)

                assert(hpy.GetNbinsX() > 20)
                fitxmin, fitxmax = hpy.GetBinCenter(10+1), hpy.GetBinCenter(hpy.GetNbinsX()-1-10+1)
                r = hpy.Fit("gaus", "INSQ", "", fitxmin, fitxmax)
                p1, e1, p2, e2 = r.Parameter(1), r.ParError(1), r.Parameter(2), r.ParError(2)
                gr1.SetPoint(ibin, x, p1)
                gr1.SetPointError(ibin, 0, 0, e1, e1)
                gr2.SetPoint(ibin, x, p2)
                gr2.SetPointError(ibin, 0, 0, e2, e2)

            hpx2 = h.ProjectionX(hname + "_px2"); hpx2.Reset(); setattr(h, "px2", hpx2)
            hpx2.SetYTitle("p_{T} resolution/p_{T}")
            hpx2.SetMinimum(0); hpx2.SetMaximum(ymax)
            hpx2.SetStats(0); hpx2.Draw()
            gr2.Draw("p")

        CMS_label()
        save(options.outdir, "%s_%s" % (hname, options.ss), dot_root=True)
    return


def drawer_sitrep(histos, options):
    print "--- SITREP ---------------------------------------------------------"
    print "--- Using tt{0}, pu{1}, ss={2}, maxPatterns={3}".format(options.tower, options.pu, options.ss, options.maxPatterns)

    prefix = "efficiency_"
    h = histos[prefix + "pt"]
    print "eff={0}/{1}={2:.5f} +/- {3:.5f}".format(h.numEntries, h.denomEntries, float(h.numEntries)/h.denomEntries, (1./h.denomEntries) * sqrt(float(h.numEntries) * (1. - float(h.numEntries)/h.denomEntries)))


# ______________________________________________________________________________
# Main function
def main(options):

    # Init
    drawerInit = DrawerInit()
    tchain = TChain("ntupler/tree", "")
    tchain.AddFileInfoList(options.tfilecoll.GetList())

    gStyle.SetNdivisions(510, "Y")
    gStyle.SetEndErrorSize(0)
    gStyle.SetPadGridX(True)
    gStyle.SetPadGridY(True)
    use_TEfficiency()

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
    parser.add_argument("ss", help="short name of superstrip definition (e.g. sf1_nz8)")
    parser.add_argument("--maxPatterns", type=int, default=999999999, help="number of patterns to reach the desired coverage")
    parser.add_argument("--minPt", type=float, default=3, help="min pT for gen particle (default: %(default)s)")
    parser.add_argument("--low-stat", action="store_true", help="low statistics (default: %(default)s)")
    parser.add_argument("--low-low-stat", action="store_true", help="low low statistics (default: %(default)s)")

    # Parse default arguments
    options = parser.parse_args()
    parse_drawer_options(options)
    options.ptmin = options.minPt

    # Call the main function
    main(options)
