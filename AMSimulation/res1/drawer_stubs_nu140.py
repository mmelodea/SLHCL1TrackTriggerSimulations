#!/usr/bin/env python

from rootdrawing import *
from parser import *

# Configurations
col  = TColor.GetColor("#e31a1c")  # nu140
fcol = TColor.GetColor("#fb9a99")  # nu140

#col  = TColor.GetColor("#6a3d9a")  # tttt140
#fcol = TColor.GetColor("#cab2d6")  # tttt140

# Load trigger tower map
ttmap = json.load(open("../data/trigger_sector_map.json"), object_pairs_hook=convert_key_to_int)
ttrmap = get_reverse_map(ttmap)

connectmap = json.load(open("../data/module_connections.json"), object_pairs_hook=convert_key_to_int)


# ______________________________________________________________________________
def get_ds_cuts(ds):
    if ds == 0:
        # Default
        return [0, 2.5, 2.5, 3.0, 4.5, 5.5, 6.5]
    if ds == 1:
        # Default
        return [0, 2.5, 2.5, 3.0, 4.5, 5.5, 6.5]
    if ds == 2:
        # Taken from https://github.com/sviret/HL_LHC/blob/master/DataProduction/test/base/2GevTune.txt
        return [0, 2.0, 2.0, 2.5, 4.0, 5.5, 6.5]
    if ds == 3:
        # Taken from https://github.com/sviret/HL_LHC/blob/master/DataProduction/test/base/3GevTune.txt
        return [0, 1.5, 1.5, 2.0, 3.0, 3.5, 4.5]

# ______________________________________________________________________________
def drawer_book():
    histos = {}

    # TH1F
    hname = "ntrkParts_pt0p5"
    nbinsx, xmin, xmax = 100, 0., 100.
    histos[hname] = TH1F(hname, "; # trkParts/tower/BX {pt>0.5}", nbinsx, xmin, xmax)

    hname = "ntrkParts_pt1"
    nbinsx, xmin, xmax = 50, 0., 50.
    histos[hname] = TH1F(hname, "; # trkParts/tower/BX {pt>1}", nbinsx, xmin, xmax)

    hname = "ntrkParts_pt2"
    nbinsx, xmin, xmax = 50, 0., 50.
    histos[hname] = TH1F(hname, "; # trkParts/tower/BX {pt>2}", nbinsx, xmin, xmax)

    hname = "ntrkParts_pt3"
    nbinsx, xmin, xmax = 50, 0., 50.
    histos[hname] = TH1F(hname, "; # trkParts/tower/BX {pt>3}", nbinsx, xmin, xmax)

    hname = "ntrkParts_pt5"
    nbinsx, xmin, xmax = 50, 0., 50.
    histos[hname] = TH1F(hname, "; # trkParts/tower/BX {pt>5}", nbinsx, xmin, xmax)

    hname = "ntrkParts_pt10"
    nbinsx, xmin, xmax = 50, 0., 50.
    histos[hname] = TH1F(hname, "; # trkParts/tower/BX {pt>10}", nbinsx, xmin, xmax)

    hname = "nstubs_per_layer_%i" % 99
    nbinsx, xmin, xmax = 1200, 0., 1200.
    histos[hname] = TH1F(hname, "; # stubs/layer/tower/BX"        , nbinsx, xmin, xmax)

    for i in xrange(6):
        hname = "nstubs_per_layer_%i" % i
        nbinsx, xmin, xmax = 400, 0., 400.
        histos[hname] = TH1F(hname, "; # stubs/layer/tower/BX {L%i}" % i, nbinsx, xmin, xmax)

    hname = "nconnections"
    histos[hname] = TH1F(hname, "; # connections", 6, 0, 6)

    # Style
    for hname, h in histos.iteritems():
        h.SetLineWidth(2); h.SetMarkerSize(0)
        h.SetLineColor(col); h.SetFillColor(fcol)
        if h.ClassName() == "TH1F":
            binwidth = (h.GetXaxis().GetXmax() - h.GetXaxis().GetXmin())/h.GetNbinsX()
            h.SetYTitle("Entries / %.1f" % binwidth)
    donotdelete.append(histos)
    return histos


def drawer_project(tree, histos, options):
    tree.SetBranchStatus("*", 0)
    tree.SetBranchStatus("trkParts_pt"     , 1)
    tree.SetBranchStatus("trkParts_eta"    , 1)
    tree.SetBranchStatus("trkParts_phi"    , 1)
    tree.SetBranchStatus("trkParts_charge" , 1)
    tree.SetBranchStatus("trkParts_primary", 1)
    tree.SetBranchStatus("TTStubs_modId"   , 1)
    tree.SetBranchStatus("TTStubs_trigBend", 1)

    moduleIds_tt = set()
    for moduleId in ttmap[options.tower]:
        moduleIds_tt.add(moduleId)

    # Delta_s cuts
    dscuts = get_ds_cuts(3)
    dscuts = dscuts[1:]

    # Loop over events
    for ievt, evt in enumerate(tree):
        if (ievt == options.nentries):  break
        if (ievt % 100 == 0):  print "Processing event: %i" % ievt

        # ntrkParts
        pt_map = {}

        for ipart, (pt, eta, phi, charge, primary) in enumerate(izip(evt.trkParts_pt, evt.trkParts_eta, evt.trkParts_phi, evt.trkParts_charge, evt.trkParts_primary)):

            if not (charge!=0 and primary):
                continue

            if not (options.etamin < eta < options.etamax and options.phimin < phi < options.phimax):
                continue

            pt_map[ipart] = pt

        histos["ntrkParts_pt0p5" ].Fill(count_if(pt_map.values(), lambda x: x>0.5))
        histos["ntrkParts_pt1"   ].Fill(count_if(pt_map.values(), lambda x: x>1))
        histos["ntrkParts_pt2"   ].Fill(count_if(pt_map.values(), lambda x: x>2))
        histos["ntrkParts_pt3"   ].Fill(count_if(pt_map.values(), lambda x: x>3))
        histos["ntrkParts_pt5"   ].Fill(count_if(pt_map.values(), lambda x: x>5))
        histos["ntrkParts_pt10"  ].Fill(count_if(pt_map.values(), lambda x: x>10))

        # nstubs
        lay_map = {}

        for istub, (moduleId, ds) in enumerate(izip(evt.TTStubs_modId, evt.TTStubs_trigBend)):
            if moduleId not in moduleIds_tt:
                continue

            lay = decodeLayer(moduleId)

            if lay-5 < len(dscuts) and not (abs(ds) <= dscuts[lay-5]):
                continue

            lay_map[istub] = lay

        hname = "nstubs_per_layer_%i" % 99
        histos[hname].Fill(count_if(lay_map.values(), lambda x: 1))

        for i in xrange(6):
            hname = "nstubs_per_layer_%i" % i
            histos[hname].Fill(count_if(lay_map.values(), lambda x: (x-5) == i))

    hname = "nconnections"
    for k, v in connectmap.iteritems():
        n = 0
        if k in ttrmap:
            assert(len(v) == len(ttrmap[k]))
            n = len(v)
        histos[hname].Fill(n)

    tree.SetBranchStatus("*", 1)
    return


def drawer_draw(histos, options):
    def display_quantiles(h, in_quantiles=[0.95,0.99], scalebox=(1.,1.)):
        # Display one-sided confidence intervals, a.k.a quantiles
        n = len(in_quantiles)
        in_quantiles = array('d', in_quantiles)
        quantiles = array('d', [0. for i in xrange(n)])
        h.GetQuantiles(n, quantiles, in_quantiles)

        gPad.Modified(); gPad.Update()
        ps = h.FindObject("stats")
        ps.SetName("mystats")

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

    for hname, h in histos.iteritems():
        h.SetMaximum(h.GetMaximum() * 1.4); h.SetMinimum(0)
        h.Draw("hist")
        display_quantiles(h)

        if hname == "nconnections":
            h.SetMarkerSize(1.5)
            h.SetNdivisions(510, "XYZ")
            h.Draw("same hist text0")

        CMS_label()
        save(options.outdir, hname, dot_root=True)


def drawer_sitrep(histos, options):
    print "--- SITREP ---------------------------------------------------------"


# ______________________________________________________________________________
# Main function
def main(options):

    # Init
    drawerInit = DrawerInit()
    tchain = TChain("ntupler/tree", "")
    tchain.AddFileInfoList(options.tfilecoll.GetList())

    # Process
    histos = drawer_book()
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
    #parser.add_argument("ss", help="short name of superstrip definition (e.g. ss256)")
    #parser.add_argument("npatterns", type=int, help="number of patterns to reach the desired coverage")

    # Parse default arguments
    options = parser.parse_args()
    parse_drawer_options(options)

    # Call the main function
    main(options)
