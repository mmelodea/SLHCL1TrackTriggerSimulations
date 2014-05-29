#!/usr/bin/env python

from rootdrawing import *
import os

# ______________________________________________________________________________
# Global

latex.SetTextSize(0.04)

line1 = TLine()
line1.SetLineColor(40)
line1.SetLineStyle(2)

# ______________________________________________________________________________
# Configurations

sections = {}
sections["coverage_vs_pt"] = False
sections["coverage_vs_eta"] = False
sections["coverage_vs_eta_phi"] = True

drawerInit = DrawerInit()
chain = TChain("ntupler/tree", "")

EOS = "/eos/uscms/store/user/l1upgrades/SLHC/GEN/620_SLHC10_results/"
RES = "SingleMuPlusMinus_OneOverPt_vz0_%s"
settings = [
    "ss128_20140525",
    "ss64_20140525",
    "ss32_20140525",
    "ss16_20140525",
    "ss8_20140525",
    "ss128_m1_20140525",
    "ss64_m1_20140525",
    "ss32_m1_20140525",
    "ss16_m1_20140525",
    "ss8_m1_20140525",
]

imgdir = "figures_amsim/"
if not imgdir.endswith("/"):  imgdir += "/"
if gSystem.AccessPathName(imgdir):
    gSystem.mkdir(imgdir)


# ______________________________________________________________________________
# Profiles

histos = {}

tricolors = [kBlack, kRed, kBlue]
def bookProf(xtitle="", nbinsx=50, xlow=0, xup=50, ylow=0, yup=1, option=""):
    for i, j in enumerate(["prof0", "prof1", "prof2"]):
        h = TProfile(j, xtitle, nbinsx, xlow, xup, ylow, yup, option)
        h.SetFillStyle(0)
        h.SetMarkerStyle(24)
        h.SetMarkerSize(0.9)
        h.SetLineWidth(2)
        h.SetMarkerColor(tricolors[i])
        h.SetLineColor(tricolors[i])
        histos[j] = h

def drawProf():
    h = histos["prof0"]
    h.SetStats(0)
    h.SetMaximum(1.2)
    h.SetMinimum(0)
    h.Draw()
    histos["prof1"].Draw("same")
    histos["prof2"].Draw("same")

    xmin, xmax = h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax()
    line.DrawLine(xmin, 1.0, xmax, 1.0)
    line1.DrawLine(xmin, 0.4, xmax, 0.4)

    CMS_label()

def bookProf2D(xtitle="", nbinsx=50, xlow=0, xup=50, nbinsy=50, ylow=0, yup=50, zlow=0, zup=1, option=""):
    for i, j in enumerate(["prof0", "prof1", "prof2"]):
        h = TProfile2D(j, xtitle, nbinsx, xlow, xup, nbinsy, ylow, yup, zlow, zup, option)
        histos[j] = h


# ______________________________________________________________________________
# Coverage

if sections["coverage_vs_pt"]:

    for sett in settings[0:]:
        chain.Reset()

        bookProf("; p_{T} [GeV]; coverage", 100, 0., 100., 0., 1., "")

        result = RES % sett
        infiles = map(lambda x: EOS+"/"+result+"/"+x, os.listdir(EOS+"/"+result))
        print sett, len(infiles)

        for f in infiles:
            chain.Add(f)
        tree = chain

        tree.Project("prof0", "(@AMTTRoads_hash.size()!=0 && Alt$(AMTTRoads_hash[0],-1)>0):genParts_pt[0]", "abs(genParts_eta[0])<1.0", "prof goff")
        tree.Project("prof1", "(@AMTTRoads_hash.size()!=0 && Alt$(AMTTRoads_hash[0],-1)>0):genParts_pt[0]", "abs(genParts_eta[0])<1.7", "prof goff")
        tree.Project("prof2", "(@AMTTRoads_hash.size()!=0 && Alt$(AMTTRoads_hash[0],-1)>0):genParts_pt[0]", "abs(genParts_eta[0])<2.2", "prof goff")

        drawProf()

        leg1 = TLegend(0.18,0.78,0.60,0.92)
        leg1.SetFillStyle(0); leg1.SetLineColor(0); leg1.SetShadowColor(0); leg1.SetBorderSize(0)
        leg1.AddEntry(histos["prof0"], "Single muon |#eta| < 1.0")
        leg1.AddEntry(histos["prof1"], "Single muon |#eta| < 1.7")
        leg1.AddEntry(histos["prof2"], "Single muon |#eta| < 2.2")
        leg1.Draw()

        latex.DrawLatex(0.2, 0.74, "%s [3.2M bank]" % (sett[:-9].replace("_", " ")))

        save(imgdir, "coverage_vs_pt_%s" % sett)


if sections["coverage_vs_eta"]:

    for sett in settings[0:]:
        chain.Reset()

        bookProf("; #eta; coverage", 100, -2.5, 2.5, 0., 1., "")

        result = RES % sett
        infiles = map(lambda x: EOS+"/"+result+"/"+x, os.listdir(EOS+"/"+result))
        print sett, len(infiles)

        for f in infiles:
            chain.Add(f)
        tree = chain

        tree.Project("prof0", "(@AMTTRoads_hash.size()!=0 && Alt$(AMTTRoads_hash[0],-1)>0):genParts_eta[0]", "genParts_pt[0]>20", "prof goff")
        tree.Project("prof1", "(@AMTTRoads_hash.size()!=0 && Alt$(AMTTRoads_hash[0],-1)>0):genParts_eta[0]", "genParts_pt[0]>5", "prof goff")
        tree.Project("prof2", "(@AMTTRoads_hash.size()!=0 && Alt$(AMTTRoads_hash[0],-1)>0):genParts_eta[0]", "genParts_pt[0]>2", "prof goff")

        drawProf()

        leg1 = TLegend(0.18,0.78,0.60,0.92)
        leg1.SetFillStyle(0); leg1.SetLineColor(0); leg1.SetShadowColor(0); leg1.SetBorderSize(0)
        leg1.AddEntry(histos["prof0"], "Single muon p_{T} > 20")
        leg1.AddEntry(histos["prof1"], "Single muon p_{T} > 5")
        leg1.AddEntry(histos["prof2"], "Single muon p_{T} > 2")
        leg1.Draw()

        latex.DrawLatex(0.2, 0.74, "%s [3.2M bank]" % (sett[:-9].replace("_", " ")))

        save(imgdir, "coverage_vs_eta_%s" % sett)

# ______________________________________________________________________________
# Coverage (2D)

if sections["coverage_vs_eta_phi"]:

    for sett in settings[0:]:
        chain.Reset()

        bookProf2D("; #eta; #phi; coverage", 50, -2.5, 2.5, 64, -3.2, 3.2, 0., 1., "")

        result = RES % sett
        infiles = map(lambda x: EOS+"/"+result+"/"+x, os.listdir(EOS+"/"+result))
        print sett, len(infiles)

        for f in infiles:
            chain.Add(f)
        tree = chain

        tree.Project("prof0", "(@AMTTRoads_hash.size()!=0 && Alt$(AMTTRoads_hash[0],-1)>0):genParts_phi[0]:genParts_eta[0]", "genParts_pt[0]>20", "prof goff")
        #tree.Project("prof1", "(@AMTTRoads_hash.size()!=0 && Alt$(AMTTRoads_hash[0],-1)>0):genParts_phi[0]:genParts_eta[0]", "genParts_pt[0]>5", "prof goff")
        #tree.Project("prof2", "(@AMTTRoads_hash.size()!=0 && Alt$(AMTTRoads_hash[0],-1)>0):genParts_phi[0]:genParts_eta[0]", "genParts_pt[0]>2", "prof goff")

        h = histos["prof0"]
        h.SetStats(0)
        h.SetMaximum(1.2)
        h.SetMinimum(0)
        draw2D([h], stats=False)

        latex.DrawLatex(0.2, 0.78, "Single muon p_{T} > 20")
        latex.DrawLatex(0.2, 0.74, "%s [3.2M bank]" % (sett[:-9].replace("_", " ")))

        save(imgdir, "coverage_vs_eta_phi_%s" % sett)
