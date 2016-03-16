#!/usr/bin/env python

from rootdrawing import *
from parser import *


# ______________________________________________________________________________
sf = [
0.7,
1.0,
1.26,
2.0,
4.0,
6.0,
]

fountain = [
146000,
58100,
33800,
12300,
2900,
1300,
]

optimized = [
151700,
60100,
34700,
12400,
3000,
1300,
]

optimized_1 = [
149200,
59300,
34200,
12300,
3000,
1300,
]

optimized_2 = [
147400,
58500,
33700,
12100,
2900,
1300,
]

optimized_3 = [
142400,
56500,
32500,
11700,
2800,
1200,
]

optimized_4 = [
140000,
55300,
32000,
11500,
2700,
1200,
]


# ______________________________________________________________________________
def drawer_draw():

    hcanvas = TH1F("hcanvas", "; SF; # patterns @ 95% cov", 100, 0., 8.)
    gr_fountain = TGraph(len(sf))
    gr_optimized = TGraph(len(sf))
    gr_optimized_1 = TGraph(len(sf))
    gr_optimized_2 = TGraph(len(sf))
    gr_optimized_3 = TGraph(len(sf))
    gr_optimized_4 = TGraph(len(sf))
    donotdelete.append([hcanvas, gr_fountain, gr_optimized, gr_optimized_1, gr_optimized_2, gr_optimized_3, gr_optimized_4])

    gr_fountain.SetLineWidth(2); gr_fountain.SetLineColor(palette[0]); gr_fountain.SetMarkerColor(palette[0])
    gr_optimized.SetLineWidth(2); gr_optimized.SetLineColor(palette[1]); gr_optimized.SetMarkerColor(palette[1])
    gr_optimized_1.SetLineWidth(2); gr_optimized_1.SetLineColor(palette[2]); gr_optimized_1.SetMarkerColor(palette[2])
    gr_optimized_2.SetLineWidth(2); gr_optimized_2.SetLineColor(palette[3]); gr_optimized_2.SetMarkerColor(palette[3])
    gr_optimized_3.SetLineWidth(2); gr_optimized_3.SetLineColor(palette[4]); gr_optimized_3.SetMarkerColor(palette[4])
    gr_optimized_4.SetLineWidth(2); gr_optimized_4.SetLineColor(palette[5]); gr_optimized_4.SetMarkerColor(palette[5])

    for i in xrange(len(sf)):
        gr_fountain.SetPoint(i, sf[i], fountain[i])
        gr_optimized.SetPoint(i, sf[i], optimized[i])
        gr_optimized_1.SetPoint(i, sf[i], optimized_1[i])
        gr_optimized_2.SetPoint(i, sf[i], optimized_2[i])
        gr_optimized_3.SetPoint(i, sf[i], optimized_3[i])
        gr_optimized_4.SetPoint(i, sf[i], optimized_4[i])

    hcanvas.SetStats(0)
    hcanvas.SetMaximum(200000)
    hcanvas.Draw()

    for gr in [gr_fountain, gr_optimized, gr_optimized_1, gr_optimized_2, gr_optimized_3, gr_optimized_4]:
        gr.Draw("lp")

    moveLegend(0.64,0.64,0.94,0.94); tlegend.Clear()
    tlegend.AddEntry(gr_fountain, "fountain", "lp")
    tlegend.AddEntry(gr_optimized, "bugia", "lp")
    tlegend.AddEntry(gr_optimized_1, "L5x3", "lp")
    tlegend.AddEntry(gr_optimized_2, "L5x3 L10x2", "lp")
    tlegend.AddEntry(gr_optimized_3, "L5x2 L10x2", "lp")
    tlegend.AddEntry(gr_optimized_4, "bugia v2", "lp")
    tlegend.Draw()

    save(".", "npatterns_vs_sf")


# ______________________________________________________________________________
# Main function
def main(options):

    # Init
    drawerInit = DrawerInit()

    # Process
    drawer_draw()


# ______________________________________________________________________________
if __name__ == '__main__':

    # Setup argument parser
    parser = argparse.ArgumentParser()

    # Parse default arguments
    options = parser.parse_args()

    # Call the main function
    main(options)
