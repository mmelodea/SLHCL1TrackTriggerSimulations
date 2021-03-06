import unittest
import sys
from math import pi
from ROOT import TFile, TTree, gROOT, vector


class TestNTuple(unittest.TestCase):
    """Test the ntuple produced by NTupleTools"""

    infile = "test_ntuple.root"

    def setUp(self):
        gROOT.SetBatch(True)
        self.runOnMC = True
        self.tfile = TFile.Open(self.infile)
        self.ttree = self.tfile.Get("ntupler/tree")
        self.nevents = self.ttree.GetEntries()
        self.fvec1 = vector('float')()
        self.fvec2 = vector('float')()
        self.fvec3 = vector('float')()
        self.fvec4 = vector('float')()
        self.fvec5 = vector('float')()
        self.fvec6 = vector('float')()
        self.ivec1 = vector('int')()
        self.ivec2 = vector('int')()
        self.ivec3 = vector('int')()
        self.uvec1 = vector('unsigned int')()
        self.uvec2 = vector('unsigned int')()
        self.uvec3 = vector('unsigned int')()
        self.uvec4 = vector('unsigned int')()
        self.uvec5 = vector('unsigned int')()
        self.uvec6 = vector('unsigned int')()
        self.uvec7 = vector('unsigned int')()
        self.uvec8 = vector('unsigned int')()
        self.uvec9 = vector('unsigned int')()
        self.bvec1 = vector('bool')()
        self.bvec2 = vector('bool')()
        self.bvec3 = vector('bool')()
        self.bvec4 = vector('bool')()


    def tearDown(self):
        pass

    def test_nevents(self):
        self.assertEqual(self.nevents, 100)

    def test_eventInfo(self):
        tree = self.ttree
        for i in xrange(self.nevents):
            tree.GetEntry(i)
            run = getattr(tree, "run")
            event = getattr(tree, "event")
            lumi = getattr(tree, "lumi")
            print "Processing Run %i, Event %i, LumiSection %i." % (run, event, lumi)

    def test_genParticles(self):
        tree = self.ttree
        tree.SetBranchAddress("genParts_pt", self.fvec1)
        tree.SetBranchAddress("genParts_phi", self.fvec2)
        for i in xrange(self.nevents):
            tree.GetEntry(i)
            n = getattr(tree, "genParts_size")
            self.assertTrue(n > 0)
            for j in xrange(n):
                self.assertTrue(self.fvec1[j] > 0)
                self.assertTrue(abs(self.fvec2[j]) < pi+1e-6)
            with self.assertRaises(IndexError):
                print self.fvec1[n]

    #def test_genJets(self):
    #    tree = self.ttree
    #    tree.SetBranchAddress("genJets_pt", self.fvec1)
    #    tree.SetBranchAddress("genJets_phi", self.fvec2)
    #    for i in xrange(self.nevents):
    #        tree.GetEntry(i)
    #        n = getattr(tree, "genJets_size")
    #        self.assertTrue(n >= 0)  # can be zero
    #        for j in xrange(n):
    #            self.assertTrue(self.fvec1[j] > 0)
    #            self.assertTrue(abs(self.fvec2[j]) < pi+1e-6)
    #        with self.assertRaises(IndexError):
    #            print self.fvec1[n]

    #def test_genMET(self):
    #    tree = self.ttree
    #    for i in xrange(self.nevents):
    #        tree.GetEntry(i)
    #        fval1 = getattr(tree, "genMET_pt")
    #        fval2 = getattr(tree, "genMET_phi")
    #        self.assertTrue(fval1 > 0)
    #        self.assertTrue(abs(fval2) < pi+1e-6)

    def test_genEventInfo(self):
        tree = self.ttree
        for i in xrange(self.nevents):
            tree.GetEntry(i)
            fval1 = getattr(tree, "gen_trueNPV")
            fval2 = getattr(tree, "gen_weightMC")
            self.assertTrue(fval1 >= 0)
            self.assertTrue(fval2 >= 0)

    def test_simTracks(self):
        tree = self.ttree
        tree.SetBranchAddress("simTracks_pt", self.fvec1)
        tree.SetBranchAddress("simTracks_phi", self.fvec2)
        for i in xrange(self.nevents):
            tree.GetEntry(i)
            n = getattr(tree, "simTracks_size")
            self.assertTrue(n > 0)
            for j in xrange(n):
                self.assertTrue(self.fvec1[j] > 0)
                self.assertTrue(abs(self.fvec2[j]) < pi+1e-6)
            with self.assertRaises(IndexError):
                print self.fvec1[n]

    def test_simVertices(self):
        tree = self.ttree
        tree.SetBranchAddress("simVertices_vx", self.fvec1)
        tree.SetBranchAddress("simVertices_vz", self.fvec2)
        for i in xrange(self.nevents):
            tree.GetEntry(i)
            n = getattr(tree, "simVertices_size")
            self.assertTrue(n > 0)
            for j in xrange(n):
                self.assertTrue(abs(self.fvec1[j]) < 1000)
                self.assertTrue(abs(self.fvec2[j]) < 1200)
            with self.assertRaises(IndexError):
                print self.fvec1[n]

    def test_trackingParticles(self):
        tree = self.ttree
        tree.SetBranchAddress("trkParts_pt", self.fvec1)
        tree.SetBranchAddress("trkParts_phi", self.fvec2)
        for i in xrange(self.nevents):
            tree.GetEntry(i)
            n = getattr(tree, "trkParts_size")
            self.assertTrue(n > 0)
            for j in xrange(n):
                self.assertTrue(self.fvec1[j] > 0)
                self.assertTrue(abs(self.fvec2[j]) < pi+1e-6)
            with self.assertRaises(IndexError):
                print self.fvec1[n]

    def test_trackingVertices(self):
        tree = self.ttree
        tree.SetBranchAddress("trkVertices_vx", self.fvec1)
        tree.SetBranchAddress("trkVertices_vz", self.fvec2)
        for i in xrange(self.nevents):
            tree.GetEntry(i)
            n = getattr(tree, "trkVertices_size")
            self.assertTrue(n > 0)
            for j in xrange(n):
                self.assertTrue(abs(self.fvec1[j]) < 1000)
                self.assertTrue(abs(self.fvec2[j]) < 1200)
            with self.assertRaises(IndexError):
                print self.fvec1[n]

    def test_TTClusters(self):
        tree = self.ttree
        tree.SetBranchAddress("TTClusters_x", self.fvec1)
        tree.SetBranchAddress("TTClusters_y", self.fvec2)
        tree.SetBranchAddress("TTClusters_z", self.fvec3)
        tree.SetBranchAddress("TTClusters_coordx", self.fvec4)
        tree.SetBranchAddress("TTClusters_coordy", self.fvec5)
        tree.SetBranchAddress("TTClusters_nhits", self.uvec1)
        tree.SetBranchAddress("TTClusters_modId", self.uvec2)
        tree.SetBranchAddress("TTClusters_isGenuine", self.bvec1)
        tree.SetBranchAddress("TTClusters_isCombinatoric", self.bvec2)
        tree.SetBranchAddress("TTClusters_isUnknown", self.bvec3)

        barrel_z_divisions = [63, 55, 54, 24, 24, 24]
        barrel_phi_divisions = [16, 24, 34, 48, 62, 76]
        endcap_ring_divisions = [20, 24, 28, 28, 32, 36, 36, 40, 40, 52, 56, 64, 68, 76, 80]

        for i in xrange(self.nevents):
            tree.GetEntry(i)
            n = getattr(tree, "TTClusters_size")
            self.assertTrue(n > 0)
            for j in xrange(n):
                self.assertTrue(abs(self.fvec1[j]) < 110)
                self.assertTrue(abs(self.fvec2[j]) < 110)
                self.assertTrue(abs(self.fvec3[j]) < 300)
                self.assertTrue(self.fvec4[j] >= 0 and self.fvec4[j] < 1016)
                self.assertTrue(self.fvec5[j] >= 0 and self.fvec5[j] < 32)
                self.assertTrue(self.uvec1[j] >= 1)
                moduleId = self.uvec2[j]
                lay = moduleId / 10000
                lad = (moduleId / 100) % 100
                mod = moduleId % 100
                self.assertTrue(lay in [5,6,7,8,9,10,11,12,13,14,15,18,19,20,21,22])
                if moduleId < 110000:  # barrel
                    self.assertTrue(lad < barrel_phi_divisions[lay-5] and mod < barrel_z_divisions[lay-5])
                else:                  # endcap
                    self.assertTrue(lad < len(endcap_ring_divisions) and mod < endcap_ring_divisions[lad])
                mc = 0 + self.bvec1[j] + self.bvec2[j] + self.bvec3[j]
                self.assertTrue(1 == mc)
                mc = 0 + 1 * self.bvec1[j] + 2 * self.bvec2[j] + 3 * self.bvec3[j]
                self.assertTrue(1 <= mc <= 3)

                # test associated tp
                if self.bvec1[j]:
                    simPt1 = getattr(tree, "TTClusters_simPt")[j]
                    simPt2 = getattr(tree, "trkParts_pt")[getattr(tree, "TTClusters_tpId")[j]]
                    self.assertEqual(simPt1, simPt2)
                else:
                    simPt1 = getattr(tree, "TTClusters_simPt")[j]
                    self.assertEqual(simPt1, -99.)
            with self.assertRaises(IndexError):
                print self.fvec1[n]

    def test_TTStubs(self):
        tree = self.ttree
        tree.SetBranchAddress("TTStubs_x", self.fvec1)
        tree.SetBranchAddress("TTStubs_y", self.fvec2)
        tree.SetBranchAddress("TTStubs_z", self.fvec3)
        tree.SetBranchAddress("TTStubs_roughPt", self.fvec4)
        tree.SetBranchAddress("TTStubs_trigDist", self.fvec5)
        tree.SetBranchAddress("TTStubs_nhits", self.uvec1)
        tree.SetBranchAddress("TTStubs_modId", self.uvec2)
        tree.SetBranchAddress("TTStubs_isGenuine", self.bvec1)
        tree.SetBranchAddress("TTStubs_isCombinatoric", self.bvec2)
        tree.SetBranchAddress("TTStubs_isUnknown", self.bvec3)

        barrel_z_divisions = [63, 55, 54, 24, 24, 24]
        barrel_phi_divisions = [16, 24, 34, 48, 62, 76]
        endcap_ring_divisions = [20, 24, 28, 28, 32, 36, 36, 40, 40, 52, 56, 64, 68, 76, 80]

        for i in xrange(self.nevents):
            tree.GetEntry(i)
            n = getattr(tree, "TTStubs_size")
            self.assertTrue(n > 0)
            for j in xrange(n):
                self.assertTrue(abs(self.fvec1[j]) < 110)
                self.assertTrue(abs(self.fvec2[j]) < 110)
                self.assertTrue(abs(self.fvec3[j]) < 300)
                self.assertTrue(self.fvec4[j] > 0)
                self.assertTrue(abs(self.fvec5[j]) < 7.5)
                self.assertTrue(self.uvec1[j] >= 1)
                moduleId = self.uvec2[j]
                lay = moduleId / 10000
                lad = (moduleId / 100) % 100
                mod = moduleId % 100
                self.assertTrue(lay in [5,6,7,8,9,10,11,12,13,14,15,18,19,20,21,22])
                if moduleId < 110000:  # barrel
                    self.assertTrue(lad < barrel_phi_divisions[lay-5] and mod < barrel_z_divisions[lay-5])
                else:                  # endcap
                    self.assertTrue(lad < len(endcap_ring_divisions) and mod < endcap_ring_divisions[lad])
                mc = 0 + self.bvec1[j] + self.bvec2[j] + self.bvec3[j]
                self.assertTrue(1 == mc)
                mc = 0 + 1 * self.bvec1[j] + 2 * self.bvec2[j] + 3 * self.bvec3[j]
                self.assertTrue(1 <= mc <= 3)

                # test associated clusters
                geoId01 = getattr(tree, "TTStubs_geoId0")[j]
                geoId02 = getattr(tree, "TTClusters_geoId")[getattr(tree, "TTStubs_clusId0")[j]]
                geoId11 = getattr(tree, "TTStubs_geoId1")[j]
                geoId12 = getattr(tree, "TTClusters_geoId")[getattr(tree, "TTStubs_clusId1")[j]]
                stackId = getattr(tree, "TTStubs_stackId")[j]
                stackId0 = getattr(tree, "TTClusters_stackId")[getattr(tree, "TTStubs_clusId0")[j]]
                stackId1 = getattr(tree, "TTClusters_stackId")[getattr(tree, "TTStubs_clusId1")[j]]
                self.assertEqual(geoId01, geoId02)
                self.assertEqual(geoId11, geoId12)
                self.assertEqual(stackId, stackId0)
                self.assertEqual(stackId, stackId1)

                # test associated tp
                if self.bvec1[j]:
                    simPt1 = getattr(tree, "TTStubs_simPt")[j]
                    simPt2 = getattr(tree, "trkParts_pt")[getattr(tree, "TTStubs_tpId")[j]]
                    self.assertEqual(simPt1, simPt2)
                else:
                    simPt1 = getattr(tree, "TTStubs_simPt")[j]
                    self.assertEqual(simPt1, -99.)
            with self.assertRaises(IndexError):
                print self.fvec1[n]

    #def test_TTTracks(self):
    #    tree = self.ttree
    #    tree.SetBranchAddress("TTTracks_pt", self.fvec1)
    #    tree.SetBranchAddress("TTTracks_phi", self.fvec2)
    #    tree.SetBranchAddress("TTTracks_isGenuine", self.bvec1)
    #    for i in xrange(self.nevents):
    #        tree.GetEntry(i)
    #        n = getattr(tree, "TTTracks_size")
    #        self.assertTrue(n >= 0)  # can be zero
    #        for j in xrange(n):
    #            self.assertTrue(self.fvec1[j] > 0)
    #            self.assertTrue(abs(self.fvec2[j]) < pi+1e-6)

    #            # test associated tp
    #            if self.bvec1[j]:
    #                simPt1 = getattr(tree, "TTTracks_simPt")[j]
    #                simPt2 = getattr(tree, "trkParts_pt")[getattr(tree, "TTTracks_tpId")[j]]
    #                self.assertEqual(simPt1, simPt2)
    #            else:
    #                simPt1 = getattr(tree, "TTTracks_simPt")[j]
    #                self.assertEqual(simPt1, -99.)
    #        with self.assertRaises(IndexError):
    #            print self.fvec1[n]

    def test_simPixelDigis(self):
        tree = self.ttree
        tree.SetBranchAddress("simPixelDigis_x", self.fvec1)
        tree.SetBranchAddress("simPixelDigis_y", self.fvec2)
        tree.SetBranchAddress("simPixelDigis_z", self.fvec3)
        tree.SetBranchAddress("simPixelDigis_modId", self.uvec1)
        tree.SetBranchAddress("simPixelDigis_chan", self.ivec1)
        for i in xrange(self.nevents):
            tree.GetEntry(i)
            n = getattr(tree, "simPixelDigis_size")
            self.assertTrue(n > 0)
            for j in xrange(n): 
                self.assertTrue(abs(self.fvec1[j]) < 110)
                self.assertTrue(abs(self.fvec2[j]) < 110)
                self.assertTrue(abs(self.fvec3[j]) < 300)
                self.assertTrue(self.uvec1[j] < 230000)
                self.assertTrue(self.ivec1[j] < ((1023<<11)|31)+1)
            with self.assertRaises(IndexError):
                print self.fvec1[n]

    #@unittest.skipIf(process != "particleGun", "Not particle gun")
    def test_particleGun(self):
        tree = self.ttree
        for i in xrange(self.nevents):
            tree.GetEntry(i)
            self.assertEqual(getattr(tree, "genParts_size"), 1)
            self.assertEqual(abs(getattr(tree, "genParts_pdgId")[0]), 13)

            nsimtracks = 0
            for j in xrange(getattr(tree, "simTracks_size")):
                if abs(getattr(tree, "simTracks_pdgId")[j]) == 13:
                    nsimtracks += 1
            self.assertEqual(nsimtracks, 1)

            #nsimvertices = 0
            #for j in xrange(getattr(tree, "simVertices_size")):
            #    if abs(getattr(tree, "simVertices_vz")[j]) < 20:
            #        nsimvertices += 1
            #self.assertEqual(nsimvertices, 1)

            ntrkparts = 0
            for j in xrange(getattr(tree, "trkParts_size")):
                if abs(getattr(tree, "trkParts_pdgId")[j]) == 13:
                    ntrkparts += 1
            self.assertEqual(ntrkparts, 1)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        TestNTuple.infile = sys.argv.pop()

    unittest.main()

