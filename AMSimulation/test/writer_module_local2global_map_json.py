#!/usr/bin/env python

import json
from math import pi, sqrt, atan2

# CUIDADO: Endcap modules are not validated!

# ______________________________________________________________________________
# Functions

convert_key_to_int = lambda pairs: dict([(int(k),v) for (k,v) in pairs])

quadsum = lambda x,y: sqrt(x*x + y*y)

average = lambda x: sum(x, 0.0) / len(x)


# ______________________________________________________________________________
# Load existing maps
#ttmap = json.load(open("../data/trigger_sector_map.json"), object_pairs_hook=convert_key_to_int)
vertexmap = json.load(open("../data/module_vertices.json"), object_pairs_hook=convert_key_to_int)

def decodeLayer(moduleId):
    return int(moduleId / 10000)

def decodeLadder(moduleId):
    return int(moduleId / 100) % 100

def decodeModule(moduleId):
    return int(moduleId) % 100

def isPSModule(moduleId):
    lay = decodeLayer(moduleId)
    if 5 <= lay <= 7:
        return True
    lad = decodeLadder(moduleId)
    if 11 <= lay <= 22 and lad <= 8:
        return True
    return False

def isBarrelModule(moduleId):
    lay = decodeLayer(moduleId)
    if 5 <= lay <= 10:
        return True
    return False

def get_nrows_ncols(moduleId):
    if isPSModule(moduleId):
        nrows = 960 * 2  # half-strip unit
        ncols = 32
    else:
        nrows = 1016 * 2  # half-strip unit
        ncols = 2
    return (nrows, ncols)

class LinearRegression:
    def __init__(self):
        self.Sx = 0.
        self.Sy = 0.
        self.Sxx = 0.
        self.Sxy = 0.
        self.Syy = 0.
        self.n = 0

    def fill(self, x, y):
        self.Sx += x
        self.Sy += y
        self.Sxx += x * x
        self.Sxy += x * y
        self.Syy += y * y
        self.n += 1

    def compute(self):
        alpha = (self.Sy * self.Sxx - self.Sx * self.Sxy) / (self.n * self.Sxx - self.Sx * self.Sx)
        beta = (self.n * self.Sxy - self.Sx * self.Sy) / (self.n * self.Sxx - self.Sx * self.Sx)
        return (alpha, beta)

def get_phi_conversion(moduleId, chipId, xyz):
    (nrows, ncols) = get_nrows_ncols(moduleId)
    n = 256  # number of strips per chip

    x0 = average([xyz[0+0], xyz[9+0]])
    y0 = average([xyz[0+1], xyz[9+1]])
    x1 = average([xyz[3+0], xyz[6+0]])
    y1 = average([xyz[3+1], xyz[6+1]])

    regression = LinearRegression()
    for i in xrange(n):
        j = (chipId * n) + i
        x = x0 + (x1 - x0)/nrows * (j+0.5)
        y = y0 + (y1 - y0)/nrows * (j+0.5)
        phi = atan2(y, x)
        regression.fill(i, phi)
    (alpha, beta) = regression.compute()
    return (alpha, beta)

def get_r_conversion(moduleId, chipId, xyz):
    (nrows, ncols) = get_nrows_ncols(moduleId)

    if isBarrelModule(moduleId):
        n = 256  # number of strips per chip

        x0 = average([xyz[0+0], xyz[9+0]])
        y0 = average([xyz[0+1], xyz[9+1]])
        x1 = average([xyz[3+0], xyz[6+0]])
        y1 = average([xyz[3+1], xyz[6+1]])
    else:
        chipId = 0  # disregard chipId
        nrows = ncols
        n = ncols

        x0 = average([xyz[0+0], xyz[3+0]])
        y0 = average([xyz[0+1], xyz[3+1]])
        x1 = average([xyz[6+0], xyz[9+0]])
        y1 = average([xyz[6+1], xyz[9+1]])

    regression = LinearRegression()
    for i in xrange(n):
        j = (chipId * n) + i
        x = x0 + (x1 - x0)/nrows * (j+0.5)
        y = y0 + (y1 - y0)/nrows * (j+0.5)
        r = quadsum(x, y)
        regression.fill(i, r)
    (alpha, beta) = regression.compute()
    return (alpha, beta)

def get_z_conversion(moduleId, chipId, xyz):
    (nrows, ncols) = get_nrows_ncols(moduleId)
    chipId = 0  # disregard chipId
    n = ncols

    z0 = average([xyz[0+2], xyz[3+2]])
    z1 = average([xyz[6+2], xyz[9+2]])

    regression = LinearRegression()
    for i in xrange(n):
        j = (chipId * n) + i
        z = z0 + (z1 - z0)/ncols * (j+0.5)
        regression.fill(i, z)
    (alpha, beta) = regression.compute()
    return (alpha, beta)

def get_conversions():
    conversions = {}
    for moduleId, xyz in vertexmap.iteritems():
        if moduleId < 0:  continue

        m = 8  # number of chips
        for chipId in xrange(m):
            (x_phi0, x_phi) = get_phi_conversion(moduleId, chipId, xyz)
            (x_z0, x_z) = get_z_conversion(moduleId, chipId, xyz)
            (x_r0, x_r) = get_r_conversion(moduleId, chipId, xyz)

            conversions[(moduleId, chipId)] = (x_phi0, x_phi, x_z0, x_z, x_r0, x_r)
    return conversions

# ______________________________________________________________________________
# Main
writeme = []
conversions = get_conversions()

for (key, values) in sorted(conversions.iteritems()):
    (moduleId, chipId) = key
    (x_phi0, x_phi, x_z0, x_z, x_r0, x_r) = values
    s = "%i,%i,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e" % (moduleId, chipId, x_phi0, x_phi, x_z0, x_z, x_r0, x_r)
    writeme.append(s)

# ______________________________________________________________________________
# Write .csv file
with open("../data/module_local2global_map.csv", "w") as f:
    writeme = ["moduleId/I, chipId/I, x_phi0/D, x_phi/D, x_z0/D, x_z/D, x_r0/D, x_r/D"] + writeme
    f.write("\n".join(writeme))

# Write .json file
mymap = {}
with open("../data/module_local2global_map.csv", "r") as f:
    for line in f:
        if not line[0].isdigit():
            continue
        values = line.split(",")
        assert(len(values) == 8)

        # Convert to int or float
        #values = [round(float(x),12) if "." in x else int(x) for x in values]
        values = [float(x) if "e" in x else int(x) for x in values]

        #key = (values[0], values[1])
        key = values[0]*100 + values[1]
        values = values[2:]
        mymap[key] = values

json.dump(mymap, open("../data/module_local2global_map.json", "w"), sort_keys=True)

