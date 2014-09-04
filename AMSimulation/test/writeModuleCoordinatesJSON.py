#!/usr/bin/env python

import json

#class Coord:
#    def __init__(self, z, r, eta, phi, sensorSpacing):
#        self.z = z
#        self.r = r
#        self.eta = eta
#        self.phi = phi
#        self.sensorSpacing = sensorSpacing

mymap = {}
with open("../data/module_coordinates.csv", "r") as f:
    for line in f:
        if not line[0].isdigit():
            continue
        values = line.split(",")
        assert(len(values) == 6)

        # Convert to int or float
        values = [float(x) if "." in x else int(x) for x in values]

        key = values[0]
        values = values[1:]
        mymap[key] = values

assert(len(mymap) == 15428)

json.dump(mymap, open("../data/module_coordinates.json", "w"), sort_keys=True)

