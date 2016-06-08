#!/usr/bin/env python

import json

mymap = {}
with open("../data/module_localtoglobal_map.csv", "r") as f:
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

json.dump(mymap, open("../data/module_localtoglobal_map.json", "w"), sort_keys=True)

