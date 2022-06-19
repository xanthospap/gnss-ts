#! /usr/bin/python

import parsers as tparse
import time_series as tsp
import os, sys, datetime

if len(sys.argv) < 2:
    print('ERROR. Need to provide a CTS filename')
    sys.exit(1)

ts = tparse.parse_cts(sys.argv[1])
print('Size of TS = {:}'.format(ts.size()))

ts = ts.drop_coordinate_type(tsp.CoordinateType.Ellipsoidal)
ts = ts.topocentric().drop_coordinate_type(tsp.CoordinateType.Cartesian)
# ts = ts.drop_coordinate_type(tsp.CoordinateType.Cartesian)

# ts.print()
print(ts.dump())