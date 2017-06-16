#! /usr/bin/python

import sys
import datetime

import gts.readers
import gts.timeseries as ts
import gts.tsploters

## first cmd is station cts file
try:
    cts_file = sys.argv[1]
except:
    print >> sys.stderr, '[ERROR] Usage: plotcts.py <cts file>'
    sys.exit(1)

## read from dso cts format
cts = gts.readers.read_cts(sys.argv[1])

## trasnform to topocentric
cts = cts.transform(ts.CoordinateType.Topocentric)

## try a simple, linear fit
cts.dummy_lin_fit()

## plot the topocentric time-series
gts.tsploters.ts_plot(cts, y_erbar=True)

sys.exit(0)
