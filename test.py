#! /usr/bin/python

import sys
import datetime

import gts.readers
import gts.timeseries as ts
import gts.tsploters

## read from dso cts format
dion_ts = gts.readers.read_cts(sys.argv[1] + ".cts", 'repro16')

## trasnform to topocentric
dion_ts = dion_ts.transform( ts.CoordinateType.Topocentric )

## let's try splitting the ts
# dion_a, dion_b = dion_ts.split(datetime.datetime(2012, 05, 01))

#print 'Missing dates:'
#dion_ts.report_missing_days()

## try a simple, linear fit
#dion_ts.dummy_lin_fit()

#gts.tsploters.ts_plot( dion_ts, y_erbar=True )
dion_ts.toJson()

sys.exit(0)
