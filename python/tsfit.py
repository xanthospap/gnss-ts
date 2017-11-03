#! /usr/bin/python2.7

import datetime
import argparse
import os
import sys
import gts.readers
import gts.timeseries as ts
import gts.tsploters
import gts.geodesy
import numpy as np

EXIT_SUCCESS = 0
EXIT_FAILURE = 1

parser = argparse.ArgumentParser(
    description = 'This tool is used to read in a DSO formated time-series file  '
    'and clear and sort its contents.'
)

parser.add_argument('-i', '--input-file',
    action = 'store',
    required = True,
    help = 'The name of the time-series file',
    metavar = 'TIME_SERIES',
    dest = 'ts_file'
)
parser.add_argument('-c', '--comment',
    action = 'store',
    required = False,
    help = 'If specified, only records containing this string as description will '
    'be read and exported',
    metavar = 'COMMENT',
    dest = 'reccom',
    default = None
)

##  Parse command line arguments
args = parser.parse_args()

try:
##  read in the time-series file
    ots = gts.readers.read_cts(args.ts_file, args.reccom)

    xmean = ots.average(0)
    ymean = ots.average(1)
    zmean = ots.average(2)
    phi, lam, hgt = gts.geodesy.cartesian2ellipsoidal(xmean,ymean,zmean)

##  Transform time-series to topocentric (from cartesian)
    ots = ots.transform(ts.CoordinateType.Topocentric)

## Construct the model to be applied
    model = ts.Model()
    time_span_in_years = ots.time_span().days / 365.25e0;
    if time_span_in_years < 1.0e0:
        print '[DEBUG] Time-series span too small for estimation:', time_span_in_years
        sys.exit(EXIT_FAILURE)
    if time_span_in_years >= 1.5e0:
        model.add_periods(365.25/2e0)
    if time_span_in_years >= 2.0e0:
        model.add_periods(365.25)
    #model.add_offsets(datetime.datetime(2014,02,01)) #Kephalonia
    #model.add_offsets(datetime.datetime(2015,11,01)) #Kephalonia
    #model.add_offsets(datetime.datetime(2015,11,17)) #Lefkada
    #model.add_offsets(datetime.datetime(2014,07,1)) #mol2 and trp2
    model.add_offsets(datetime.datetime(2008,02,14)) # methoni, pyl?

##  Fit the model to the data (per component)
    xmdl, _, _ = ots.fit_model(0, model, True, 120, False)
    ymdl, _, _ = ots.fit_model(1, model, True, 120, False)
    zmdl, _, _ = ots.fit_model(2, model, True, 120, False)

## Plot the timeseries
    gts.tsploters.ts_plot(ots, False, True, xmdl, ymdl, zmdl)

##
    station = os.path.basename(args.ts_file)[0:4]
    print station, np.rad2deg(lam), np.rad2deg(phi), hgt, xmdl.__vx__[0]*1e3, 0.0, ymdl.__vx__[0]*1e3, 0.0, zmdl.__vx__[0]*1e3, 0.0
    status = EXIT_SUCCESS
except:
    status = EXIT_FAILURE

sys.exit(status)
