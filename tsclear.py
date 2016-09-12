#! /usr/bin/python

import sys, os
import datetime
import argparse

import gts.readers
import gts.timeseries as ts
import gts.tsploters

EXIT_SUCCESS = 0
EXIT_FAILURE = 1

parser = argparse.ArgumentParser(
    description = 'This tool is used to read in a DSO formated time-series file  '
    'and clear and sort its contents.'
)

parser.add_argument('-s', '--station',
    action = 'store',
    required = True,
    help = 'The name of the station; the time-series file will be assumed to be: '
    '\"<station-name>.cts\"',
    metavar = 'STATION_NAME',
    dest = 'tsf'
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

parser.add_argument('-t', '--topocentric',
    action = 'store_true',
    required = False,
    help = 'Transform the time-series from (geocentric) cartesian to a topocentric '
    'reference frame',
    #metavar = 'TRANSFORM_TOPOCENTRIC',
    dest = 'totopo'
)


##  Parse command line arguments
args = parser.parse_args()

##  The time-series file
ts_file = args.tsf + '.cts'

##  read in the time-series file
ots = gts.readers.read_cts(ts_file, args.reccom)

##  trasnform to topocentric if needed
if args.totopo:
    ots = ots.transform(ts.CoordinateType.Topocentric)

##  print the time-series
ots.print_as_cts()

sys.exit(EXIT_SUCCESS)
