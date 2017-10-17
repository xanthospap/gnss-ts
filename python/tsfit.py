#! /usr/bin/python

import argparse
import gts.readers
import gts.timeseries as ts
import gts.tsploters
# import gts.tsploters

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

##  read in the time-series file
ots = gts.readers.read_cts(args.ts_file, args.reccom)

##  Transform time-series to topocentric (from cartesian)
ots = ots.transform(ts.CoordinateType.Topocentric)

## Construct the model to be applied
model = ts.Model()
time_span_in_years = ots.time_span().days / 365.25e0;
if time_span_in_years >= 1.5e0:
    model.add_periods(365.25/2e0)
if time_span_in_years >= 2.0e0:
    model.add_periods(365.25)

##  Fit the model to the data (per component)
xmdl, _, _ = ots.fit_model(0, model)
ymdl, _, _ = ots.fit_model(1, model)
zmdl, _, _ = ots.fit_model(2, model)

## Plot the timeseries
gts.tsploters.ts_plot(ots, False, xmdl, ymdl, zmdl)
