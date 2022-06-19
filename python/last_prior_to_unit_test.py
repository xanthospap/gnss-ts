#! /usr/bin/python

import parsers as tparse
import time_series as tsp
import os
import sys
import datetime

##
## TimeSeries::last_prior_to(value) returns the last element in the series,
## such that element < value, or TimeSeries::last if no such element is found.
##

if len(sys.argv) < 2:
    print('ERROR. Need to provide a CTS filename')
    sys.exit(1)

ts = tparse.parse_cts(sys.argv[1])
print('Size of TS = {:}'.format(ts.size()))

mint, maxt = ts.time_span()

# No element in the time-series is less that t!
t = mint - datetime.timedelta(days=1)
idx, tr = ts.last_prior_to(t)
assert(idx == tsp.TimeSeries.last and tr is None)

# No element in the time-series is less that t!
t = mint
idx, tr = ts.last_prior_to(t)
assert(idx == tsp.TimeSeries.last and tr is None)

# First element is less that t
t = mint + datetime.timedelta(hours=1)
idx, tr = ts.last_prior_to(t)
assert(idx == 0 and tr == ts._data[idx]['t'])

# Some element is less that t
t = mint + datetime.timedelta(days=100)
idx, tr = ts.last_prior_to(t)
idx2, tr2 = ts.last_prior_to(t, idx+1)
assert(idx2 == tsp.TimeSeries.last and tr2 is None)
idx2, tr2 = ts.last_prior_to(t, idx)
assert(idx2 == 0 and tr2 == tr)
idx2, tr2 = ts.last_prior_to(t, idx-1)
assert(idx2 == 1 and tr2 == tr)

#
t = maxt
idx, tr = ts.last_prior_to(t)
assert(idx == ts.size()-2)
t = maxt + datetime.timedelta(hours=1)
idx, tr = ts.last_prior_to(t)
assert(idx == ts.size()-1)
