#! /usr/bin/python

import sys
import gts.readers
import gts.timeseries as ts
import gts.tsploters

dion_ts = gts.readers.read_ntua_cts("dion.c.cts")
dion_ts = dion_ts.transform( ts.CoordinateType.Topocentric )
gts.tsploters.ts_plot( dion_ts, y_erbar=True )

sys.exit(0)
