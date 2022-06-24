#! /usr/bin/python

import datetime
import argparse
import os
import parsers as tp
import time_series as ts
import model_fit as tm
import periodogram as th
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(
    description='Simple plot utility for visualizing time-series files'
)
parser.add_argument('-i', '--cts-file',
                    action='store',
                    required=True,
                    help='Input .cts file',
                    metavar='INPUT_FILE',
                    dest='cts'
)
parser.add_argument('-f', '--from',
                    action='store',
                    required=False,
                    help='Only consider data later or equal to this datetime',
                    metavar='START_EPOCH',
                    dest='tstart',
                    default=datetime.datetime.min
)
parser.add_argument('-t', '--to',
                    action='store',
                    required=False,
                    help='Only consider data prior to this datetime',
                    metavar='STOP_EPOCH',
                    dest='tstop',
                    default=datetime.datetime.max
)
parser.add_argument('-o', '--output-file',
                    action='store',
                    required=False,
                    help='Save the time-series plot to an output file. The format is inferred from the extension of the given filename.',
                    metavar='OUTPUT_FILE',
                    dest='saveas',
                    default=None
)
parser.add_argument('--transparent',
                    dest='transparent',
                    help='Only relevant if we are saving the figure to an output file. If set, the figure will be transparent (see https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.savefig.html)',
                    action='store_true'
)
parser.add_argument('-q', '--quiet',
                    dest='quiet',
                    help='Do not display the produced figure.',
                    action='store_true'
)

args = parser.parse_args()

if not os.path.isfile(args.cts):
    print('ERROR Failed to locate cts file {:}'.format(args.cts))
    sys.exit(1)

s = tp.parse_cts(args.cts)
s = s.topocentric().drop_coordinate_type(
    ts.CoordinateType.Cartesian).drop_coordinate_type(ts.CoordinateType.Ellipsoidal)

## simple, linear model
lmodels = []
for ct in ts.coordinate_type_keys(ts.CoordinateType.Topocentric):
    mdl = tm.TsModel()
    lmodels.append(mdl.fit(s.get('t'), s.get(ct))[0])

## de-trended    
residuals = s.residuals(model=lmodels, coordinate_type=ts.CoordinateType.Topocentric)

lslist = th.lomb_scargle(residuals, coordinate_type=ts.CoordinateType.Topocentric, model=lmodels)

## plot (try north)
component = 'north'
m = lslist[0]
periods, power = m.periodogram_auto(nyquist_factor=100)
fig, ax = plt.subplots()
ax.plot(periods, power)
ax.set(xlim=(0.2, 2.0), ylim=(0, 0.99), xlabel='period (days)', ylabel='Lomb-Scargle Power')
plt.show()
