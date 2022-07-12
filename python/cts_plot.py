#! /usr/bin/python

import parsers as tparse
import time_series as tsp
import outlier_detection as tod
import os
import sys
import datetime
import argparse
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
parser.add_argument('-s', '--site-name',
                    action='store',
                    required=False,
                    help='Site name to be added to the plot title',
                    metavar='SITE_NAME',
                    dest='site_name',
                    default=''
)
parser.add_argument('--no-outliers',
                    dest='remove_outliers',
                    help='Find and remove outliers before ploting.',
                    action='store_true'
)

args = parser.parse_args()

if not os.path.isfile(args.cts):
    print('ERROR Failed to locate cts file {:}'.format(args.cts))
    sys.exit(1)

ts = tparse.parse_cts(args.cts)
ts = ts.topocentric().drop_coordinate_type(
    tsp.CoordinateType.Cartesian).drop_coordinate_type(tsp.CoordinateType.Ellipsoidal)

if args.remove_outliers:
    ts = tod.threeSigma(ts, tsp.CoordinateType.Topocentric, True)

fig = plt.figure()
gs = fig.add_gridspec(3, hspace=0)
axs = gs.subplots(sharex=True, sharey=True)
fig.suptitle('Time-Series for Station ' + args.site_name)

axs[0].scatter(ts.get('t'), ts.get('east'), s=10, facecolors='black', edgecolors='r')
axs[0].set(ylabel='East [m]')
axs[1].scatter(ts.get('t'), ts.get('north'), s=10, facecolors='black', edgecolors='r')
axs[1].set(ylabel='North [m]')
axs[2].scatter(ts.get('t'), ts.get('up'), s=10, facecolors='black', edgecolors='r')
axs[2].set(ylabel='Up [m]')

# Hide x labels and tick labels for all but bottom plot.
for ax in axs: ax.label_outer()

# Do we need to save the figure?
if args.saveas is not None:
    plt.savefig(args.saveas, transparent=args.transparent)

if not args.quiet: plt.show()
