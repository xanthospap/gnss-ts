#! /usr/bin/python

import parsers as tparse
import time_series as tsp
import outlier_detection as tod
import model_fit as tm
import os
import sys
import datetime
import argparse

parser = argparse.ArgumentParser(
    description='Simple plot utility for visualizing time-series files'
)
parser.add_argument('-i', '--cts-file',
                    action='store',
                    required=True,
                    help='Input .cts file',
                    metavar='INPUT_FILE',
                    dest='fcts'
)
parser.add_argument('-r', '--cts-urapid-file',
                    action='store',
                    required=True,
                    help='Input .cts file',
                    metavar='INPUT_FILE',
                    dest='ucts'
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

args = parser.parse_args()

if not os.path.isfile(args.fcts):
    print('ERROR Failed to locate cts file {:}'.format(args.cts))
    sys.exit(1)
if not os.path.isfile(args.ucts):
    print('ERROR Failed to locate cts file {:}'.format(args.cts))
    sys.exit(1)

try:
    args.tstart < datetime.datetime.max
except:
    tstart = datetime.datetime.strptime(args.tstart, "%Y-%m-%d")
    args.tstart = tstart

fts = tparse.parse_cts(args.fcts, start=args.tstart)
uts = tparse.parse_cts(args.ucts, start=args.tstart)

allts = fts.append_if_missing(uts).sort()

allts.dump_as_cts()
