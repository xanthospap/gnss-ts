#! /usr/bin/python

from __future__ import print_function
import argparse
import datetime

parser = argparse.ArgumentParser(
    description='Given a range (aka xmin, xmax) and a number of ticks, compute annotation interval.'
)
parser.add_argument('min',
#    type     = float,
    help     = 'Min tick on axis'
)
parser.add_argument('max',
#    type     = float,
    help     = 'Max tick on axis.'
)
parser.add_argument('ticks',
    type     = int,
    help     = 'Number of ticks on axis.',
    nargs    = '?'
)
parser.add_argument('-t', '--time',
    action   = 'store_true',
    required = False,
    help     = 'Treat min/max values as dates.',
    dest     = 'is_date',
)

args = parser.parse_args()

if not args.is_date:
    args.min, args.max = float(args.min), float(args.max)
    assert args.max > args.min
    assert args.ticks > 0

    anev = (args.max - args.min) / float(args.ticks)

    if anev < 1e-3: anev = 1e-3

    print('{:5.3f}'.format(anev))
else:
    mind   = datetime.datetime.strptime(args.min, '%Y-%b-%dT%H:%M:%S')
    maxd   = datetime.datetime.strptime(args.max, '%Y-%b-%dT%H:%M:%S')
    assert maxd > mind
    dt     = maxd - mind
    dyears = (dt.days+dt.seconds/86400)/365.25
    if dyears > 15:
        ant = 0
    elif dyears > 10:
        ant = 6
    elif dyears > 7:
        ant = 3
    elif dyears > 4:
        ant = 2
    else:
        ant = 1
    print('{:02d}'.format(ant))
        
