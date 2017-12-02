#! /usr/bin/python

from __future__ import print_function
import argparse

parser = argparse.ArgumentParser(
    description='Given a range (aka xmin, xmax) and a number of ticks, compute annotation interval.'
)
parser.add_argument('min',
    type     = float,
    help     = 'Min tick on axis'
)
parser.add_argument('max',
    type     = float,
    help     = 'Max tick on axis.'
)
parser.add_argument('ticks',
    type     = int,
    help     = 'Number of ticks on axis.'
)

args = parser.parse_args()

assert args.max > args.min
assert args.ticks > 0

anev = (args.max - args.min) / float(args.ticks)

if anev < 1e-3: anev = 1e-3

print('{:5.3f}'.format(anev))
