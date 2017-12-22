#! /usr/bin/python

"""  Given a .cts or .xyz file (actually any file starting with a date field of
     the type '%Y-%m-%d %H:%M:%S', this program will output the records (lines)
     of the file in chronological order. If any date is found more than once,
     then an error will result and no output will be provided.
     Output is directed to STDOUT.
     Example: sortcts.py <.cts/.xyz file>
"""

from __future__ import print_function
import argparse
import sys
import operator
import datetime

parser = argparse.ArgumentParser(
    description='''Rearange any file starting with a date field of type
    \'%Y-%m-%d %H:%M:%S\' in chronological order. The input file is left as is,
    while the output is directed to STDOUT''',
    formatter_class = argparse.RawTextHelpFormatter
)
parser.add_argument('-i', '--input-file',
    action   = 'store',
    required = True,
    help     = 'The input file.',
    metavar  = 'INPUT_FILE',
    dest     = 'data_file'
)

args = parser.parse_args()
ifl = args.data_file
dct = {}

with open(ifl, 'r') as fin:
    for line in fin.readlines():
        l = line.split()
        t = datetime.datetime.strptime(l[0]+' '+l[1], '%Y-%m-%d %H:%M:%S')
        if t in dct:
            print('[ERROR] Date {:} already read!'.format(t.strftime('%Y-%m-%d %H:%M:%S')))
            sys.exit(1)
        dct[t] = line

for l in sorted(dct.items(), key=operator.itemgetter(0)):
    print('{:}'.format(l[1].strip()))
