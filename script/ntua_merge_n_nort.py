#! /usr/bin/python

from __future__ import print_function
import argparse
import datetime
import collections

def strip_comment(line):
    if len(line)>196: return line[196:].strip()
    raise ValueError()

parser = argparse.ArgumentParser(
    description='Given two cts files, this program will merge their records, sort them and filter them.'
)
parser.add_argument('-i', '--input',
    help     = 'Specify one or more input cts files',
    nargs    = '+',
    required = True,
    dest     = 'ifiles'
)
parser.add_argument('-c', '--comment',
    help     = 'Only keep records with this string as comment',
    required = False,
    dest     = 'comment',
    default  = None
)

args = parser.parse_args()

records = []

##  Collect all records
for fl in args.ifiles:
    with open(fl, 'r') as fin:
        if args.comment:
            records += [ ln for ln in fin.readlines() if strip_comment(ln) == args.comment ]
        else:
            records += [ ln for ln in fin.readlines() ]

##  Sort records
##  records = sorted(records, key=lambda trec: datetime.datetime.strptime(trec[0:19], '%Y-%m-%d %H:%M:%S'))

##  Now we need to remove duplicates, based on their time stamp
dct = {}
for record in records:
    r_eph = datetime.datetime.strptime(record[0:19], '%Y-%m-%d %H:%M:%S')
    # !note! Seconds must be an integer number for the "S" identifier to work
    r_tst = datetime.datetime.strptime(record[169:169+record[169:195].rfind(".")], '%Y-%m-%d %H:%M:%S')
    if r_eph not in dct:
        dct[r_eph] = [record, r_tst]
    else:
        if dct[r_eph][1] < r_tst:
            dct[r_eph] = [record, r_tst]

##  Print results in sorted order
od = collections.OrderedDict(sorted(dct.items()))
for k, v in od.items(): print("{:}".format(v[0].rstrip()))
