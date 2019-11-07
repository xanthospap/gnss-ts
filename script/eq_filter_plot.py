#! /usr/bin/python2.7
#-*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import argparse

def fun(c1, c2, x):
  return c1 + c2*np.log10(x)

distances = np.arange(1e0, 500e0, 10e0)
# c2s = np.arange(3e0, 6e0, 5e-1)
C1  = -5.0E0
C2  = '3e0/6e0/5e-1'

parser = argparse.ArgumentParser(
  description='Plot funtion: c1+c2*log10(x) for various c1, c2 and x values.',)
parser.add_argument('--c1',
  default = C1,
  metavar = 'C1',
  dest    = 'c1',
  type = float,
  required = False,
  help = 'Value of c1 coefficient')
parser.add_argument('--c2',
  default = C2,
  metavar = 'C2',
  dest    = 'c2',
  required = False,
  help = 'Range or value of c2 coefficient. If range use the format: from/to/step')

args  = parser.parse_args()

if args.c2.find('/') >= 0:
  c2s, c2e, c2t = [ float(i) for i in args.c2.split('/') ]
else:
  c2s = c2e = c2t = float(args.c2)
  c2e += 1e-10
c2_ar = np.arange(c2s, c2e, c2t)

for c2 in c2_ar:
  y = fun(args.c1, c2, distances);
  plt.plot(distances, y, label='c2=%3.1f'%c2)

plt.legend()
plt.title('Earthquake Filtering Formula')
plt.xlabel('Distance in Km')
plt.ylabel('%5.2f+c2*log10(distance_km)'%args.c1)
plt.grid(True)

plt.savefig('eqres.pdf', bbox_inches='tight')
