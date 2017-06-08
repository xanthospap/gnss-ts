#! /bin/python

from math import pow, exp, log
import matplotlib.pyplot as plt

numfreq = 0
ofac    = 4
sig_lev = [.5, .1, .05, .01, .005, .001]
USE_YLOG_SCALE  = True
PLOT_SIG_LEVELS = True
PLOT_PERIOD     = False

def read_pxy(filename):
    px = []
    py = []
    with open(filename) as fn:
        for line in fn.readlines():
            if len(line) > 1:
                l = map(float, line.split())
                px.append(l[0])
                py.append(l[1])
    return px, py

px, py   = read_pxy('lomb.out')
numfreq  = len(px)
# sig_prob = [ 1-pow(1-exp(-z), 2*numfreq/ofac) for z in sig_lev ]
sig_prob = [ -log(pow(p,1e0/(2e0*numfreq/ofac))) for p in sig_lev ]

fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.set_title("Lomb-Scargle Periodogram")
ax1.set_xlabel('Frequency')
ax1.set_ylabel('Power')

if USE_YLOG_SCALE: ax1.set_yscale('log')

ax1.plot(px, py, c='r', label='the data')
if PLOT_SIG_LEVELS:
    for idx, val in enumerate(sig_prob):
        print 'max={}'.format(max(px))
        ax1.axhline(y=val, xmin=min(px), xmax=max(px), linestyle='dashed', label='{}'.format(sig_lev[idx]))
        print('ploting sign level for z = {} at {}'.format(sig_lev[idx], val))

leg = ax1.legend()

plt.show()
