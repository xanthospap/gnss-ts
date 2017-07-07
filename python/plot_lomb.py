#! /usr/bin/python2.7

from math import pow, exp, expm1, log
import matplotlib.pyplot as plt
import argparse

numfreq = 0
ofac    = 4
sig_lev = [.5, .1, .05, .01, .005, .001]
#USE_YLOG_SCALE  = True
#PLOT_SIG_LEVELS = True
#PLOT_PERIOD     = False

##############################################################################
parser = argparse.ArgumentParser(
    description='Plot Lomb-Scargle periodogram using a file with cols: \'freq power\''
    )
parser.add_argument('-i', '--input-file',
    action   = 'store',
    required = True,
    help     = 'The input file',
    metavar  = 'INPUT_FILE',
    dest     = 'data_file'
    )
parser.add_argument('-f', '--ofac',
    action   = 'store',
    required = False,
    help     = 'The ofac parameter of the periodogram (for computing significance levels)',
    metavar  = 'OFAC',
    dest     = 'ofac',
    default  = 4,
    type     = float
    )
parser.add_argument('-l', '--log-yaxis',
    action   = 'store_true',
    required = False,
    help     = 'Set the y axis to logarithmic scale',
    dest     = 'use_ylog_scale',
    default  = False
    )
parser.add_argument('-p', '--use-periods',
    action   = 'store_true',
    required = False,
    help     = 'Set the x axis to period (not frequency)',
    dest     = 'plot_period',
    default  = False
    )
parser.add_argument('-s', '--show-sig-level',
    action   = 'store_true',
    required = False,
    help     = 'Plot horizontal lines for significance levels',
    dest     = 'plot_sig_levels',
    default  = False
    )
##############################################################################
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

def significance(z, M):
    try:
        return 1e0 - pow(-expm1(z), M)
    except:
        print 'Overflow for z={}'.format(z)
        return float('Inf')
##############################################################################

args = parser.parse_args()

px, py   = read_pxy(args.data_file)
numfreq  = len(px)
sig_prob = [ -log(pow(p,1e0/(2e0*numfreq/ofac))) for p in sig_lev ]

fig = plt.figure()
ax1 = fig.add_subplot(211)

ax1.set_title("Lomb-Scargle Periodogram")
ax1.set_xlabel('Frequency')
ax1.set_ylabel('Power')

if args.use_ylog_scale: ax1.set_yscale('log')

ax1.plot(px, py, c='r', label='the data')
if args.plot_sig_levels:
    for idx, val in enumerate(sig_prob):
        print 'max={}'.format(max(px))
        ax1.axhline(y=val, xmin=min(px), xmax=max(px), linestyle='dashed', label='{}'.format(sig_lev[idx]))
        print('ploting sign level for z = {} at {}'.format(sig_lev[idx], val))

leg = ax1.legend()

if args.plot_sig_levels:
    ## M    = 2.0 * nout / ofac
    M = 2e0 * numfreq / args.ofac
    probs = [ significance(z, M) for z in py ]
    ax2 = fig.add_subplot(212)
    ax2.plot(px, probs, c='b')

plt.show()
