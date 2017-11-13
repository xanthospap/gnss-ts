#! /usr/bin/python
# -*- coding: utf-8 -*-

import sys
import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import argparse

def read_ts(ifile):
    mjdl = []
    xl   = []
    yl   = []
    zl   = []
    with open(ifile, 'r') as fin:
        for line in fin.readlines():
            if 'o' not in line:
                mjd, x, sx, y, sy, z, sz = [ float(k) for k in line.split() ]
                mjdl.append(mjd)
                xl.append(x)
                yl.append(y)
                zl.append(z)
    return mjdl, xl, yl, zl

def lsp(t, y, freqs):
    assert len(t) == len(y)
    return signal.lombscargle(np.asarray(t), np.asarray(y), freqs)

def angular2f(angf):
    """ ω = 2πf
    """
    return angf / (2e0*np.pi)

def f2angular(f):
    return f*2e0*np.pi

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
parser.add_argument('-m', '--min-freq',
    action   = 'store',
    required = False,
    help     = 'Minimum (angular) frequency',
    metavar  = 'MIN_FREQ',
    dest     = 'min_freq',
    type     = float,
    default  = 1e-5
    )
parser.add_argument('-x', '--max-freq',
    action   = 'store',
    required = False,
    help     = 'Maximum (angular) frequency',
    metavar  = 'MAX_FREQ',
    dest     = 'max_freq',
    type     = float,
    default  = 1e0
    )
parser.add_argument('-s', '--freq-step',
    action   = 'store',
    required = False,
    help     = '(angular) Frequency step size',
    metavar  = 'FREQ_STEP',
    dest     = 'freq_step',
    default  = None
    )
parser.add_argument('-z', '--freq-size',
    action   = 'store',
    required = False,
    help     = 'Number of frequencies',
    metavar  = 'FREQ_SIZE',
    dest     = 'freq_size',
    default  = None
    )
parser.add_argument('-n', '--normalize',
    action   = 'store_false',
    required = False,
    help     = 'Normalize periodogram values',
    dest     = 'normalize'
    )
parser.add_argument('-a', '--transform-to-angular',
    action = 'store_true',
    required  = False,
    help = 'Input frequencies are NOT angular',
    dest = 'not_angular'
    )
parser.add_argument('-p', '--non-angular-periodogram',
    action = 'store_true',
    required  = False,
    help = 'Priodogram frequencies are transformed to non-angular',
    dest = 'non_angular_periodogram'
    )

args = parser.parse_args()

if not args.freq_step and not args.freq_size:
    print '[ERROR] Must specify either \'FREQ_SIZE\' of \'FREQ_STEP\'.'
    sys.exit(1)

min_freq = args.min_freq if not args.not_angular else f2angular(args.min_freq)
max_freq = args.max_freq if not args.not_angular else f2angular(args.max_freq)

if args.freq_step:
    freq_step = float(args.freq_step)
    freqs = np.arrange(min_freq, max_freq, freq_step)
    print '[DEBUG] Angular frequency range is [{:}, {:}] with a step of {:}'.format(freqs[0], freqs[len(freqs)-1], freq_step)
else:
    freq_size = int(args.freq_size)
    freqs = np.linspace(min_freq, max_freq, num=freq_size)
    print '[DEBUG] Angular frequency range is [{:}, {:}] with a total size of {:}'.format(freqs[0], freqs[len(freqs)-1], freq_size)

mjdl, xl, yl, zl = read_ts(args.data_file)
pgram = lsp(mjdl, xl, freqs)
max_idx = np.argmax(pgram)
print 'Maximum (angular) frequency is {:}, a period of {:} days.'.format(pgram[max_idx], 2e0*np.pi/pgram[max_idx])

factor = 1e0
if args.non_angular_periodogram:
    factor = 1e0/(2e0*np.pi)
    print '[DEBUG] Note that angular frequencies will be transformed to normal.'
plt.subplot(2, 1, 1)
plt.plot(mjdl, xl, 'b+')
plt.subplot(2, 1, 2)
if args.normalize:
    normval = float(len(mjdl))
    plt.plot(factor*freqs, np.sqrt(4*(pgram/normval)))
else:
    plt.plot(factor*freqs, pgram)
plt.show()
