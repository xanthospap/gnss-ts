#! /bin/python

import sys, math, datetime
import numpy as np
from scipy.signal import lombscargle as ls_scipy
import matplotlib.pyplot as plt

def read_new_input(filename):
    """ Read xxxx.res station's residuals file and return lists:
        epochs, res_north, res_east, res_up
    """
    epochs    = []
    res_north = [];
    res_east  = [];
    res_up    = [];
    index = 0

    with open(filename, 'r') as fin:
         for line in fin:
             if "o" not in line:
                 l = line.split()
                 epochs.append(float(l[0]))
                 res_north.append(float(l[1]))
                 res_east.append(float(l[3]))
                 res_up.append(float(l[1]))
    return epochs, res_north, res_east, res_up

t, n, e, u = read_new_input(sys.argv[1])
t = np.array(t)
n = np.array(n)

#Choose a period grid
periods = np.linspace(1, 365.25*2, 5000)
ang_freqs = 2* np.pi / periods

#Compute the unormalized periodogram
#note pre-centering of y values!
power = ls_scipy(t, n, ang_freqs)

#normalize the power
N = len(t)
power *= 2 / (N*n.std()**2)

#plot the results
fig, ax = plt.subplots()
ax.plot(periods, power)
ax.set(ylim=(0,0.15), xlabel='period(days)', ylabel='Lomb-scargle Power');

plt.show()
