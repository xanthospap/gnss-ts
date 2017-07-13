#! /usr/bin/python

import math
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

start = 53005e0
stop  = 55562e0
event = 54693e0
# event = start
t0    = start + (stop-start)/2e0

## Model is: y = x0 + Vx*(t - t0) A * log(1 + (t - t_eq) / tau)
## tmjd <- t
## x0   <- t0
## vx   <- Vx
## a1   <- A1
## t1   <- tau
## Note that detla-times should always be in decimal years
def value_at(tmjd, x0, vx, a1, t1):
    d    = x0 + vx * ( (tmjd - t0)/365.25e0 )
    D    = d
    dtq1 = (tmjd - event) / 365.25e0
    dtq  = dtq1[np.where(dtq1 >= 0)[0]]
    arg  = dtq / t1
    d[tmjd.size-arg.size :] = d[tmjd.size-arg.size :] + (a1 * np.log(1.0e0 + arg))
    return d

## Model is: y = x0 + Vx*(t - t0) A * log(1 + (t - t_eq) / tau
## tmjd <- t
## a1   <- A1 (approx)
## t1   <- tau (approx)
## Note that detla-times should always be in decimal years
## dx0 = dy/d(x0) i.e. partial derivative
## ...
def partials_at(tmjd, a1, t1):
    A      = np.zeros( (tmjd.size, 4) )
    A[:,0] = A[:,0] + 1e0
    A[:,1] = (tmjd - t0) / 365.25e0
    dtq1   = (tmjd - event) / 365.25e0
    dtq    = dtq1[np.where(dtq1 >= 0)[0]]
    arg    = dtq / t1
    A[tmjd.size-arg.size :, 2] = np.log(1e0 + arg)
    A[tmjd.size-arg.size :, 3] = a1*(-arg/t1)/(1e0+arg)
    arg2   = dtq1 / t1
    return A

## Model coefficients
x0 = 0e0;    x0_aprx = 0e0;
vx = 5e-4;   vx_aprx = 0e0;
a1 = 1e-3;   a1_aprx = 0e0;
t1 = 1.3e0;  t1_aprx = .5e0;

## observation times
t = np.arange(start, stop, 1e0)
## white noise
noise = np.random.normal(0, 1e-3, t.size)
## observation values
l = value_at(t, x0, vx, a1, t1) + noise
## observed - computed
b = l - value_at(t, x0_aprx, vx_aprx, a1_aprx, t1_aprx)
## jacobian matrix
A = partials_at(t, a1_aprx, t1_aprx)
## plot the data
plt.plot(t, l, 'b-', label='data')
A=A/1000e0
b=b/1000e0

# print matrices
np.set_printoptions(precision=10, threshold=1e10, linewidth=125)
with open('design-mat.py', 'w') as fout:
    print >> fout, A
np.set_printoptions(precision=10, threshold=1e10, linewidth=25)
with open('obs-mat.py', 'w') as fout:
    print >> fout, np.transpose(b-noise)

max_iters = 3
colors = ['b--', 'r--', 'g--', 'k--', 'm--', 'y--']
for i in range(max_iters):
    print('Iteration: {}, Model Parameters:'.format(i))
    print('X0 {:+7.4f} | {:+7.4f} | {:+7.4f}'.format(x0_aprx, x0, abs(x0_aprx-x0)))
    print('Vx {:+7.4f} | {:+7.4f} | {:+7.4f}'.format(vx_aprx, vx, abs(vx_aprx-vx)))
    print('A  {:+7.4f} | {:+7.4f} | {:+7.4f}'.format(a1_aprx, a1, abs(a1_aprx-a1)))
    print('T  {:+7.4f} | {:+7.4f} | {:+7.4f}'.format(t1_aprx, t1, abs(t1_aprx-t1)))
    plt.plot(t, value_at(t, x0_aprx, vx_aprx, a1_aprx, t1_aprx), colors[i], label='iter: {}'.format(i), linewidth=2.0)
    x, _, _, _, = np.linalg.lstsq(A, b)
    x0_aprx += x[0]
    vx_aprx += x[1]
    a1_aprx += x[2]
    t1_aprx += x[3]
    b = l - value_at(t, x0_aprx, vx_aprx, a1_aprx, t1_aprx)
    A = partials_at(t, a1_aprx, t1_aprx)

plt.legend()
plt.show()
