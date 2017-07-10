#! /usr/bin/python

import math
import numpy as np
from scipy.optimize import curve_fit

start = 53005e0
stop  = 55562e0
event = 54693e0
t0    = start + (stop-start)/2e0

## Model is: y = x0 + Vx*(t - t0) A * log(1 + (t - t_eq) / tau
## tmjd <- t
## x0   <- t0
## vx   <- Vx
## a1   <- A1
## t1   <- tau
## Note that detla-times should always be in decimal years
def value_at(tmjd, x0, vx, a1, t1):
    d = x0 + vx * ( (tmjd - t0)/365.25e0 )
    #if tmjd > event:
    dtq = (tmjd - event) / 365.25e0
    # d += (a1 * math.log(dtq / t1))
    tmjd.size
    d += float(tmjd > event) * ( a1 * math.exp(dtq / t1) )
    return d

## Model is: y = x0 + Vx*(t - t0) A * log(1 + (t - t_eq) / tau
## tmjd <- t
## a1   <- A1 (approx)
## t1   <- tau (approx)
## Note that detla-times should always be in decimal years
## dx0 = dy/d(x0) i.e. partial derivative
## ...
def partials_at(tmjd, a1, t1):
    dx0 = 0e0
    dvx = (tmjd - t0) / 365.25e0
    da  = 0e0
    dt  = 0e0
    if tmjd > event:
        dtq = (tmjd - event) / 365.25e0
        arg = dtq / t1
        #da  = math.log(1e0 + arg)
        #dt  = a1*(-arg/t1)/(1+arg)
        da   = math.exp(dtq / t1)
        dt   = a1 * math.exp(dtq / t1) * (-dtq) / (t1*t1)
    return dx0, dvx, da, dt

## Model coefficients
x0 = 0e0;    x0_aprx = 0e0;
vx = 5e-3;   vx_aprx = 0e0;
a1 = 1e-3;   a1_aprx = 1e-3;
t1 = 1.2e0;  t1_aprx = 1.2e0;

## observation times
t = np.arange(start, stop, 1e0)
## white noise
noise = np.random.normal(0, 1e-3, t.size)
## observation values
l = np.asarray([value_at(t[i], x0, vx, a1, t1) + noise[i] for i in range(0, t.size)])
## observed - computed
b = l - np.asarray([value_at(t[i], x0_aprx, vx_aprx, a1_aprx, t1_aprx) for i in range(0, t.size)])
## jacobian matrix
A = np.zeros( (t.size, 4) )
for row in range(t.size):
    A[row, :] = partials_at(t[row], a1_aprx, t1_aprx)

#for i in range(t.size):
#    print('{:15.5f} {:15.10f} {:15.10f}'.format(t[i], l[i], b[i]))

max_iters = 5
for i in range(max_iters):
    x, _, _, _, = np.linalg.lstsq(A, b)
    x0_aprx += x[0]
    vx_aprx += x[1]
    a1_aprx += x[2]
    t1_aprx += x[3]
    print('Iteration: {}, Model Parameters:'.format(i))
    print('X0 {:+7.4f} | {:+7.4f} | {:+7.4f}'.format(x0_aprx, x0, abs(x0_aprx-x0)))
    print('Vx {:+7.4f} | {:+7.4f} | {:+7.4f}'.format(vx_aprx, vx, abs(vx_aprx-vx)))
    print('A  {:+7.4f} | {:+7.4f} | {:+7.4f}'.format(a1_aprx, a1, abs(a1_aprx-a1)))
    print('T  {:+7.4f} | {:+7.4f} | {:+7.4f}'.format(t1_aprx, t1, abs(t1_aprx-t1)))

xx = [ t[i] for i in range (t.size) ]
yy = [ l[i] for i in range (t.size) ]
popt, pcov = curve_fit(value_at, xx, yy)
