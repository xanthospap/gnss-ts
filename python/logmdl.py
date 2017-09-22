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

## Model is: y = x0 + Vx*(t - t0) 
##             + A * log(1 + (t - t_eq) /tau1)
##             + B * (1 - exp(-(t-t_eq)/tau2) )
## tmjd <- t
## x0   <- t0
## vx   <- Vx
## a1   <- A1
## t1   <- tau1
## a2   <- B
## t2   <- tau2
## Note that detla-times should always be in decimal years
def value_at(tmjd, x0, vx, a1, t1, a2, t2):
    d    = x0 + vx * ( (tmjd - t0)/365.25e0 )
    D    = d
    dtq1 = (tmjd - event) / 365.25e0
    dtq  = dtq1[np.where(dtq1 >= 0)[0]]
    arg1 = dtq / t1
    arg2 = dtq / t2
    d[tmjd.size-arg1.size :] = ( d[tmjd.size-arg1.size :]
                            + (a1 * np.log(1.0e0 + arg1))
                            + (a2*(1e0-np.exp(-arg2))) )
    return d

## Model is: y = x0 + Vx*(t - t0) 
##             + A * log(1 + (t - t_eq) /tau1)
##             + B * (1 - exp(-(t-t_eq)/tau2) )
## tmjd <- t
## a1   <- A1   (approx)
## t1   <- tau1 (approx)
## a2   <- B1   (approx)
## t2   <- tau2 (approx)
## Note that detla-times should always be in decimal years
## dx0 = dy/d(x0) i.e. partial derivative
## ...
## A = | dy/dx0 dy/dVx dy/da1 dy/dt1 dy/da2 dy/dt2 |
def partials_at(tmjd, a1, t1, a2, t2):
    A      = np.zeros( (tmjd.size, 6) )
    A[:,0] = A[:,0] + 1e0
    A[:,1] = (tmjd - t0) / 365.25e0
    dtq1   = (tmjd - event) / 365.25e0
    dtq    = dtq1[np.where(dtq1 >= 0)[0]]
    arg1   = dtq / t1
    arg2   = dtq / t2
    A[tmjd.size-arg1.size :, 2] = np.log(1e0 + arg1)
    A[tmjd.size-arg1.size :, 3] = a1*(-arg1/t1)/(1e0+arg1)
    A[tmjd.size-arg1.size :, 4] = 1e0-np.exp(-arg2)
    A[tmjd.size-arg1.size :, 5] = a2*np.exp(-arg2)*dtq/(t2*t2)
    return A

def jacobian(tmjd, x0, vx, a1, t1, a2, t2):
    return partials_at(tmjd, a1, t1, a2, t2)

## Model coefficients
x0 = 0e0;    x0_aprx = 0e0;
vx = 5e-4;   vx_aprx = 0e0;
a1 = 1e-3;   a1_aprx = 0e0;
t1 = 1.3e0;  t1_aprx = .5e0;
a2 = 5e-3;   a2_aprx = 0e0;
t2 = 0.7e0;  t2_aprx = .5e0;

## observation times
t = np.arange(start, stop, 1e0)
## white noise
noise = np.random.normal(0, 1e-3, t.size)
## observation values
l = value_at(t, x0, vx, a1, t1, a2, t2) + noise
## observed - computed
b = l - value_at(t, x0_aprx, vx_aprx, a1_aprx, t1_aprx, a2_aprx, t2_aprx)
## jacobian matrix
A = partials_at(t, a1_aprx, t1_aprx, a2_aprx, t2_aprx)
## plot the data
plt.plot(t, l, 'b-', label='data')

# print matrices
"""
np.set_printoptions(precision=10, threshold=1e10, linewidth=125)
with open('design-mat.py', 'w') as fout:
    print >> fout, A
np.set_printoptions(precision=10, threshold=1e10, linewidth=25)
with open('obs-mat.py', 'w') as fout:
    print >> fout, np.transpose(b-noise)
"""

#initial guess
p0 = np.asarray([x0_aprx, vx_aprx, a1_aprx, t1_aprx, a2_aprx, t2_aprx])
x_best = np.asarray([x0, vx, a1, t1, a2, t2])
min_vals = np.asarray([-np.inf, -np.inf, -1e3, 5e0/365.25, -1e3, 5e0/365.25])
max_vals = np.asarray([+np.inf, +np.inf, +1e3, 5e0, +1e3, 5e0])
popt, pcov = curve_fit(value_at, t, b, p0, jac=jacobian, bounds=(min_vals, max_vals))
print popt

max_iters = 5
colors = ['b--', 'r--', 'g--', 'k--', 'm--', 'y--']
for i in range(max_iters):
    print('Iteration: {}, Model Parameters:'.format(i))
    print('X0 {:+7.4f} | {:+7.4f} | {:+7.4f}'.format(x0_aprx, x0, abs(x0_aprx-x0)))
    print('Vx {:+7.4f} | {:+7.4f} | {:+7.4f}'.format(vx_aprx, vx, abs(vx_aprx-vx)))
    print('A1 {:+7.4f} | {:+7.4f} | {:+7.4f}'.format(a1_aprx, a1, abs(a1_aprx-a1)))
    print('T1 {:+7.4f} | {:+7.4f} | {:+7.4f}'.format(t1_aprx, t1, abs(t1_aprx-t1)))
    print('A2 {:+7.4f} | {:+7.4f} | {:+7.4f}'.format(a2_aprx, a2, abs(a2_aprx-a2)))
    print('T2 {:+7.4f} | {:+7.4f} | {:+7.4f}'.format(t2_aprx, t2, abs(t2_aprx-t2)))
    plt.plot(t, value_at(t, x0_aprx, vx_aprx, a1_aprx, t1_aprx, a2_aprx, t2_aprx), colors[i], label='iter: {}'.format(i), linewidth=2.0)
    x, _, _, _, = np.linalg.lstsq(A, b)
    x0_aprx += x[0]
    vx_aprx += x[1]
    a1_aprx += x[2]
    t1_aprx += x[3]
    a2_aprx += x[4]
    t2_aprx += x[5]
    b = l - value_at(t, x0_aprx, vx_aprx, a1_aprx, t1_aprx, a2_aprx, t2_aprx)
    A = partials_at(t, a1_aprx, t1_aprx, a2_aprx, t2_aprx)

plt.legend()
plt.show()
