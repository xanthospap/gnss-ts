#! /usr/bin/python

import numpy    as np
import datetime as dt

def make_model(t, earthqk_at):
    a = 0e0
    v = 0e0
    gaus_mean = 0e0
    gaus_std  = 0.05

    y = np.zeros(t.size)
    np.random.normal(loc=gaus_mean, scale=gaus_std, size=t.size)
    i = 0
    for mjd in np.nditer(t):
        tmp  = a
        tmp += v*(mjd - 
        if mjd >= earthqk_at:

start_mjd = 53005e0
stop_mjd  = 55562e0
t         = np.arange(start_mjd, stop_mjd, 1e0)
