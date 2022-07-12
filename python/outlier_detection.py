#! /usr/bin/python

import parsers as tparse
import time_series as tsp
import datetime
import numpy as np

def fractional_days(startt, endt):
    difference = endt - startt
    return difference.total_seconds() / datetime.timedelta(days=1).total_seconds()

def threeSigma(ts, ct=tsp.CoordinateType.Ellipsoidal, delete=False):

    num_pts = 15

    for ctype in tsp.coordinate_type_keys(ct):
        if not ts.includes_key(ctype):
            errmsg = 'ERROR. Coordinate type missing!'
            raise RuntimeERror(errmsg)

    more_outliers = 1
    iteration = 1
    max_iterations = 10
    while more_outliers and iteration < max_iterations:
        more_outliers = 0
        iteration += 1

        for count, entry in enumerate(ts._data):
           middle_idx = count
           start_idx = 0 if (middle_idx < num_pts) else middle_idx - num_pts
           stop_idx  = ts.size() if (middle_idx + num_pts > ts.size()) else middle_idx + num_pts

           if start_idx == 0:
               while stop_idx < ts.size() and stop_idx - start_idx < 2*num_pts:
                   stop_idx+=1

           if stop_idx == ts.size():
               while start_idx >=0 and stop_idx - start_idx < 2*num_pts:
                   start_idx -= 1

           t0 = ts._data[count]['t']
           # tarray = [ fractional_days(t0, _t['t']) for _t in ts._data[start_idx:stop_idx] ]
           for component in tsp.coordinate_type_keys(ct):
               val = ts._data[count][component]
               yarray = [ _t[component] for _t in ts._data[start_idx:stop_idx] ]
               mean = np.mean(yarray)
               stddev = np.std(yarray)
               if abs(mean-val) > 3e0 * stddev:
                   print('Marking epoch {:} as outlier; component {:}'.format(t0, component))
                   ts._data[count]['ignore'] = True
                   more_outliers += 1

    if delete:
        print('Delteting outliers from time-series ...')
        return ts.drop_ignored()
    else:
        return ts
