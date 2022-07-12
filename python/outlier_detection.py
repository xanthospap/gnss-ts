#! /usr/bin/python

import parsers as tparse
import time_series as tsp
import datetime
import numpy as np

def fractional_days(startt, endt):
    difference = endt - startt
    return difference.total_seconds() / datetime.timedelta(days=1).total_seconds()

def threeSigma(ts, **kwargs):
    if 'coordinate_key' not in kwargs and 'coordinate_type' not in kwargs:
        errmsg = 'ERROR. Must provide a \'coordinate_[type|key]\' to compute residuals'
        raise RuntimeError(errmsg)

    if 'coordinate_key' in kwargs and 'coordinate_type' in kwargs:
        errmsg = 'ERROR. Must provide one of \'coordinate_[type|key]\' to compute residuals'
        raise RuntimeError(errmsg)

    if 'coordinate_type' in kwargs:
        ccomponents = tsp.coordinate_type_keys(kwargs['coordinate_type'])
    else:
        ccomponents = [kwargs['coordinate_key']]
    

    delete  = True if 'delete' in kwargs and kwargs['delete'] == True else False
    verbose = True if 'verose' in kwargs and kwargs['verbose'] == True else False

    num_pts = 15
    max_percent = 6e0
    coef = 3e0

    more_outliers = 1
    iteration = 0
    max_iterations = 5
    while more_outliers and iteration < max_iterations:
        more_outliers = 0
        iteration += 1
        ## size of ts, excluding outliers
        starting_size = ts.size(False)
        ## marked indexes of outliers for this iteration
        outlier_indexes = []

        for count, entry in enumerate(ts._data):
           if ts._data[count]['ignore'] == False:
               middle_idx = count
               start_idx = 0 if (middle_idx < num_pts) else middle_idx - num_pts
               stop_idx  = ts.size() if (middle_idx + num_pts > ts.size()) else middle_idx + num_pts

               if start_idx == 0:
                   while stop_idx < ts.size() and stop_idx - start_idx < 2*num_pts:
                       stop_idx+=1

               if stop_idx == ts.size():
                   while start_idx >=0 and stop_idx - start_idx < 2*num_pts:
                       start_idx -= 1

               ## keep indexes of outliers in the outlier_indexes list
               t0 = ts._data[count]['t']
               # tarray = [ fractional_days(t0, _t['t']) for _t in ts._data[start_idx:stop_idx] ]
               for component in ccomponents:
                   val = ts._data[count][component]
                   # yarray = [ _t[component] for _t in ts._data[start_idx:stop_idx] if _t['ignore']==False ]
                   yarray = [ _t[component] for _t in ts._data[start_idx:stop_idx] ]
                   mean = np.mean(yarray)
                   stddev = np.std(yarray)
                   if abs(mean-val) > coef * stddev:
                       if verbose: print('Marking epoch {:} as outlier; component {:}'.format(t0, component))
                       # ts._data[count]['ignore'] = True
                       outlier_indexes.append(count)

        ## should we mark the data ?
        if verbose: print('Iteration {:}, Outliers marked: {:}/{:}'.format(iteration, len(outlier_indexes), starting_size))
        if len(outlier_indexes) > 0:
            more_outliers = True
            num_outliers = len(outlier_indexes)
            percent = (100e0 * num_outliers) / starting_size
            if percent < max_percent:
                for index in outlier_indexes:
                    ts._data[index]['ignore'] = True
                coef = 3e0
            else:
                print('Refusing to mark {:}% of observations. Criterion too stif, lowering 3-sigma rule'.format(percent))
                coef = coef + coef / 10e0
                iteration -= 1

    if delete:
        if verbose: print('Deleting outliers from time-series ...')
        return ts.drop_ignored()
    else:
        return ts
