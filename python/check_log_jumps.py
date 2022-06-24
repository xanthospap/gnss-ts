#! /usr/bin/python

import parsers as tparse
import time_series as tsp
import model_fit as mdl
import datetime
import os
import math
import pybern.products.formats.igs_log_file as blog

def parse_logInstance_date(tstr):
    tstr = tstr[0:-1] if tstr[-1] == 'Z' else tstr
    return datetime.datetime.strptime(tstr, '%Y-%m-%dT%H:%M')


def get_log_possible_jumps(logfn):
    """ Example return
        [
        {'t': datetime.datetime(2007, 5, 22, 0, 0), 'comment': 'Antenna/Receiver change (span00grc_20211117.log'}, 
        {'t': datetime.datetime(2008, 3, 27, 0, 0), 'comment': 'Receiver change (span00grc_20211117.log'}, 
        {'t': datetime.datetime(2012, 12, 12, 0, 0), 'comment': 'Receiver change (span00grc_20211117.log'}, 
        {'t': datetime.datetime(2019, 6, 24, 0, 0), 'comment': 'Receiver change (span00grc_20211117.log'}
        ]
    """
    log = blog.IgsLogFile(logfn)
    posjumps = []

    # if a change is at the day of installement, it is skipped (not a change)
    name, domes, installed_at = log.site_name()
    # date to datetime
    installed_at = datetime.datetime.combine(
        installed_at, datetime.datetime.min.time())

    # parse receiver changes
    dct = log.parse_block(3)
    """
        {'3.1': {'Receiver Type': 'LEICA GRX1200PRO', ...}, 
         '3.2': {'Receiver Type': 'LEICA GRX1200PRO', ...},
        }
    """
    for key, value in dct.items():
        t = parse_logInstance_date(value['Date Installed'])
        if not installed_at.date() == t.date():
            posjumps.append(
                {'t': t, 'comment': 'Receiver change ({:})'.format(os.path.basename(logfn))})

    # parse antenna changes; if an antenna change is also a receiver change,
    # just change the corresponding entr. Else, append the antenna change
    dct = log.parse_block(4)
    for key, value in dct.items():
        t = parse_logInstance_date(value['Date Installed'])
        if not installed_at.date() == t.date():
            isReceiverChange = False
            for c, entry in enumerate(posjumps):
                if entry['t'] == t:
                    posjumps[c]['comment'] = 'Antenna/Receiver change ({:})'.format(
                        os.path.basename(logfn))
                    isReceiverChange = True
            if not isReceiverChange:
                posjumps.append(
                    {'t': t, 'comment': 'Antenna change ({:})'.format(os.path.basename(logfn))})

    return posjumps


def check_jump(ts, jumpdct, component, max_data_points=30, min_data_points=15):
    detected_jumps = []  # to be populated and returned ...

    if not ts.includes_key(component):
        errmsg = 'Error. Time-Series has no component named \'{:}\''.format(
            component)
        raise RuntimeError(errmsg)

    t = jumpdct['t']  # datetime of event
    print(
        '-> Checking possible jump at {:} reason: {:}'.format(t, jumpdct['comment']))
    tstart, tend = ts.time_span()
    if t <= tstart or t >= tend:
        errmsg = 'Error. Checking for jump for date outside Time-Series limits! t={:}'.format(
            t)
        raise RuntimeError(errmsg)

    # find the jump (index) in the time-series
    tidx, tts = ts.last_prior_to(t)

    # ------------------------------------
    # check the interval prior to the jump
    # ------------------------------------

    # indexes ...
    idx_end = tidx
    if idx_end < min_data_points:
        errmsg = 'Error. Too few points to check possible jump at t={:}'.format(
            t)
        raise RuntimeError(errmsg)
    idx_start = idx_end - min_data_points
    while (idx_end - idx_start < max_data_points) and (idx_start > 0):
        idx_start -= 1
    # cut time-series
    #print('\tCutting time-series from {:} to {:}'.format(idx_start, idx_end))
    ts_cut = ts.icut(idx_start, idx_end+1)
    #print('\tCut time-series from {:} to {:}'.format(ts_cut.time_span()[0], ts_cut.time_span()[1]))
    # fit linear model and get slope
    lmd = mdl.TsModel(t)
    x, res = lmd.fit(ts_cut.get('t'), ts_cut.get(component))
    estimate_prior = x.value(t)
    rms2_prior = res[0] / (idx_end-idx_start)
    #print('\tslope before is: {:}, res is {:}, estimate is {:}'.format(x.slope,res,estimate_prior))

    # ------------------------------------
    # check the interval after the jump
    # ------------------------------------

    # indexes ...
    idx_start = tidx + 1
    if idx_start > ts.size() - min_data_points:
        errmsg = 'Error. Too few points to check possible jump at t={:}'.format(
            t)
        raise RuntimeError(errmsg)
    idx_end = idx_start + min_data_points
    while (idx_end - idx_start < max_data_points) and (idx_end < ts.size() - 1):
        idx_end += 1
    # cut time-series
    #print('\tCutting time-series from {:} to {:}'.format(idx_start, idx_end))
    ts_cut = ts.icut(idx_start, idx_end)
    #print('\tCut time-series from {:} to {:}'.format(ts_cut.time_span()[0], ts_cut.time_span()[1]))
    # fit linear model and get slope
    lmd = mdl.TsModel(t)
    x, res = lmd.fit(ts_cut.get('t'), ts_cut.get(component))
    estimate_after = x.value(t)
    rms2_after = res[0] / (idx_end-idx_start)
    #print('\tslope after is: {:}, res is {:} estimate is {:}'.format(x.slope,res, estimate_after))

    # check if the posibble jump, is indeed a jump ...
    slope_std = math.sqrt(rms2_prior + rms2_after)
    print(
        '\tDy = {:.6f}, Std = +-{:.7f}'.format(abs(estimate_after - estimate_prior), slope_std))
    if abs(estimate_after - estimate_prior) > 3*slope_std:
        print('WARNING! Slope detected!')
        detected_jumps.append(jumpdct)

    return detected_jumps

# Test Program
# ---------------------------
#import sys
#ts = tparse.parse_cts(sys.argv[2])
# ts = ts.topocentric().drop_coordinate_type(
#    tsp.CoordinateType.Cartesian).drop_coordinate_type(tsp.CoordinateType.Ellipsoidal)
#d = get_log_possible_jumps(sys.argv[1])
#for j in d: check_jump(ts, j)
