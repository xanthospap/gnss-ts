import datetime
import os
import time_series as ts

##
## Parse DSO .cts files
##

def resolve_line(line):
    ## quick return if comment line ...
    if line.lstrip().startswith('#'): return None

    l = line.split()
    assert(len(l) >= 16)

    t = datetime.datetime.strptime(' '.join(l[0:2]), '%Y-%m-%d %H:%M:%S')
    x, sx, y, sy, z, sz, lat, slat, lon, slon, hgt, shgt = [
        float(c) for c in l[2:14]]
    
    ## time-stamp may have fractional seconds or not ...
    try :
        tstamp = datetime.datetime.strptime(
            ' '.join(l[14:16]),'%Y-%m-%d %H:%M:%S.%f')
    except:
        tstamp = datetime.datetime.strptime(
            ' '.join(l[14:16]),'%Y-%m-%d %H:%M:%S')

    comments = '' if len(l) == 16 else ' '.join(l[16:])

    return {'t': t, 'tstamp': tstamp, 'x': x, 'y': y, 'z': z, 'sx': sx, 'sy': sy, 'sz': sz, 'lat': lat, 'lon': lon, 'hgt': hgt, 'slat': slat, 'slon': slon, 'shgt': shgt, 'comments': comments}

## append line parsed from cts (for a single epoch) to a dictionary with 
## elements cts-lines (stored in dictionary) and keys their respective epochs
def dappend(dct, dentry):
    t = dentry['t']
    ## update
    if t in dct:
        dct[t] = dentry if dentry['tstamp'] > dct[t]['tstamp'] else dct[t]
    ## append
    else:
        dct[t] = dentry
    return dct

##
def comment_is_valid(allowed_comments, dct):
    return True if (dct['comments'] in allowed_comments or allowed_comments==['all']) else False

##
## kwargs allowed arguments:
##           from=t1 : datetime.datetime, default = datetime.datetime.min
##           to=t2   : datetime.datetime, default = datetime.datetime.max
##           unique  : bool, default = True
##  allowed_comments : list of strings, default=['all']
def parse_cts(filename, **kwargs):
    if not os.path.isfile(filename):
        ermsg = 'ERROR Failed to locate cts file {:}'.format(filename)
        raise RuntimeError(ermsg)

    tstart = kwargs['from'] if 'from' in kwargs else datetime.datetime.min
    tstop  = kwargs['to'] if 'to' in kwargs else datetime.datetime.max
    unique = kwargs['unique'] if 'unique' in kwargs else True
    allowed_comments = kwargs['allowed_comments'] if 'allowed_comments' in kwargs else ['all']

    ret = {} if unique else []

    with open(filename, 'r') as cts:
        for line in cts.readlines():
            dct = resolve_line(line)
            #print('resolved line to {:}'.format(dct))
            if dct:
                if dct['t'] >= tstart and dct['t'] < tstop and comment_is_valid(allowed_comments, dct):
                    dct['ignore'] = False
                    if not unique: 
                        ret.append(dct)
                    else:
                        ret = dappend(ret, dct)

    ## if needed drop keys and make list (from dictionary)
    if unique: ret = [ret[key] for key in ret]

    ## sort based on epoch
    return ts.TimeSeries('', sorted(ret, key=lambda d: d['t']))
