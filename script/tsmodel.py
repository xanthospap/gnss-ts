#! /usr/bin/python

"""  
"""

from __future__ import print_function
import argparse
import sys
import operator
import datetime
import julian
import math

tmax = datetime.datetime.max
tmin = datetime.datetime.min

def p2mjd(t):
    return julian.to_jd(t) - 2400000.5

class Harmonic:

    def __init__(self, period_in_days, in_phase=0e0, out_of_phase=0e0, start_at=tmin, stop_at=tmax):
        self.angfreq = (math.pi*2e0)/period_in_days
        self.inphase   = in_phase
        self.outphse   = out_of_phase
        self.fro       = p2mjd(start_at)
        self.to        = p2mjd(stop_at)

    def value(self, t_mjd, tref_mjd):
        d = 0e0
        dt = t_mjd - tref_mjd
        if t_mjd >= self.fro and t_mjd <= self.to:
            d += self.inphase * math.cos(self.angfreq*dt)
            d += self.outphse * math.sin(self.angfreq*dt)
        return d

class Jump:
    def __init__(self, start_at=tmin, offset_in_meters=0e0):
        self.t_mjd = p2mjd(start_at)
        self.val   = offset_in_meters

    def value(self, t_mjd):
        return self.val if t_mjd > self.t_mjd else 0e0

class VelocityChange:

    def __init__(self, new_vel=0e0, start_at=tmin, stop_at=tmax):
        self.val = new_vel
        self.fro = p2mjd(start_at)
        self.to = p2mjd(stop_at)

    def value(self, t_mjd):
        d = 0e0
        if t_mjd >= self.fro and t <= self.to:
            dt = t_mjd - self.fro
            d += (dt/365.25e0)*self.val
        return d

class Earthquake:

    def __init__(self, start_at=tmin, model_nr=0, a1=0e0, t1=1e0, a2=0e0, t2=1e0):
        self.t = p2mjd(start_at)
        self.model = model_nr
        self.a1 = a1
        self.t1 = t1
        self.a2 = a2
        self.t2 = t2

    def value(self, t_mjd):
        if t_mjd < self.t: return 0e0
        if self.model == 0:
            return self.a1
        else:
            raise RuntimeError


class Model:

    def __init__(self, central_epoch=tmin):
        self.tref = p2mjd(central_epoch)
        self.harmonics = []
        self.jumps = []
        self.velocity_changes = []
        self.earthquakes = []

    def linear_terms(self, a, b):
        self.a = a
        self.b = b

    def get_event_list(self, use_pydatetime=False):
        def dtf(mjd, pydt=False):
            return mjd if not pydt else julian.from_jd(mjd, fmt='mjd')
        event_list = {'j':[], 'v':[], 'e':[]}
        for i in self.jumps:
            event_list['j'].append(dtf(i.t_mjd, use_pydatetime))
        for i in self.velocity_changes:
            event_list['v'].append(dtf(i.fro, use_pydatetime))
        for i in self.earthquakes:
            event_list['e'].append(dtf(i.t, use_pydatetime))
        return event_list

    def add_harmonic(self, period_in_days, in_phase=0e0, out_of_phase=0e0, start_at=tmin, stop_at=tmax):
        self.harmonics.append(Harmonic(period_in_days, in_phase, out_of_phase, start_at, stop_at))

    def add_jump(self, start_at=tmin, offset_in_meters=0e0):
        self.jumps.append(Jump(start_at, offset_in_meters))

    def add_velocity_change(self, new_vel=0e0, start_at=tmin, stop_at=tmax):
        self.velocity_changes.append(VelocityChange(new_vel, start_at, stop_at))

    def add_earthquake(self, start_at=tmin, model_nr=0, a1=0e0, t1=1e0, a2=0e0, t2=1e0):
        self.earthquakes.append(Earthquake(start_at, model_nr, a1, t1, a2, t2))

    def value(self, t_mjd):
        d = 0e0
        dt = (t_mjd - self.tref)/365.25e0
        d += self.a * dt
        d += self.b
        for i in self.harmonics:
            d += i.value(t_mjd, self.tref)
        for i in self.jumps:
            d += i.value(t_mjd)
        for i in self.velocity_changes:
            d += i.value(t_mjd)
        for i in self.earthquakes:
            d += i.value(t_mjd)
        return d

    def make_line(self, from_, to_, use_pydatetime=False):
        t0 = p2mjd(from_)
        t1 = p2mjd(to_)
        t  = []
        v  = []
        while t0 <= t1:
            t.append(t0)
            v.append(self.value(t0))
            t0 += 1
        if not use_pydatetime:
            return t, v
        else:
            d = [ julian.from_jd(mjd, fmt='mjd') for mjd in t ]
            return d, v

def model_from_ascii(filename):
    models = []
    with open(filename, 'r') as fin:
        line = fin.readline()
        while len(line) > 1:
            l = line.split()
            if l[0] != 'Central' or l[1] != 'Epoch':
                sys.stderr.write("[Error] Expected line: \"Central Epoch\"")
                raise RuntimeError 
            tref = datetime.datetime.strptime(l[3] + " " + l[4], "%Y-%m-%d %H:%M:%S")
            model = Model(tref)
            line = fin.readline()
            if line.strip() != 'Linear Coefficients:':
                sys.stderr.write("[Error] Expected line: \"Linear Coefficients:\"")
                raise RuntimeError
            model.b = float(fin.readline().strip())
            model.a = float(fin.readline().strip())
            line = fin.readline()
            while len(line) > 1:
                if line.strip() == 'Harmonic Coefficients:':
                    line = fin.readline()
                    while len(line.split()) > 2 and line.split()[1] == '(in-phase)':
                        l = line.split()
                        assert l[1] == '(in-phase)' and l[2] == 'From:'
                        if l[3] == 'start':
                            start = tmin
                            nxt = 4
                        else:
                            start = datetime.datetime.strptime(l[3] + " " + l[4], "%Y-%m-%d %H:%M:%S")
                            nxt = 5
                        assert l[nxt] == 'To:'
                        if l[nxt+1] == 'end':
                            stop = tmax
                            nxt  = nxt+2
                        else:
                            stop = datetime.datetime.strptime(l[nxt+1] + " " + l[nxt+2],"%Y-%m-%d %H:%M:%S")
                            nxt = nxt+3
                        assert l[nxt] == 'Period:'
                        period = float(l[nxt+1])
                        in_phase = float(l[0])
                        out_of_phase = float(fin.readline().split()[0])
                        model.add_harmonic(period, in_phase, out_of_phase, start, stop)
                        ##print('adding harmonic with period={:}'.format(period))
                        line = fin.readline()
                elif line.strip() == 'Jumps/Offsets:':
                    line = fin.readline()
                    while len(line.split()) > 2 and line.split()[1] == 'at':
                        l = fin.readline().split()
                        assert l[1] == 'at'
                        t = datetime.datetime.strptime(l[2] + " " + l[3], "%Y-%m-%d %H:%M:%S")
                        model.add_jump(t, float(l[0]))
                        ##print('adding jump with value={:}'.format(float(l[0])))
                        line = fin.readline()
                elif line.strip() == 'Velocity Changes:':
                    line = fin.readline()
                    while len(line.split()) > 2 and line.split()[1] == 'From':
                        l = line.split()
                        val = float(l[0])
                        assert l[1] == 'From:'
                        if l[2] == 'start':
                            start = tmin
                            nxt = 3
                        else:
                            start = datetime.datetime.strptime(l[2] + " " + l[3], "%Y-%m-%d %H:%M:%S")
                            nxt = 4
                        assert l[nxt] == 'To:'
                        if l[nxt+1] == 'end':
                            stop = tmax
                        else:
                            stop = datetime.datetime.strptime(l[nxt+1] + " " + l[nxt+2], "%Y-%m-%d %H:%M:%S")
                        model.add_velocity_change(val, start. stop)
                        line = fin.readline()
                elif line.strip() == 'PSD (Earthquake Post Deformation):':
                    line = fin.readline()
                    while len(line.split()) > 2 and line.split()[2] == 'earthquake:':
                        l = line.split()
                        assert l[0] + ' ' + l[1] + ' ' + l[2] == 'Time of earthquake:'
                        t = datetime.datetime.strptime(l[3] + " " + l[4], "%Y-%m-%d %H:%M:%S")
                        l = fin.readline().split()
                        assert l[0] == 'Offset:'
                        val = float(l[1])
                        model.add_earthquake(t, 0, val)
                        ##print('adding earthquake with val={:}'.format(val))
                        line = fin.readline()
                elif len(line.split()) > 2 and (line.split()[0] == 'Central' and line.split()[1] == 'Epoch'):
                    models.append(model)
                    break
                else:
                    raise RuntimeError
    models.append(model)
    return models

parser = argparse.ArgumentParser(
    description='Given a model description (in ascii), compute (x,y) value pairs.'
)
parser.add_argument('-f', '--file',
    action   = 'store',
    required = True,
    help     = 'The input file.',
    metavar  = 'INPUT_FILE',
    dest     = 'input_file'
)
parser.add_argument('-s', '--from',
    action   = 'store',
    required = True,
    help     = 'Starting date as \'yyyy-mm-ddThh:mm:ss. Month can be a 3-char string',
    metavar  = 'START_DATE',
    dest     = 'start_t'
)
parser.add_argument('-e', '--to',
    action   = 'store',
    required = True,
    help     = 'Ending date as \'yyyy-mm-ddThh:mm:ss. Month can be a 3-char string.',
    metavar  = 'END_DATE',
    dest     = 'end_t'
)
parser.add_argument('-m', '--mjd',
    action   = 'store_true',
    required = False,
    help     = 'Output dates written as MJD.',
    dest     = 'use_mjd',
)
parser.add_argument('-v', '--event-file',
    action   = 'store',
    required = False,
    help     = '(output) event file.',
    metavar  = 'EVENT_OUT_FILE',
    dest     = 'event_file'
)

##  Parse command line arguments
args = parser.parse_args()

status = 0
try:
    models = model_from_ascii(args.input_file)
    try:
        start_t = datetime.datetime.strptime(args.start_t, "%Y-%m-%dT%H:%M:%S")
    except:
        try:
            start_t = datetime.datetime.strptime(args.start_t, "%Y-%b-%dT%H:%M:%S")
        except:
            raise RuntimeError
    try:
        end_t = datetime.datetime.strptime(args.end_t, "%Y-%m-%dT%H:%M:%S")
    except:
        try:
            end_t = datetime.datetime.strptime(args.end_t, "%Y-%b-%dT%H:%M:%S")
        except:
            raise RuntimeError
    if len(models) != 3:
        sys.stderr.write("[Error] Failed to read 3 models from input file\n")
        raise RuntimeError
    yvals = []
    for mdl in models:
        t, y = mdl.make_line(start_t, end_t, not args.use_mjd)
        yvals.append(y)
    if not args.use_mjd:
        t = [ i.strftime('%Y-%m-%dT%H:%M:%S') for i in t ]
    for idx, eph in enumerate(t):
        print('{:} {:10.5f} {:10.5f} {:10.5f}'.format(eph, yvals[0][idx], yvals[1][idx], yvals[2][idx]))
    if args.event_file:
        event_dicts = []
        for i in models:
            event_dicts.append(i.get_event_list(not args.use_mjd))
        with open(args.event_file, 'w') as fout:
            for key in ['j', 'v', 'e']:
                for t in list(set(event_dicts[0][key] + event_dicts[1][key] + event_dicts[2][key])):
                    print('{:} {:}'.format(t, key), end="\n", file=fout)
except:
    status = 1

sys.exit(status)
