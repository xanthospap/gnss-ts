#! /usr/bin/python

from __future__ import print_function
import datetime
import numpy as np
import math
import argparse

class Ellipsoid:
    ''' A class to represent a Reference Ellipsoid. It can be a "standard" ellipsoid,
          e.g. WGS84, GRS80, etc., or a user-defined one. This depends on the instance
          construction (see the constructor or whatever Pythons calls constructors by)
    '''
    def __init__(self, *args, **kwargs):
        ''' Constructing an ellipsoid can be done in one of two ways :
            * Construct a standard ellipsoid, e.g. ::
            refel = Ellipsoid('GRS80')
            use ``Ellipsoid('NAME')``, where ``'NAME'`` is one of:
            #. 'GRS80'
            #. 'WGS84'
            #. 'PZ90'
            * Construct a user-defined ellipsoid, e.g ::
            refel = Ellipsoid(a=6378137.0, f=.00335281, name='my ell'),
            where the parameter ``name`` is optional.
        '''
        if len(args) == 1: ## Construct from ellipsoid name
            if args[0] == 'GRS80':
                self.__a = 6378137.0e0
                self.__f = 1.0e00/298.257222101e0
                self.__name = 'GRS80'
            elif args[0] == 'WGS84':
                self.__a = 6378137.0e0
                self.__f = 1.0e00/298.257223563e0
                self.__name = 'WGS84'
            elif args[0] == 'PZ90':
                self.__a = 6378135.0e0
                self.__f = 1.0e00/298.257839303e0
                self.__name = 'PZ90'
            else:
                raise RuntimeError('Invalid ellipsoid name [%s]' %args[0])
        elif len(kwargs) >= 1: ## Construct a user-define ellipsoid
            self.__a = kwargs['a']
            self.__f = kwargs['f']
            if 'name' in kwargs:
                self.__name = kwargs['name']
            else:
                self.__name = 'user-defined'
        else:
            raise RuntimeError('Invalid ellipsoid initialization')

    def semiMajorAxis(self):
        ''' return the Semi-Major Axis (i.e. parameter ``a``)
        '''
        return self.__a

    def flattening(self):
        ''' return the Flattening (i.e. parameter ``f``)
        '''
        return self.__f

    def name(self):
        ''' return the ellipsoid's name
        '''
        return self.__name

    def eccentricitySquared(self):
        ''' return the Eccentricity squared (i.e. parameter ``e**2``)
        '''
        return ( 2.0e0 - self.__f ) * self.__f

    def semiMinorAxis(self):
        ''' return the Semi-Minor Axis (i.e. parameter ``b``)
        '''
        return self.__a - self.eccentricitySquared() * self.__a

def cartesian2ellipsoidal(x, y, z, ellipsoid=None):
    ''' Given a set of geocentric, cartesian coordinates and optionaly a reference
        ellispoid, transform the set to ellipsoidal coordinates, i.e. longtitude,
        latitude and (ellipsoidal) height.
    '''

    if ellipsoid == None: ellipsoid = Ellipsoid('GRS80')

    ## Functions of ellipsoid parameters.
    a     = ellipsoid.semiMajorAxis()
    f     = ellipsoid.flattening()
    aeps2 = a*a*1e-32
    e2    = (2.0e0-f)*f
    e4t   = e2*e2*1.5e0
    ep2   = 1.0e0-e2
    ep    = math.sqrt(ep2)
    aep   = a*ep

    ''' Compute Coefficients of (Modified) Quartic Equation
        Remark: Coefficients are rescaled by dividing by 'a'
    '''

    ## Compute distance from polar axis squared.
    p2 = x*x + y*y

    ## Compute longitude lon
    if p2 != .0e0:
        lon = math.atan2(y, x)
    else:
        lon = .0e0

    ## Ensure that Z-coordinate is unsigned.
    absz = abs(z)

    if p2 > aeps2: ## Continue unless at the poles
        ## Compute distance from polar axis.
        p = math.sqrt(p2)
        ## Normalize.
        s0 = absz/a
        pn = p/a
        zp = ep*s0
        ## Prepare Newton correction factors.
        c0  = ep*pn
        c02 = c0*c0
        c03 = c02*c0
        s02 = s0*s0
        s03 = s02*s0
        a02 = c02+s02
        a0  = math.sqrt(a02)
        a03 = a02*a0
        d0  = zp*a03 + e2*s03
        f0  = pn*a03 - e2*c03
        ## Prepare Halley correction factor.
        b0 = e4t*s02*c02*pn*(a0-ep)
        s1 = d0*f0 - b0*s0
        cp = ep*(f0*f0-b0*c0)
        ## Evaluate latitude and height.
        lat = math.atan(s1/cp)
        s12 = s1*s1
        cp2 = cp*cp
        hgt = (p*cp+absz*s1-a*math.sqrt(ep2*s12+cp2))/math.sqrt(s12+cp2);
    else: ## Special case: pole.
        lat = math.pi / 2e0;
        hgt = absz - aep;
    ## Restore sign of latitude.
    if z < 0.e0:
        lat = -lat

    return lon, lat, hgt

def resolve_cts_line(line):
    l = line.split()
    t = datetime.datetime.strptime(l[0]+' '+l[1], '%Y-%m-%d %H:%M:%S')
    x, sx, y, sy, z, sz = [ float(c) for c in l[2],l[3],l[4],l[5],l[6],l[7] ]
    tstamp = datetime.datetime.strptime(l[8]+' '+l[9], '%Y-%m-%d %H:%M:%S')
    comment = l[10].strip() if len(l)>9 else ""
    return {'t':t, 'x':x, 'sx':sx, 'y':y, 'sy':sy, 'z':z, 'sz':sz, 'stamp':tstamp, 'comment':comment}

def to_GITSA(infile, comment=None):
    data = []
    with open(infile, 'r') as fin:
        for line in fin.readlines():
            if line[0] != '#':
                if not comment:
                    data.append(resolve_cts_line(line))
                else:
                    infod = resolve_cts_line(line)
                    if infod['comment'] == comment: data.append(infod)
    mean_x = sum(k['x'] for k in data)/float(len(data))
    mean_y = sum(k['y'] for k in data)/float(len(data))
    mean_z = sum(k['z'] for k in data)/float(len(data))
    lon0, lat0, hgt0 = cartesian2ellipsoidal(mean_x, mean_y, mean_z)
    t0     = data[0]['t']
    sinl   = math.sin(lon0)
    cosl   = math.cos(lon0)
    sinp   = math.sin(lat0)
    cosp   = math.cos(lat0)
    R = np.matrix([[-sinl, cosl, 0e0],
        [-cosl*sinp, -sinl*sinp, cosp],
        [cosl*cosp, sinl*cosp, sinp]])
    for p in data:
        Dx = np.matrix([[p['x']-mean_x],[p['y']-mean_y],[p['z']-mean_z]])
        enu = R * Dx
        e = enu.item((0,0))
        n = enu.item((1,0))
        u = enu.item((2,0))
        decimal_years = float(p['t'].strftime('%j')) / 365.25e0 + float(p['t'].year)
        print('{:9.4f} {:4d} {:3s} {:13.4f} {:13.4f} {:13.4f} {:13.5e} {:13.5e} {:13.5e} {:7.4f} {:7.4f} {:7.4f} {:1d} {:1d} {:1d}'.format(decimal_years, p['t'].year, p['t'].strftime('%j'), n, e, u, p['sx'], p['sy'],p['sz'], 0, 0, 0, 0, 0, 0))

parser = argparse.ArgumentParser(
    description='''Rearange any file starting with a date field of type
    \'%Y-%m-%d %H:%M:%S\' in chronological order. The input file is left as is,
    while the output is directed to STDOUT''',
    formatter_class = argparse.RawTextHelpFormatter
)
parser.add_argument('-i', '--input-file',
    action   = 'store',
    required = True,
    help     = 'The input file (<station>.xyz).',
    metavar  = 'INPUT_FILE',
    dest     = 'data_file'
)
parser.add_argument('-c', '--comment',
    action   = 'store',
    required = False,
    help     = 'Only keep records with the given comment column.',
    metavar  = 'COMMENT',
    dest     = 'comment',
    default  = None
)

args = parser.parse_args()
to_GITSA(args.data_file, args.comment)
