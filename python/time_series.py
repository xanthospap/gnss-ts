from car2ell import cartesian2ellipsoidal
from topocentric_matrix import topocentric_matrix
from enum import Enum
import datetime
import numpy as np
import json

## dictionary.has_key(key) not available in Python 3
def p3_has_key(dct, key): return key in dct

class CoordinateType(Enum):
    Cartesian = 1
    Ellipsoidal = 2
    Topocentric = 3

def coordinate_type_keys(ct):
    if ct == CoordinateType.Cartesian:
        cnames = ['x', 'y', 'z']
    elif ct == CoordinateType.Ellipsoidal:
        cnames = ['lon', 'lat', 'hgt']
    else:
        cnames = ['east', 'north', 'up']
    return cnames

class TimeSeries:

    last = -999

    def __init__(self, site, entries=[]):
        self._site = site
        ## list of dictionaries ...
        self._data = entries
    
    def get(self, key):
        return [ entry[key] for entry in self._data ]

    def cut(self, tstart, tend=datetime.datetime.max):
        """ cut on time, aka tstart and tend are datetime.datetime 
        """
        entries = []
        for entry in self._data:
            if entry[t] >= tstart and entry[t] < tend:
                entries.append(entry)
            if entry[t] > tend:
                break
        return TimeSeries(self._site, entries)
    
    def icut(self, istart, iend=None):
        """ cut on index
        """
        if iend is None:
            iend = len(self._data)
        assert(istart >= 0 and iend <= len(self._data) and istart < iend)
        return TimeSeries(self._site, self._data[istart:iend])

    def includes_coordinate_type(self, ct):
        cnames = coordinate_type_keys(ct)
        #return all([self._data[0].has_key(k) for k in cnames])
        return all([ p3_has_key(self._data[0],k) for k in cnames ])

    def includes_key(self, key):
        return p3_has_key(self._data[0], key)
    
    def drop_coordinate_type(self,ct):
        if not self.includes_coordinate_type(ct):
            return TimeSeries(self._site, self._data)
        
        cnames = coordinate_type_keys(ct)
        snames = ['s'+ct for ct in cnames]
        unused = cnames + snames
        return TimeSeries(self._site, [{k: v for k, v in entry.items() if k not in unused} for entry in self._data])

    def sort(self, onkey='t'):
        ## onkey can be any of the keys included in the self._data list 
        ## elements
        return TimeSeries(self._site, sorted(self._data, key=lambda d: d[onkey]))

    def average_crd(self):
        ax = ay = az = 0e0
        n = 0
        for entry in self._data:
            if not entry['ignore']:
                ax += (entry['x'] - ax) / (n+1)
                ay += (entry['y'] - ay) / (n+1)
                az += (entry['z'] - az) / (n+1)
        return ax,ay,az

    def topocentric(self):
        if not self.includes_coordinate_type(CoordinateType.Cartesian):
            raise RuntimeError('Time-Series missing cartesian components\nCannot transform to topocentric\n')

        x, y, z = self.average_crd()
        lon, lat, hgt = cartesian2ellipsoidal(x,y,z)

        R = np.array(topocentric_matrix(lon, lat))
        
        newentries = []
        for i,entry in enumerate(self._data):
            dr = np.array([entry['x']-x, entry['y']-y, entry['z']-z])
            enu = R @ dr
            #newentries.append(dict(
            #    entry, {'east': enu[0], 'north': enu[1], 'up': enu[2]}))
            newentries.append({**entry, **{'east': enu[0], 'north': enu[1], 'up': enu[2]}})

        return TimeSeries(self._site, newentries)
    
    def size(self): return len(self._data)

    def time_span(self):
        assert len(self._data) > 1
        return self._data[0]['t'], self._data[-1]['t']

    def dump(self):
        return json.dumps(self._data, indent=4, sort_keys=True, default=str)
    
    def print(self, ct_2print=None):

        if ct_2print and not self.includes_coordinate_type(ct_2print):
            ermsg = 'Time-Series missing {:} components\nCannot print!\n'.format(
                ct_2print)
            raise RuntimeError(ermsg)

        unprint = []
        unprint += [coordinate_type_keys(ct) for ct in [CoordinateType.Cartesian,
                                                        CoordinateType.Ellipsoidal, CoordinateType.Topocentric] if self.includes_coordinate_type(ct) and ct != ct_2print]
        for entry in self._data:
            print('{:}'.format([entry[k] for k in entry if k not in unprint]))
    
    def last_prior_to(self, t, from_index=0, to_index=None):
        if to_index is None:
            to_index = len(self._data)
        
        if self._data[from_index]['t'] >= t:
            return TimeSeries.last, None

        last_i = 0
        last_e = self._data[from_index]
        for idx,entry in enumerate(self._data[from_index:to_index]):
            if entry['t'] >= t:
                break
            last_i = idx
            last_e = entry
        
        return last_i,last_e['t']
