import numpy    as np
import datetime as dt
from functools import partial

import geodesy  as geo
import tsflags

##  dummy Enumeration type for holding coordinate types
def enum(**enums):
    return type('Enum', (), enums)
CoordinateType = enum(Cartesian=1, Ellipsoidal=2, Topocentric=3, Unknown=4)

class TimeSeries:
    ##  A simple TimeSeries class
    ##  All arrays are NumPy arrays except from epoch_array which a a list of
    ##  datetime instances.
    ##  Warning All arrays should be sorted
    def __init__(self, **kwargs):
        ##  all possible arguments set to None if not given
        for i in ['name',   'type',     'epoch_array',
                'x_array',  'y_array',  'z_array',
                'sx_array', 'sy_array', 'sz_array']:
            try:    kwargs[i]
            except: kwargs[i] = None
        self.station     = kwargs['name']
        self.crd_type    = kwargs['type']
        self.epoch_array = kwargs['epoch_array']
        ##  make sure all arrays are (or can be trasnformed to) floats !
        self.x_array     = kwargs['x_array']
        if self.x_array is not None: self.x_array = self.x_array.astype(float)
        self.y_array     = kwargs['y_array']
        if self.y_array is not None: self.y_array = self.y_array.astype(float)
        self.z_array     = kwargs['z_array']
        if self.z_array is not None: self.z_array = self.z_array.astype(float)
        self.sx_array    = kwargs['sx_array']
        if self.sx_array is not None: self.sx_array = self.sx_array.astype(float)
        self.sy_array    = kwargs['sy_array']
        if self.sy_array is not None: self.sy_array = self.sy_array.astype(float)
        self.sz_array    = kwargs['sz_array']
        if self.sz_array is not None: self.sz_array = self.sz_array.astype(float)
        assert self.check_sizes()
        ## the flag array, one flag per epoch
        self.flags       = len(self.epoch_array)*[ tsflags.TsFlag() ]

    def check_sizes(self):
        """ Check if all not None arrays of the instance are of equal size
        """
        sizes = []
        def foo(x):
            if x is not None: return x.size
        sizes = [ foo(x) for x in [self.x_array, self.y_array, self.z_array,
                                self.sx_array, self.sy_array, self.sz_array, 
                                self.flags]
                                if foo(x) is not None ]
        if self.epoch_array is not None: sizes.append(len(self.epoch_array))
        return sizes[1:] == sizes[:-1]

    def size(self):
        """ Return number of epochs, i.e. the size of the Timeseries
        """
        return len(self.epoch_array)

    def min_epoch(self):
        return self.epoch_array[0] if self.epoch_array is not None else None

    def max_epoch(self):
        return self.epoch_array[-1] if self.epoch_array is not None else None

    def array_w_index(self, idx):
        if   idx == 0: return self.x_array, self.sx_array
        elif idx == 1: return self.y_array, self.sy_array
        elif idx == 2: return self.z_array, self.sz_array
        else         : raise RuntimeError

    def split(self, epoch):
        """ Split the TimeSeries into two seperate TimeSeries at point epoch
            (this should be a datetime instance. The two TimeSeries are returned
            as a list (of two elements).
        """
        if self.epoch_array is None:
            raise RuntimeError
        if epoch < self.min_epoch():
            return [None, self]
        if epoch > self.max_epoch():
            return [self, None]
        # TODO
    
    def average(self, component=None, sigma=1.0):
        """ Compute and return the average values for all or for an individual
            component. If component is set to None, then averages values for all
            components will be computed and returned in a list; to compute the
            average for an individual component, provide its index [0-3) as value
            for the component argument.
            Foe every component that has a respective sigma array, the average will
            be the weighted average using a weight computed as:
            weight(i) = sigma / sigma_array(i).
        """
        averages = []
        if component is not None: comps = [ component ]
        else                    : comps = [ x for x in [0, 1, 2] ]
        for i in comps:
            x, s = self.array_w_index(i)
            if s is None:
                averages.append(np.average(x))
            else:
                weight_array = np.array([sigma/j for j in s])
                averages.append(np.average(x, weights=weight_array))
        return averages[0] if len(averages) == 1 else averages

    def select_tansform_function(self, crd_enum):
        """ Return a pointer to (or whatever the fuck this is named in python)
            the appropriate function to transform this TimeSeries to the given
            coordinate type (i.e. crd_enum is of type CoordinateType). If such a 
            function is not available just throw. If no transformation is needed,
            return None.
            The transformation function returned, should be of type:
            x_new, y_new, z_new = function(x_old, y_old, z_old)
        """
        if self.crd_type == crd_enum:
            return None
        assert ( self.crd_type is not CoordinateType.Topocentric and
            self.crd_type is not CoordinateType.Unknown )
        if ( self.crd_type is CoordinateType.Cartesian and
            crd_enum is CoordinateType.Ellipsoidal ):
            return geo.cartesian2ellipsoidal
        elif ( self.crd_type is CoordinateType.Cartesian and
            crd_enum is CoordinateType.Topocentric ):
            [ x_mean, y_mean, z_mean ] = self.average()
            ## bind the first three arguments to mean
            return partial(geo.cartesian2topocentric, x_mean, y_mean, z_mean)
        elif ( self.crd_type is CoordinateType.Ellipsoidal and
            crd_enum is CoordinateType.Cartesian ):
            return geo.ellipsoidal2cartesian
        elif ( self.crd_type is CoordinateType.Ellipsoidal and
            crd_enum is CoordinateType.Topocentric ):
            raise RuntimeError
        return None

    def transform(self, crd_enum):
        """ Transform a TimeSeries to the given CoordinateType. Note that:
            1. Only non-None arrays of the instance are transformed,
            2. Sigmas are NOT transformed
            TODO: Fix point 2. above.
            Warning. This function will **NOT** affect the self instance; a
            new, transformed TimeSeries will be returned.
        """
        ##  get the transformation function; this should be of type:
        ##  x_new, y_new, z_new = function(x_old, y_old, z_old)
        fun_ptr = self.select_tansform_function(crd_enum)
        if not fun_ptr:
            return self
        new_list = [ fun_ptr(x[0], x[1], x[2]) for x in 
            [ y for y in zip(self.x_array, self.y_array, self.z_array) ] ]
        return TimeSeries(name=self.station, type=crd_enum,
                    epoch_array = self.epoch_array,
                    x_array     = np.array([ x[0] for x in new_list ]),
                    y_array     = np.array([ x[1] for x in new_list ]),
                    z_array     = np.array([ x[2] for x in new_list ]),
                    sx_array    = self.sx_array,
                    sy_array    = self.sy_array,
                    sz_array    = self.sz_array)

'''
    def to_json(json_file):
        """ Output the timeseries instance to a file named json_file in
            json format.
        """
        with open(json_file, 'w') as jf:
'''
