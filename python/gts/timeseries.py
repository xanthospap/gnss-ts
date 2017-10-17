import numpy    as np
import datetime
#from datetime import timedelta
from functools import partial
from bisect import bisect_right
import copy

import geodesy  as geo
import tsflags
import gpstime

## for debuging
## from scipy import stats

def angular_frequency(period): return 2.0e0*np.pi / period;

class Model:

    def __init__(self):
        self.__t0__      = None
        self.__x0__      = 0e0; self.__vx__          = 0e0
        self.__periods__ = [];  self.__period_vals__ = [];
        self.__offsets__ = [];  self.__offset_vals__ = [];

    def set_central_epoch(self, t0):
        assert isinstance(t0, datetime.datetime)
        self.__t0__ = t0

    def add_periods(self, *args):
        """ Add one or more periods to the model.
            Parameters:
            -----------
            *args : any number of numeric (i.e. int or float) values. These
                    shall represent period in days.
        """
        for period in args:
            assert float(period)
            self.__periods__.append(period)

    def periods(self):
        return self.__periods__

    def add_offsets(self, *args):
        """ Add one or more offsets to the model.
            Parameters:
            -----------
            *args: any number of datetime.datetime instances to represent the
                   epochs the offsets occured.
        """
        for offset in args:
            assert isinstance(offset, datetime.datetime)
            self.__offsets__.append(offset)
    
    def offsets(self):
        return self.__offsets__

    def parameters(self):
        return 2 + 2*len(self.__periods__) + len(self.__offsets__)

    def make_model(self, start, stop, dt_in_days=1):
        """ Construct an 'artificial' time-series, using the model. The timeseries
            will span the interval [start,stop) with a step of dt_in_days days.
        """
        t = start
        dt_in_days = datetime.timedelta(days=dt_in_days)
        vals   = []
        epochs = []
        while t < stop:
            dt  = (t-self.__t0__).days / 365.25e0
            dtd = float( (t-self.__t0__).days )
            y = self.__x0__ + self.__vx__*dt
            for idx, fst in enumerate(self.__offsets__):
                if t >= fst:
                    y += self.__offset_vals__[idx]
            for idx, per in enumerate(self.__periods__):
                y += self.__period_vals__[idx*2] *np.cos(angular_frequency(per)*dtd)
                y += self.__period_vals__[idx*2+1]*np.sin(angular_frequency(per)*dtd)
            vals.append(y)
            epochs.append(t)
            t += dt_in_days
        return epochs, vals

    def assign(self, xestim):
        assert self.parameters() == xestim.shape[0]
        self.__x0__ = xestim[0]
        self.__vx__ = xestim[1]
        idx = 2
        if len(self.__offset_vals__): self.__offset_vals__ = []
        for fst in self.__offsets__:
            self.__offset_vals__.append(xestim[idx])
            idx += 1
        if len(self.__period_vals__): self.__period_vals__ = []
        for per in self.__periods__:
            self.__period_vals__.append(xestim[idx])
            idx += 1
            self.__period_vals__.append(xestim[idx])
            idx += 1
        return

    def assign_design_mat_row(self, t, y, w=1e0):
        t0     = self.__t0__
        dt     = (t-t0).days / 365.25e0
        dtd    = float( (t-t0).days )
        row    = [0]*self.parameters()
        row[0] = w*1.0e0
        row[1] = w*dt
        idx    = 2;
        for fst in self.offsets():
            row[idx] = 1e0*w if fst >= t else 0e0
            idx += 1
        for prd in self.periods():
            row[idx]   = w*np.cos(angular_frequency(prd)*dtd)
            row[idx+1] = w*np.sin(angular_frequency(prd)*dtd)
            idx += 2
        return row

    def print_model_info(self):
        print "Central Epoch: ", self.__t0__
        print "X0           : ", self.__x0__
        print "Velocity     : ", self.__vx__,"m/yr"
        for idx, ofs in enumerate(self.__offsets__):
            print "Offset at :", ofs," displacement: ", self.__offset_vals__[idx], "m"
        for idx, prd in enumerate(self.__periods__):
            print "Period of:", prd, "days", "out-of-phase: ", self.__period_vals__[idx*2], "in-phase: ",self.__period_vals__[idx*2+1], "m"
        return

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
                'sx_array', 'sy_array', 'sz_array',
                'time_stamp', 'comment']:
            try:
                kwargs[i]
            except:
                kwargs[i] = None
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
        ## the flag array, one flag per epoch
        self.flags       = [ tsflags.TsFlag() for i in range(0, len(self.epoch_array)) ]
        self.time_stamps = kwargs['time_stamp']
        self.comments    = kwargs['comment']
        assert self.check_sizes()

    def average_epoch(self):
        return self.epoch_array[0] + (self.epoch_array[self.size()-1] - self.epoch_array[0]) / 2

    def check_sizes(self):
        """ Check if all not None arrays of the instance are of equal size
        """
        sizes = []
        def foo(x):
            if type(x) is list:
                return len(x)
            if x is not None:
                return x.size
        sizes = [ foo(x) for x in [ self.x_array, self.y_array, self.z_array,
                                self.sx_array, self.sy_array, self.sz_array, 
                                self.flags, self.time_stamps, self.comments ]
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

    def time_span(self):
        """" Return the timespan of the time-series as a 'datetime.timedelta'
             instance, from min_epoch to max_epoch.
        """
        if self.epoch_array is None:
            return datetime.timedelta(0)
        return self.max_epoch() - self.min_epoch()

    def array_w_index(self, idx):
        if   idx == 0: return self.x_array, self.sx_array
        elif idx == 1: return self.y_array, self.sy_array
        elif idx == 2: return self.z_array, self.sz_array
        else         : raise RuntimeError

    def fit_model(self, cmp, mdl):
        """ cmp is [0,2]
        """
        model = copy.deepcopy(mdl)
        if not model.__t0__:
            model.set_central_epoch(self.average_epoch())
        x_array, sx_array = self.array_w_index(cmp)
        yls = np.zeros([self.size(), 1])
        Als = np.zeros([self.size(), model.parameters()])
        for i in range(0, self.size()):
            if not self.flags[i].check(tsflags.TsFlagOption.outlier):
                yls[i,0] = x_array[i]
                Als[i,:] = model.assign_design_mat_row(self.epoch_array[i], x_array[i])
        x, sos, rank, s = np.linalg.lstsq(Als, yls)
        rmse = np.asscalar(np.sqrt(sos / yls.size))
        residuals = yls - Als.dot(x)
        model.assign(x)
        model.print_model_info()
        return model, rmse, residuals

    def dummy_lin_fit(self, perform_outlier_filtering=False):
        new_outlier = True
        iteration   = 0
        while new_outlier and iteration <= 5:
            print 'Iteration {:2d}'.format(iteration)
            new_outlier = False
            for cmp in self.x_array, self.y_array, self.z_array:
                yls = []
                Als = []
                t0  = self.average_epoch()
                rel_pos = [] ## original index of each element in the system
                for i in range(0, self.size()):
                    if not self.flags[i].check(tsflags.TsFlagOption.outlier):
                        yls.append(cmp[i])
                        Als.append([1.0, (self.epoch_array[i] - t0).days/365.25])
                        rel_pos.append(i)
                y = np.array(yls) #.astype(float)
                A = np.array(Als) #.astype(float)
                x, sos, rank, s = np.linalg.lstsq(A, y)
                rmse = np.asscalar(np.sqrt(sos / y.size))
                print '\tV  = {:7.2f} +- {:7.2f} mm/y # of points {:4d} from {:} to {:}'.format(x[1]*1000.0, rmse*1000.0, y.size, min(self.epoch_array).date(), max(self.epoch_array).date())
                residuals = y - A.dot(x)
                for i,r in enumerate(residuals):
                    if abs(r) >= 3.0 * rmse:
                        self.flags[rel_pos[i]].set(tsflags.TsFlagOption.outlier)
                        print '\t\tOutlier at epoch {:} residual = {:+8.5f}mm 3*sigma = {:8.4f}'.format(self.epoch_array[rel_pos[i]], r, 3.0 * rmse)
                        new_outlier = True
            iteration += 1

    def split(self, epoch):
        """ Split the TimeSeries into two seperate TimeSeries at point epoch
            (this should be a datetime instance. The two TimeSeries are returned
            as a list (of two elements).
        """
        if self.epoch_array is None:
            raise RuntimeError
        if epoch < self.min_epoch():
            return None, self
        if epoch > self.max_epoch():
            return self, None
        'Find leftmost value greater than x'
        i = bisect_right(self.epoch_array, epoch)
        if i == len(self.epoch_array):
            raise RuntimeError
        ts_left = TimeSeries(name=self.station, type=self.crd_type, epoch_array=self.epoch_array[0:i],
            x_array=self.x_array[0:i], y_array=self.y_array[0:i],z_array=self.z_array[0:i],
            sx_array=self.sx_array[0:i], sy_array=self.sy_array[0:i],sz_array=self.sz_array[0:i],
            time_stamp=self.time_stamps[0:i], comment=self.comments[0:i])
        sz = len(self.epoch_array)
        ts_right= TimeSeries(name=self.station, type=self.crd_type, epoch_array=self.epoch_array[i:sz],
            x_array=self.x_array[i:sz], y_array=self.y_array[i:sz],z_array=self.z_array[i:sz],
            sx_array=self.sx_array[i:sz], sy_array=self.sy_array[i:sz],sz_array=self.sz_array[i:sz],
            time_stamp=self.time_stamps[i:sz], comment=self.comments[i:sz])
        return ts_left, ts_right
    
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

            Parameter:
            ----------
            crd_enum : string
                       Can be any of:
                            * 'Cartesian'
                            * 'Ellipsoidal'
                            * 'Topocentric'

            Returns:
            --------
            A (pointer to a) function. This function (lets call it fun), should
            be called like: x_new, y_new, z_new = fun(x_old, y_old, z_old)

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
            
            Parameter:
            ----------
            crd_enum : string
                       Can be any of:
                            * 'Cartesian'
                            * 'Ellipsoidal'
                            * 'Topocentric'
            Returns:
            --------
            A TimeSeries instance. This will be a copy of the calling instance
                                   with its [x,y,z] coordinates transformed to
                                   crd_enum type.
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

    def report_missing_days(self):
        dt  = self.epoch_array[0]
        idx = 0
        while dt <= self.epoch_array[self.size()-1]:
            if dt.date() != self.epoch_array[idx].date():
                week, sow = gpstime.pydt2gps(dt)
                print 'Missing date: {:} GPSt: {:04d}{:01d}'.format(dt, week, int(sow/86400.0))
            else:
                idx += 1
            dt += timedelta(days=1)

    def collect_outliers(self):
        ea = []
        ra = []
        fa = []
        for i in range(0, self.size()):
            if self.flags[i].check(tsflags.TsFlagOption.outlier):
                ea.append(self.epoch_array[i])
                ra.append([self.x_array[i], self.y_array[i], self.z_array[i], self.sx_array[i], self.sy_array[i], self.sz_array[i]])
                fa.append(tsflags.TsFlag(tsflags.TsFlagOption.outlier))
        return TimeSeries(name=self.station, type=self.crd_type,
                    epoch_array = ea,
                    x_array     = np.array(list(zip(*ra)[0] )),
                    y_array     = np.array(list(zip(*ra)[1] )),
                    z_array     = np.array(list(zip(*ra)[2] )),
                    sx_array    = np.array(list(zip(*ra)[3] )),
                    sy_array    = np.array(list(zip(*ra)[4] )),
                    sz_array    = np.array(list(zip(*ra)[5] )) )

    def toJson(self):
        jlst = []
        for i in range(0, self.size()):
            if not self.flags[i].check(tsflags.TsFlagOption.outlier):
                s = '{{ \"date\":\"{:}\", \"north\":{:10.5f}, \"ln\":{:10.5f}, \"un\":{:10.5f}, \"east\":{:10.5f}, \"le\":{:10.5f}, \"ue\":{:10.5f}, \"up\":{:10.5f}, \"lu\":{:10.5f}, \"uu\":{:10.5f}, \"flag\":\"{:}\" }}'.format(
                self.epoch_array[i].strftime('%Y-%m-%d'),
                self.x_array[i] * 1000.0,
                (abs(self.sx_array[i])*-1.0 + self.x_array[i]) * 1000.0,
                (abs(self.sx_array[i]) + self.x_array[i]) * 1000.0,
                self.y_array[i] * 1000.0,
                (abs(self.sy_array[i])*-1.0 + self.y_array[i]) * 1000.0,
                (abs(self.sy_array[i]) + self.y_array[i]) * 1000.0,
                self.z_array[i] * 1000.0,
                (abs(self.sz_array[i])*-1.0 + self.z_array[i]) * 1000.0,
                (abs(self.sz_array[i]) + self.z_array[i]) * 1000.0,
                self.flags[i])
                jlst.append(s)
        print ',\n'.join(jlst)

    def print_as_cts(self):
        for i in range(0, self.size()):
            print '{:s} {:+15.5f} {:9.5f} {:+15.5f} {:9.5f} {:+15.5f} {:9.5f} {:s} {:s}'.format(self.epoch_array[i].strftime('%Y-%m-%d %H:%M:%S'), self.x_array[i], self.sx_array[i], self.y_array[i], self.sy_array[i], self.z_array[i], self.sz_array[i], self.time_stamps[i].strftime('%Y-%m-%d %H:%M:%S'), self.comments[i])
