import datetime
import numpy as np
import copy

def fractional_days(startt, endt):
    difference = endt - startt
    return difference.total_seconds() / datetime.timedelta(days=1).total_seconds()


class TsModel:

    def __init__(self, t0=None, include_linear_terms=True):
        self.t0 = t0
        ## [{'t':datetime.datetime, 'comment': 'receiver change'}, {...}]
        self.jumps = None 
        ## [{'dy': float , 'sdy': float (std-deviation)}, {...}]
        self.jump_estimates = None

        if include_linear_terms:
            self.intercept = None
            self.slope = None

    def append_jumps(self, jump_dict_list):
        if self.jumps is None:
            self.jumps = []
        for jump in jump_dict_list: self.jumps.append(jump)
        self.jumps = sorted(self.jumps, key=lambda d: d['t'])

    def value(self, t):
        linear = self.intercept + self.slope * fractional_days(self.t0, t)

        jump = 0e0;
        if self.jumps is not None:
            jump = sum([self.jump_estimates[i]['dy']
                        for i in range(len(self.jumps)) if self.jumps[i]['t'] > t])

        return linear + jump

    ## todo:: jumps are not added
    def fit(self, t, y, **kwargs):
        """
            kwargs can be:
            sy = [sigma_y1, sigma_y2, ...] aka the std.devaition values 
               corresponding to the input data
            t0 = datetime.datetime The central epoch. If not provided, then
               if the model includes a central epoch this is used. If this 
               is None, then the median epoch is used
        """

        # central epoch
        t0 = kwargs['t0'] if 't0' in kwargs else None
        if t0 is None:
            t0 = self.t0 if self.t0 is not None else t[0] + (t[-1] - t[0]) / 2

        numobs = len(t)
        assert(numobs == len(y))

        # least squares matrices
        dt = [fractional_days(t0, ti) for ti in t]
        A = np.vstack([dt, np.ones(numobs)]).T

        # least squares fit
        coefs, res, _, _ = np.linalg.lstsq(A, y, rcond=None)

        # make new (copy) model
        # mdl = TsModel()
        mdl = copy.deepcopy(self)
        mdl.intercept = coefs[0]
        mdl.slope = coefs[1]
        mdl.t0 = t0

        return mdl, res

    def synthetic(self, tstart, tend, tstep=datetime.timedelta(days=1), crd_name='x'):
        data = []
        t = tstart
        while t < tend:
            data.append({'t': t, 'x': self.intercept +
                         self.slope * fractional_days(self.t0, t)})
            t += tstep
        return TimeSeries('', data)
