#! /usr/bin/python

import datetime
import time_series as ts
import model_fit as tm
from gatspy.periodic import LombScargleFast

def fractional_days(startt, endt):
    difference = endt - startt
    return difference.total_seconds() / datetime.timedelta(days=1).total_seconds()

##
##  **kwargs can be:
##    - when wanting to compute the periodogram only for one component, use
##      model = None or TsModel
##      coordinate_key = string (e.g. 'east', 'x', etc)
##    - when computing the periodogram for all three of a coordinate type (e.g.
## Topocentric, Cartesian, etc), then provide:
##    - model = [TsModel, TsModel, TsModel] or None, or []
##    - coordinate_type = CoordinateType (e.g. CoordinateType.Topocentric, etc)
##
#def lomb_scargle(timesr, models=None, ct=ts.CoordinateType.Topocentric):
def lomb_scargle(timesr, **kwargs):
    ## I guess we need to transform the datetimes to fractional days from 
    ## some start epoch. Let's suppose that the start epoch is the first epoch 
    ## in the time-series, minus a day
    tfrom, tto = timesr.time_span();
    t0 = tfrom - datetime.timedelta(days=1)

    ## We need to fix any jump included in the time-series and apply (aka
    ## remove) the trend. Hence, get the model's residuals ...
    if 'model' not in kwargs or (kwargs['model'] is None or kwargs['model'] == []):
        residuals = timesr
    else:
        if type(kwargs['model']) is list and (type(kwargs['model']) is list and len(kwargs['model'])==3):
            if 'coordinate_type' not in kwargs:
                errmsg = 'ERROR. Expecting key \'coordinate_type\' to compute periodogram (given multiple models)'
                raise RuntimeError(errmsg)
            residuals = timesr.residuals(model=kwargs['model'], coordinate_type=kwargs['coordinate_type'])

        else:
            if 'coordinate_key' not in kwargs:
                errmsg = 'ERROR. Expecting key \'coordinate_key\' to compute periodogram (given one model)'
                raise RuntimeError(errmsg)
            residuals = timesr.residuals(model=kwargs['model'], coordinate_key=kwargs['coordinate_key'])

    coordinate_components = [kwargs['coordinate_key']] if 'coordinate_key' in kwargs else ts.coordinate_type_keys(kwargs['coordinate_type'])

    ret = []
    ft = [ fractional_days(t0,t) for t in residuals.get('t') ]
    for c in coordinate_components:
        model = LombScargleFast().fit(ft, residuals.get(c), [1e-3]*residuals.size())
        ret.append(model)
    return ret
    # return model.periodogram_auto(nyquist_factor=100)
    # periods, power = model.periodogram_auto(nyquist_factor=100)
