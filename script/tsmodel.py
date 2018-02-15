#! /usr/bin/python

"""  
"""

from __future__ import print_function
import argparse
import sys
import operator
import datetime
import julian

tmax = datetime.datetime.max
tmin = datetime.datetime.min

def p2mjd(t):
    return julian.to_jd(t) - 2400000.5

class Harmonic:

    def self.__init__(self, period_in_days, in_phase=0e0, out_of_phase=0e0, start_at=tmin, stop_at=tmax):
        self.__angfreq__ = (math.pi*2e0)/period_in_days
        self.__inphase   = in_phase
        self.__outphse   = out_of_phase
        self.__from      = p2mjd(start_at)
        self.__to        = p2mjd(stop_at)

    def value(self, t_mjd, tref_mjd):
        d = 0e0
        dt = t_mjd - tref_mjd
        if t_mjd >= self.__from and t =< self.__to:
            d += self.__inphase * math.cos(self.__angfreq*dt)
            d += self.__outphse * math.sin(self.__angfreq*dt)
        return d

class Jump:
    def self.__init__(self, start_at=tmin, offset_in_meters):
        self.__t_mjd = p2mjd(start_at)
        self.__val   = offset_in_meters

    def value(self, t_mjd):
        return self.__val if t_mjd > self.__t_mjd else 0e0

class VelocityChange:

    def self.__init__(self, new_vel=0e0, start_at=tmin, stop_at=tmax):
        self.__val = new_vel
        self.__from = p2mjd(start_at)
        self.__to = p2mjd(stop_at)

    def value(self, t_mjd):
        d = 0e0
        if t_mjd >= self.__from and t =< self.__to:
            dt = t_mjd - self.__from
            d += (dt/365.25e0)*self.__val
        return d

class Earthquake:

    def self.__init__(self, start_at=tmin, model_nr=0, a1=0e0, t1=1e0, a2=0e0, t2=1e0):
        self.__t = p2mjd(start_at)
        self.__model = model_nr
        self.__a1 = a1
        self.__t1 = t1
        self.__a2 = a2
        self.__t2 = t2

    def value(self, t_mjd):
       if t_mjd < self.__t: return 0e0
        if self.__model == 0:
            return self.__a1
        else:
            raise RuntimeError


class Model:
    pass

    def linear_terms(self, central_epoch, a, b):
        self.__tref = p2mjd(central_epoch)
        self.__a = a
        self.__b = b
        self.__harmonics = []
        self.__jumps = []
        self.__velocity_changes = []
        self.__earthquakes = []

    def add_harmonic(self, period_in_days, in_phase=0e0, out_of_phase=0e0, start_at=tmin, stop_at=tmax):
        self.__harmonics.append(Harmonic(period_in_days, in_phase, out_of_phase, start_at, stop_at))

    def add_jump(self, start_at=tmin, offset_in_meters):
        self.__jumps.append(Jump(start_at, offset_in_meters))

    def add_velocity_change(self, new_vel=0e0, start_at=tmin, stop_at=tmax):
        self.__velocity_changes.append(VelocityChange(new_vel, start_at, stop_at))

    def add_earthquake(self, start_at=tmin, model_nr=0, a1=0e0, t1=1e0, a2=0e0, t2=1e0):
        self.__earthquakes.append(Earthquake(start_at, model_nr, a1, t1, a2, t2))

    def value(self, t_mjd):
        d = 0e0
        dt = (t_mjd - self.__tref)/365.25e0
        d += self.__a * dt
        d += self.__b
        for i in self.__harmonics:
            d += i.value(t_mjd, self.__tref)
        for i in self.__jumps:
            d += i.value(t_mjd)
        for i in self.__velocity_changes:
            d += i.value(t_mjd)
        for i in self.__earthquakes:
            d += i.value(t_mjd)
        return d
