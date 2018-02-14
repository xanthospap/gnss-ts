#! /usr/bin/python

"""  
"""

from __future__ import print_function
import argparse
import sys
import operator
import datetime

tmax = datetime.datetime.max
tmin = datetime.datetime.min

class Harmonic:

    def self.__init__(self, period, fromt=tmin, tot=tmax):
        self.__period__ = period
        self.__from__   = fromt
        self.__to__     = tot

    def value(self, t):
        

class Model:
    pass

    def linear_terms(a, b):
        self.__a = a
        self.__b = b

    def add_harmonic_term()
