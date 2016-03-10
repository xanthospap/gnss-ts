#! /usr/bin/python

from enum import Enum, unique

@unique
class TsFlagOption(Enum):
    ''' An enumeration type containing (unique) flags with which to characterise
        an observation/measurement
    '''
    jump      = 1
    outlier   = 2
    vel_chg   = 4
    earthqk   = 8

class TsFlag(object):
    ''' A flag class to mark individual time-series epoch records
    '''
    
    def __init__(self, f):
        self.flag = f

    @property
    def flag(self): return self.__flag

    @flag.setter
    def flag(self, val):
        # val must be either a TsFlagOption instance or a valid integer
        if isinstance(val, TsFlagOption):
            self.__flag = val.value
            return
        if val >= len(TsFlagOption):
            raise RuntimeError('Invalid flag value!')
        self.__flag = val

    def set(self, *options):
        ''' Set (or add) one or more TsFlagOption(s)
        '''
        for i in options: self.__flag = self.__flag | i.value

    def clear(self, *options):
        ''' Clear one or more TsFlagOption(s)
        '''
        for i in options: self.__flag = self.__flag & (~i.value)

    def check(self, option):
        ''' Check if a specific TsFlagOption is set
        '''
        return self.__flag & option.value

    def __repr__(self):
        return 'flag:' + ' & '.join(filter(lambda a: a is not None, [ i.name if self.check(i) else None for i in TsFlagOption ]))


if __name__=="__main__":
    f1 = TsFlag(TsFlagOption.jump)
    f1.set(TsFlagOption.outlier, TsFlagOption.vel_chg)
    print f1

    f2 = TsFlag(f1.flag())
    print f2

    f3 = TsFlag(f2.flag()+1000)
    print f3
