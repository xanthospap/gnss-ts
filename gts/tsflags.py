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
    
    def __init__(self, f=0):
        self.flag = f

    @property
    def flag(self): return self.__flag

    @flag.setter
    def flag(self, val):
        # val must be either a TsFlagOption instance or a valid integer
        if isinstance(val, TsFlagOption):
            self.__flag = val.value
            return
        if val >= pow(2, len(TsFlagOption)):
            raise RuntimeError('Invalid flag value!')
        self.__flag = val

    def set(self, *options):
        ''' Set (or add) one or more TsFlagOption(s)
        '''
        for i in options: self.__flag = self.__flag | i.value

    def clear(self, *options):
        ''' Clear one or more TsFlagOption(s). If options is empty, clear all
            flags.
        '''
        if len(options) == 0:
            self.__flag = 0
            return
        for i in options: self.__flag = self.__flag & (~i.value)

    def check(self, option):
        ''' Check if a specific TsFlagOption is set
        '''
        return self.__flag & option.value

    def __repr__(self):
        return 'flag:' + ' & '.join(filter(lambda a: a is not None,
            [ i.name if self.check(i) else None for i in TsFlagOption ]))

'''
    driver routine
'''

if __name__=="__main__":
    f1 = TsFlag(TsFlagOption.jump)
    f1.set(TsFlagOption.outlier, TsFlagOption.vel_chg)
    print f1

    f2 = TsFlag(f1.flag)
    print f2

    f2.clear()
    print f2
    
    # the following should not work
    try:
        f3 = TsFlag(f2.flag+1000)
        print f3
    except:
        print 'Exception caught!'
    
    fl = [TsFlag() for i in range(0,5)] * 5
    fl.append(f1)
    fl.append(f2)
    for i in fl: print i
    
    print type(fl), type(fl[0])
    (fl[0]).set(TsFlagOption.jump)
    for i in fl: print i
    
    f = fl[0]
    f.set(TsFlagOption.vel_chg)
    for i in fl: print i
