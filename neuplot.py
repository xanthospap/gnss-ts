#! /bin/python

import sys, math
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.dates as mdates

month_day = [
    [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365],
    [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
]

def mjd2pydt(mjd_as_float):
    JAN11901          = 15385
    fmjd, mjd         = math.modf(mjd_as_float)
    mjd               = int(mjd)
    days_fr_jan1_1901 = int(mjd) - JAN11901
    num_four_yrs      = int(days_fr_jan1_1901)/1461
    years_so_far      = 1901 + 4*num_four_yrs
    days_left         = days_fr_jan1_1901 - 1461*num_four_yrs
    delta_yrs         = int(days_left)/365 - int(days_left)/1460
    year              = years_so_far + delta_yrs
    yday              = days_left - 365*delta_yrs + 1
    hour              = int(fmjd*24.0e0)
    minute            = int(fmjd*1440.0e0 - hour*60.0e0)
    second            = int(fmjd*86400.0e0 - hour*3600.0e0 - minute*60.0e0)
    leap              = int(year%4 == 0)
    guess             = int(yday*0.032)
    more              = int((yday - month_day[leap][guess+1]) > 0)
    month             = guess + more + 1
    mday              = yday - month_day[leap][guess+more]
    return datetime.datetime(year, month, mday, hour, minute, second)

def mjd2gpsw(mjd_as_float):
    JAN61980          = 44244
    SEC_PER_DAY       = 86400.0e0
    fmjd, mjd         = math.modf(mjd_as_float)
    mjd               = int(mjd)
    gps_week          = int((mjd - JAN61980)/7)
    sec_of_week       = ((mjd - JAN61980)-gps_week*7+fmjd)*SEC_PER_DAY
    return gps_week, sec_of_week

def sow2dow(sec_of_week):
    return sec_of_week/(7*SEC_PER_DAY)

def read_new_input(filename):
    ##  fill up the arrays ...
    epochs = []
    north  = []; sigma_north = []; north_flag = []
    east   = []; sigma_east  = []; east_flag = []
    up     = []; sigma_up    = []; up_flag = []
    ##  read the file, line by line ...
    ##  skip lines starting with '#'
    index = 0
    with open(filename, 'r') as fin:
        for line in fin:
            if line[0] != '#':
                l = line.split()
                index = 0
                epochs.append(float(l[index]))
                index += 1
                north.append(float(l[index]))
                index += 1
                sigma_north.append(float(l[index]))
                index += 1
                ## do we have a flag for north?
                if l[index].isalpha():
                    north_flag.append(l[index])
                    index += 1
                else:
                    north_flag.append('a')
                east.append(float(l[index]))
                index += 1
                sigma_east.append(float(l[index]))
                index += 1
                ## do we have a flag for east?
                if l[index].isalpha():
                    east_flag.append(l[index])
                    index += 1
                else:
                    east_flag.append('a')
                up.append(float(l[index]))
                index += 1
                sigma_up.append(float(l[index]))
                index += 1
                ## do we have a flag for east?
                if index < len(l) and l[index].isalpha():
                    up_flag.append(l[index])
                else:
                    up_flag.append('a')
    return epochs, north, sigma_north, north_flag, east, sigma_east, east_flag, up, sigma_up, up_flag

t, n, sn, fn, e, se, fe, u, su, fu = read_new_input(sys.argv[1])

plt.subplot(3, 1, 1)
plt.plot(t, n, 'yo-')
plt.title('A tale of 2 subplots')
plt.ylabel('North (m)')

plt.subplot(3, 1, 2)
plt.plot(t, e, 'r.-')
plt.xlabel('time (s)')
plt.ylabel('East (m)')

plt.subplot(3, 1, 3)
plt.plot(t, u, 'r.-')
plt.xlabel('time (s)')
plt.ylabel('Up (m)')

plt.show()
