#! /bin/python

import sys, math, datetime, argparse
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.dates as mdates

## defines for datetime algorithms
month_day = [
    [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365],
    [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
]
JAN61980          = 44244
JAN11901          = 15385
SEC_PER_DAY       = 86400.0e0

## Mapping characteristics for each possible flag
flag_dict = {
    'j': {'what': 'jump', 'color': 'k'},
    'e': {'what': 'earthquake', 'color': 'g'},
    'o': {'what': 'outlier', 'color': 'm'},
    'v': {'what': 'velocity_change', 'color': 'b'}
}

## Default x-label
x_label_str = 'Modified Julian Date'

def mjd2pydate(mjd_as_float):
    """ Modified Julian Date (as float) to a Python datetime instance """
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
    """ Modified Julian Date (as float) to gps week and seconds of week """
    fmjd, mjd         = math.modf(mjd_as_float)
    mjd               = int(mjd)
    gps_week          = int((mjd - JAN61980)/7)
    sec_of_week       = ((mjd - JAN61980)-gps_week*7+fmjd)*SEC_PER_DAY
    return gps_week, sec_of_week

def pydate2mjd(epoch):
    """ Python datetime instance to Modified Julian Date (integral and fractional part)"""
    year  = int(epoch.year)
    month = int(epoch.month)
    dom   = int(epoch.day)
    hr    = int(epoch.hour)
    mn    = int(epoch.minute)
    ss    = int(epoch.second)
    leap  = int(year%4 == 0)
    yday  = int(month_day[leap][month-1] + dom)
    mjd   = int((year-1901)/4)*1461 + ((year-1901)%4)*365+yday-1+JAN11901
    fmjd  = ((ss/60.0 + mn)/60.0 + hr)/24.0
    return mjd, fmjd

def read_evn_input(filename):
    """ Read in an event file (.evn) and return a list (of lists) as:
        [ [epoch_as_mjd, event_flag_as_char],
          [epoch_as_mjd, event_flag_as_char], ...]
    """
    ##  fill up the arrays ...
    event_list = []
    with open(filename, 'r') as fin:
        for line in fin:
            if len(line) > 2 and line[0] != '#' and line[0:2] != 'YY':
                epoch = datetime.datetime.strptime(line[0:21].rstrip(), '%Y-%m-%d %H:%M:%S')
                flag  = line[21:33].strip()
                event_list.append( [pydate2mjd(epoch)[0]+pydate2mjd(epoch)[1], flag] )
    print '## Number of events in ts: {0:3d}'.format(len(event_list))
    return event_list

def get_outliers(epochs, x, sigma_x, flag_x):
    new_epochs = []
    new_x = []
    new_sigma_x = []
    counter = 0
    for i in flag_x:
        if "o" in i:
            new_epochs.append(epochs[counter])
            new_x.append(x[counter])
            new_sigma_x.append(sigma_x[counter])
        counter += 1
    return new_epochs, new_x, new_sigma_x

def read_new_input(filename):
    """ Read in a cts file; this needs more checking ...
        Returns a list (of lists) as:
        [ [epoch_as_mjd, north, sigma_north, north_flag, east, sigma_east, east_flag, up, sigma_up, up_flag],
          [...], ... ]
    """
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


## Set the cmd parser
## --------------------------------------------------------------------------- 
parser = argparse.ArgumentParser(
    description='Simple Python time-series plots.')
parser.add_argument('-i', '--input',
    action   = 'store',
    required = True,
    help     = 'A time-series file.',
    metavar  = 'TIME_SERIES',
    dest     = 'neu_file'
)
parser.add_argument('-e', '--event-file',
    action   = 'store',
    required = False,
    help     = 'An event list file.',
    metavar  = 'EVENT_LIST',
    dest     = 'event_file',
    default  = None
)
parser.add_argument('-t', '--time-format',
    action   = 'store',
    required = False,
    help     = 'Choose the time format of the plot.',
    metavar  = 'DATETIME_FORMAT',
    dest     = 'time_format',
    choices  = ['mjd', 'gps', 'ymd'],
    default  = 'mjd'
)

## Parse command line arguments
args = parser.parse_args()

## Read the timeseries .neu file
t, n, sn, fn, e, se, fe, u, su, fu = read_new_input(args.neu_file)

## Read the .evn file (if it exists)
if args.event_file is not None:
    events = read_evn_input(args.event_file)
else:
    events = None

## Change datetime format to gps weeks
if args.time_format == 'gps':
    sec_in_week = 7 * SEC_PER_DAY
    t = [ mjd2gpsw(i)[0] + mjd2gpsw(i)[1]/sec_in_week for i in t ]
    if events is not None:
        events = [ [mjd2gpsw(i[0])[0] + mjd2gpsw(i[0])[1]/sec_in_week, i[1]] for i in events ]
    x_label_str = 'Gps Week'

## Change datetime format to y/m/d
if args.time_format == 'ymd':
    t = [ mjd2pydate(i) for i in t ]
    if events is not None:
        events = [ [mjd2pydate(i[0]), i[1]] for i in events ]
    x_label_str = 'Date (YYYY/MM/DD)'

## Plot (or subplot)
plt.subplot(3, 1, 1)
plt.plot(t, n, 'yo-')
oe, ox, osx = get_outliers(t, n, sn, fn)
plt.plot(oe, ox, 'ro')
if events:
    for event in events:
        plt.axvline(event[0], linewidth=.5, color=flag_dict[event[1]]['color'])
plt.title('A tale of 2 subplots')
plt.ylabel('North (m)')

plt.subplot(3, 1, 2)
plt.plot(t, e, 'r.-')
oe, ox, osx = get_outliers(t, e, se, fe)
plt.plot(oe, ox, 'bo')
if events:
    for event in events:
        plt.axvline(event[0], linewidth=.5, color=flag_dict[event[1]]['color'])
plt.ylabel('East (m)')

plt.subplot(3, 1, 3)
plt.plot(t, u, 'r.-')
oe, ox, osx = get_outliers(t, u, su, fu)
plt.plot(oe, ox, 'bo')
if events:
    for event in events:
        plt.axvline(event[0], linewidth=.5, color=flag_dict[event[1]]['color'])
plt.xlabel(x_label_str)
plt.ylabel('Up (m)')

## all done! show
plt.show()
