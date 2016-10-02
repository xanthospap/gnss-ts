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
    'j': {'what': 'jump', 'color': 'g'},
    'e': {'what': 'earthquake', 'color': 'b'},
    'o': {'what': 'outlier', 'color': 'y'},
    'v': {'what': 'velocity_change', 'color': 'k'}
}

## Default x-label
x_label_str = 'Modified Julian Date'

def parse_axis_limits(limit_string):
    """ The axis limits can be any of:
        xmin/xmax
        xmin/xmax,ymin/ymax
        xmin/xmax,ymin/ymax,zmin/zmax
        if xmin is negative, it can be given as "n3.23" instead of "-3.23"
    """
    if limit_string[0] == "n" :
        limit_str = '-' + limit_string[1:]
    else:
        limit_str = limit_string
    l = [ i.split("/") for i in limit_str.split(",") ]
    xmin, xmax = [float(i) for i in l[0]]
    ymin, ymax = [float(i) for i in l[1]] if len(l) > 1 else [None]*2
    zmin, zmax = [float(i) for i in l[2]] if len(l) > 2 else [None]*2
    return xmin, xmax, ymin, ymax, zmin, zmax

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
    outlier_char = ''
    for i in flag_dict:
        if flag_dict[i]['what'] == 'outlier':
            outlier_char = i
            break
    assert outlier_char != ''
    new_epochs = []
    new_x = []
    new_sigma_x = []
    for index, flag in enumerate(flag_x):
        if outlier_char in flag:
            new_epochs.append(epochs[index])
            new_x.append(x[index])
            new_sigma_x.append(sigma_x[index])
    return new_epochs, new_x, new_sigma_x

def remove_outliers(t, x, sx, fx):
    outlier_char = ''
    for i in flag_dict:
        if flag_dict[i]['what'] == 'outlier':
            outlier_char = i
            break
    assert outlier_char != ''
    ll = zip(*[ i for i in zip(t, x, sx, fx) if outlier_char not in i[3]])
    return list(ll[0]), list(ll[1]), list(ll[2]), list(ll[3])

def read_model_line_file(filename):
    ## fill up the arrays ...
    epochs = []
    x = []; y = []; z = [];
    with open(filename, 'r') as fin:
        for line in fin:
            if len(line) >= 2 and line[0] != '#':
                l = line.split()
                epochs.append(float(l[0]))
                x.append(float(l[1]))
                y.append(float(l[2]))
                z.append(float(l[3]))
    return epochs, x, y, z

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
parser.add_argument('-m', '--model-file',
    action   = 'store',
    required = False,
    help     = 'A file containing the model line.',
    metavar  = 'MODEL_LINE',
    dest     = 'model_file',
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
parser.add_argument('-s', '--smooth',
    action   = 'store_true',
    required = False,
    help     = 'Smooth the time-series, i.e. do not show outliers.',
    dest     = 'hide_outliers'
)
parser.add_argument('-l', '--axis-limits',
    action   = 'store',
    required = False,
    help     = 'Choose the axis limits. The format should be: \"xmin/xmax[,ymin/ymax[,zmin/zmax]]\" '
    'If the first argument (aka xmin) is negative, it should be given with a leading \"n\" '
    'and NOT the negative sign (e.g. \"n3.00/3.00,-5.0/5.0\" NOT \"-3.00/3.00,-5.0/5.0\". '
    'This is an argparse bug.',
    metavar  = 'AXIS_LIMITS',
    dest     = 'axis_limits',
    default  = None
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

##  Read the .mod file (if it exists)
if args.model_file is not None:
    modt, modx, mody, modz = read_model_line_file(args.model_file)
else:
    modt = None

##  Parse the axis limits (if given)
if args.axis_limits is not None:
    xmin, xmax, ymin, ymax, zmin, zmax = parse_axis_limits(args.axis_limits)

## Change datetime format to gps weeks
if args.time_format == 'gps':
    sec_in_week = 7 * SEC_PER_DAY
    t = [ mjd2gpsw(i)[0] + mjd2gpsw(i)[1]/sec_in_week for i in t ]
    if events is not None:
        events = [ [mjd2gpsw(i[0])[0] + mjd2gpsw(i[0])[1]/sec_in_week, i[1]] for i in events ]
    if modt is not None:
        modt = [ [mjd2gpsw(i)[0] + mjd2gpsw(i)[1]/sec_in_week, i[1]] for i in modt ]
    x_label_str = 'Gps Week'

## Change datetime format to y/m/d
if args.time_format == 'ymd':
    t = [ mjd2pydate(i) for i in t ]
    if events is not None:
        events = [ [mjd2pydate(i[0]), i[1]] for i in events ]
    if modt is not None:
        modt = [ mjd2pydate(i) for i in modt ]
    x_label_str = 'Date (YYYY/MM/DD)'

## Get the min and max epochs
min_epoch = t[0]
max_epoch = t[len(t)-1]

## Plot (or subplot)
axis_x = plt.subplot(3, 1, 1)
if args.hide_outliers:
    ho_t, ho_n, ho_sn, ho_fn = remove_outliers(t, n, sn, fn)
    plt.plot(ho_t, ho_n, 'c.')
else:
    plt.plot(t, n, 'c.')
    ## Plot outliers
    oe, ox, osx = get_outliers(t, n, sn, fn)
    plt.plot(oe, ox, 'yo')
## Plot events
if events:
    for event in events:
        plt.axvline(event[0], linewidth=.5, color=flag_dict[event[1]]['color'])
## Plot model line
if modt:
    plt.plot(modt, modx, '-')
## Set axis limits
axis_x.set_xlim(min_epoch, max_epoch)
if args.axis_limits is not None:
    axis_x.set_ylim([xmin, xmax])

plt.title('A tale of 2 subplots')
plt.ylabel('North (m)')

axis_y = plt.subplot(3, 1, 2)
if args.hide_outliers:
    ho_t, ho_e, ho_se, ho_fe = remove_outliers(t, e, se, fe)
    plt.plot(ho_t, ho_e, 'c.')
else:
    plt.plot(t, e, 'r.')
    oe, ox, osx = get_outliers(t, e, se, fe)
    plt.plot(oe, ox, 'yo')
if events:
    for event in events:
        plt.axvline(event[0], linewidth=.5, color=flag_dict[event[1]]['color'])
if modt:
    plt.plot(modt, mody, '-')
## Set axis limits
axis_y.set_xlim(min_epoch, max_epoch)
if args.axis_limits is not None and ymin is not None:
    axis_y.set_ylim([ymin, ymax])
plt.ylabel('East (m)')

axis_z = plt.subplot(3, 1, 3)
if args.hide_outliers:
    ho_t, ho_u, ho_su, ho_fu = remove_outliers(t, u, su, fu)
    plt.plot(ho_t, ho_u, 'c.')
else:
    plt.plot(t, u, 'g.')
    oe, ox, osx = get_outliers(t, u, su, fu)
    plt.plot(oe, ox, 'yo')
if events:
    for event in events:
        plt.axvline(event[0], linewidth=.5, color=flag_dict[event[1]]['color'])
if modt:
    plt.plot(modt, modz, '-')
if args.axis_limits is not None and zmin is not None:
    axis_y.set_ylim([zmin, zmax])
axis_z.set_xlim(min_epoch, max_epoch)
plt.xlabel(x_label_str)
plt.ylabel('Up (m)')

## all done! show
plt.show()
