#! /usr/bin/python

import datetime
import argparse

JAN61980    = 44244
JAN11901    = 15385
SEC_PER_DAY = 86400.0e0
month_day   = [
  [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365],
  [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
]

def mjd2pydt(mjd):
    ''' Convert a date given in Modified Julian Date format, to
        a Python `datetime` object.

        :param mjd: The Modified Julian Date; can be `float`, `integer` or a
                    string convertible to float.
        :returns:   A Python `datetime` object, representing the input
                    (MJD) date.
    '''
    jd   = float(mjd)
    i, d = divmod(jd, 1)
    mjd  = int(i) ## just to be sure !
    fmjd = float(d)

    days_fr_jan1_1901 = mjd - JAN11901
    num_four_yrs      = int(days_fr_jan1_1901/1461)
    years_so_far      = 1901 + 4*num_four_yrs
    days_left         = days_fr_jan1_1901 - 1461*num_four_yrs
    delta_yrs         = int(days_left/365 - days_left/1460)

    year   = years_so_far + delta_yrs
    yday   = days_left - 365*delta_yrs + 1
    hour   = int(fmjd*24.0)
    minute = int(fmjd*1440.0 - hour*60.0)
    second = int(fmjd*86400.0 - hour*3600.0 - minute*60.0)
    leap   = int(year%4 == 0)
    guess  = int(yday*0.032)
    more   = int(( yday - month_day[leap][guess+1] ) > 0)
    month  = guess + more + 1
    mday   = yday - month_day[leap][guess+more]

    return datetime.datetime(year, month, mday, hour, minute, second)

def pydt2gmt(t):
    ''' Convert a datetime.datetime instance to a string of type:
        YYYY-MM-DDTHH:MM:SS
    '''
    return t.strftime("%Y-%m-%dT%H:%M:%S")

parser = argparse.ArgumentParser(
    description='Transform a datetime from MJD to YYYY-MM-DD HH:MM:SS format.'
)
parser.add_argument('mjd',
    type     = float,
    nargs    = '?',
    help     = 'Input datetime as MJD (float or int)'
)
parser.add_argument('-f', '--file',
    action   = 'store',
    required = False,
    help     = 'The input file (if any).',
    metavar  = 'INPUT_FILE',
    dest     = 'data_file'
)
parser.add_argument('-s', '--skip-rows',
    action   = 'store',
    type     = int,
    required = False,
    help     = 'Rows to skip off from the input file.',
    metavar  = 'SKIP_ROWS',
    dest     = 'skip_rows',
    default  = 0
)
parser.add_argument('-c', '--date-col',
    action   = 'store',
    type     = int,
    required = False,
    help     = 'Column of input file to transform.',
    metavar  = 'DATE_COL',
    dest     = 'date_col',
    default  = 0
)
parser.add_argument('-d', '--col-delimeter',
    action   = 'store',
    required = False,
    help     = 'Column seperator in input file.',
    metavar  = 'COL_SEPERATOR',
    dest     = 'col_seperator',
    default  = ' '
)

##  Parse command line arguments
args = parser.parse_args()

##  Quick exit
if not args.mjd and not args.data_file: exit(0)
if args.mjd and args.data_file:
    print '[ERROR] Cannot provide both an input file and a date!'
    exit(1)

##  One mjd given at input
if args.mjd:
    print '{:}'.format(pydt2gmt(mjd2pydt(args.mjd)))
    exit(0)

##  Input file given
if args.data_file:
    with open(args.data_file, 'r') as fin:
        if args.skip_rows > 0:
            for i in range(args.skip_rows):
                line = fin.readline()
                print line.rstrip("\n")
        for line in fin.readlines():
            l      = line.split(args.col_seperator)
            gmt_dt = pydt2gmt(mjd2pydt(l[args.date_col]))
            print line.replace(l[args.date_col], gmt_dt).rstrip("\n")
    exit(0)
