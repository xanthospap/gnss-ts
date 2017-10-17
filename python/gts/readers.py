import numpy      as np
import datetime   as dt
import timeseries as ts
import bisect

JAN11901 = 15385
month_day = [
    [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365],
    [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
]

def mjd2pydt(mjd,fmjd=.0):
    jd   = mjd + fmjd
    i, d = divmod(jd, 1)
    mjd  = int(i) ## just to be sure !
    fmjd = float(d)
    days_fr_jan1_1901 = mjd - JAN11901
    num_four_yrs      = days_fr_jan1_1901/1461
    years_so_far      = 1901 + 4*num_four_yrs
    days_left         = days_fr_jan1_1901 - 1461*num_four_yrs
    delta_yrs         = days_left/365 - days_left/1460
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

    return dt.datetime(year, month, mday, hour, minute, second)

def read_plt(plt_file):
    with open(plt_file, 'r') as fin:
        epochs = []
        x_ar   = []
        y_ar   = []
        z_ar   = []

        line = fin.readline()
        while line:
            l = line.split()
            namex = l[0]
            record_nrx = int(l[1])
            residualx = float(l[2])
            mjdx = l[3]
            l = fin.readline().split()
            namey = l[0]
            record_nry = int(l[1])
            residualy = float(l[2])
            mjdy = l[3]
            l = fin.readline().split()
            namez = l[0]
            record_nrz = int(l[1])
            residualz = float(l[2])
            mjdz = l[3]
            if namex != namey or namex != namez:
                raise RuntimeError("Invalid PLT format (names)")
            if record_nrx != record_nry or record_nrx != record_nrz:
                raise RuntimeError("Invalid PLT format (records)")
            if mjdx != mjdy or mjdx != mjdz:
                raise RuntimeError("Invalid PLT format (epochs)")
            x_ar.append( residualx )
            y_ar.append( residualy )
            z_ar.append( residualz )
            epochs.append( mjd2pydt(float(mjdx)) )
            line = fin.readline()
        return ts.TimeSeries(name=namex, type=ts.CoordinateType.Topocentric,
            epoch_array = epochs,
            x_array = np.array( x_ar ),
            y_array = np.array( y_ar ),
            z_array = np.array( z_ar ))

def sorted_appender(epoch_array, rest_array, time_stamps, new_epoch, new_record, t_stamp):
    if len(epoch_array) == 0:
            epoch_array.append(new_epoch)
            rest_array.append(new_record)
            time_stamps.append(t_stamp)
            return
    while True:
        try:
            idx = epoch_array.index(new_epoch)
            if time_stamps[idx] > t_stamp:
                print '## [DEBUG] Record for date {:} already included with more recent time-stamp'.format(epoch_array[idx])
                return
            elif time_stamps[idx] == t_stamp:
                return
            else:
                epoch_array[idx] = new_epoch
                rest_array[idx]  = new_record
                time_stamps[idx] = t_stamp
        except ValueError:
            idx = bisect.bisect_right(epoch_array, new_epoch)
            epoch_array.insert(idx, new_epoch)
            rest_array.insert(idx, new_record)
            time_stamps.insert(idx, t_stamp)
            return

def read_cts(cts_file, comment=None):
    ''' Read a raw .cts file.
        Lines starting with '#' are considered comment lines and are ignored.
        If the input parameter 'comment' is used, then only lines ending with
        this string are considered (i.e. lines for which line(196:) == comment).
        If 'comment' is not used, then all lines are considered.
        Valid lines should have at least 16 fields.

        Prameters:
        ----------
        cts_file : string
                   The filename of the .cts file to be read
        comment  : string

    '''
    with open(cts_file, 'r') as fin:
        epochs_array = []
        rest_array   = []
        time_stamps  = []
        for line in fin.readlines():
            if not line[0] == '#':
                line_cmnt = line[196:].rstrip()
                if comment is None or comment.rstrip() == line_cmnt:
                    l = line.split()
                    assert len(l) >= 16
                    new_epoch  = dt.datetime.strptime(l[0]+' '+l[1], '%Y-%m-%d %H:%M:%S')
                    new_record = [x for x in l[2:14]]
                    new_record.append(line_cmnt)
                    t_stamp    = dt.datetime.strptime(l[14]+' '+l[15], '%Y-%m-%d %H:%M:%S.%f')
                    sorted_appender(epochs_array, rest_array, time_stamps, new_epoch, new_record, t_stamp) 
                else:
                    print '## [DEBUG] Skipping record with description: \"{:}\"'.format(line_cmnt)
        assert len(epochs_array) == len(rest_array) and len(epochs_array) == len(time_stamps)
        print 'Time-Series from {:} to {:}'.format(epochs_array[0], epochs_array[len(epochs_array)-1])
        return ts.TimeSeries(name=cts_file[0:4], type=ts.CoordinateType.Cartesian,
                    # following won't work. Why ?? Use classic datetime instead
                    # epoch_array=np.array([x.strftime('%Y-%m-%d:%H:%M:%S') for x in sorted_lst[0]], dtype='datetime64'))
                    # epoch_array = sorted_lst[0], ## the epochs list (NOT ARRAY)
                    epoch_array = epochs_array,
                    x_array     = np.array([i[0]  for i in rest_array]), # rest_array[0]),
                    y_array     = np.array([i[2]  for i in rest_array]),
                    z_array     = np.array([i[4]  for i in rest_array]),
                    sx_array    = np.array([i[7]  for i in rest_array]),
                    sy_array    = np.array([i[9]  for i in rest_array]),
                    sz_array    = np.array([i[11] for i in rest_array]),
                    comment     = list(zip(*rest_array)[12]),
                    time_stamp  = time_stamps)

def read_ntua_cts(cts_file):
    ##  format of .c.cts files
    ##  year(0) month(1) day(2) hour(3) doy(4) fractional_year(5) x(6) y(7) z(8) sx(9) sy(10) sz(11)
    with open(cts_file, 'r' ) as fin :
        epochs_array = [] ## a list of datetime instances
        rest_array   = [] ## a list of lists of type [x, y, z, sx, sy, sz]
        for line in fin.readlines():
            if not line[0] == '#':
                l = line.split()
                assert len(l) == 12 ## line must have 12 fields/columns
                dt_str = l[0] + '-' + l[1] + '-' + l[2] + ':' + str(int(float(l[3])))
                epochs_array.append(dt.datetime.strptime(dt_str,'%Y-%m-%d:%H'))
                rest_array.append( [x for x in l[6:]] )

        ##  sort the elements of the two lists based in date (i.e. epochs_array)
        ##  this will create a new list with len()=2, where the first list is the
        ##+ sorted epochs_array and the second is the corresponding rest_array
        sorted_lst = [list(x) for x in zip(*sorted(zip(epochs_array, rest_array), key=lambda pair: pair[0]))]
        ##  Now, sorted_lst[1] contains lists of type: [x,y,z,sx,sy,sz]. Transpose that
        ##  so that we can easily read np.array(s) (i.e. make sorted_lst[1]=[[x1,x2,...],[y1,y2,...]...]]
        rest_array = map(list, zip(*sorted_lst[1]))

        ## Return a TimeSeries object
        return ts.TimeSeries(name=cts_file[0:4], type=ts.CoordinateType.Cartesian,
                    # following won't work. Why ?? Use classic datetime instead
                    # epoch_array=np.array([x.strftime('%Y-%m-%d:%H:%M:%S') for x in sorted_lst[0]], dtype='datetime64'))
                    epoch_array = sorted_lst[0], ## the epochs list (NOT ARRAY)
                    x_array     = np.array(rest_array[0]),
                    y_array     = np.array(rest_array[1]),
                    z_array     = np.array(rest_array[2]),
                    sx_array    = np.array(rest_array[3]),
                    sy_array    = np.array(rest_array[4]),
                    sz_array    = np.array(rest_array[5]))
