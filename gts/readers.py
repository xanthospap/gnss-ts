import numpy      as np
import datetime   as dt
import timeseries as ts

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
